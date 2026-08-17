"""Per cell type expression for the indication.

The obvious atlas has no tumour dataset for this indication — every pancreatic
collection in it is normal, islet or diabetes — so a tumour series is used
instead, and the substitution is recorded in the stage 1 dataset list rather
than left implicit.

The archive is far larger than memory both compressed and expanded, so it is
streamed on download, streamed again on expansion, and read back in row blocks.
Nothing here ever holds the whole matrix.

The authors' own annotations are used as given. The single editorial act is
merging their two immune branches into one compartment, and every mapping is
printed so that act is visible rather than buried.

**A zero here never rejects a target.** This assay drops transcripts that bulk
measurement finds abundantly present; a gene reading zero across every cell type
is evidence about the assay, not about the protein. This source separates
compartments. It does not refute.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np

from car_pipeline.data.source import CacheEntry, DataSource, stream_to_file

SERIES = "GSE202051"
ARCHIVE = f"{SERIES}_totaldata-final-toshare.h5ad.gz"
URL = (
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/{SERIES}/suppl/{ARCHIVE}"
)

# Below this, a group mean is treated as a capture failure rather than absence.
# An exact-zero test misses the case entirely: genes that bulk measurement puts
# in the hundreds of transcripts read here at a ten-thousandth, not at zero.
DROPOUT_EPSILON = 0.001

IMMUNE = "immune"
MALIGNANT = "malignant"
FIBROBLAST = "fibroblast"
EPITHELIAL = "epithelial non-malignant"
ENDOTHELIAL = "endothelial"
OTHER = "other"

SUBSET_ALL = "all"
SUBSET_UNTREATED = "untreated"
UNTREATED_LABEL = "Untreated"

LEVEL1 = "Level 1 Annotation"
LEVEL3 = "Level 3 Annotation"
TREATMENT = "treatment_status"

# The authors' top-level branches, collapsed only where two of them describe the
# same compartment. Everything not named here is reported as other rather than
# forced into one of the five.
COMPARTMENT_MAP = {
    "Epithelial (malignant)": MALIGNANT,
    "Cancer-associated fibroblast": FIBROBLAST,
    "Epithelial (non-malignant)": EPITHELIAL,
    "Lymphoid": IMMUNE,
    "Myeloid": IMMUNE,
    "Endothelial": ENDOTHELIAL,
}
COMPARTMENT_ORDER = [MALIGNANT, FIBROBLAST, EPITHELIAL, IMMUNE, ENDOTHELIAL, OTHER]

# Compartments a genuine tumour antigen must rise above. Peak is taken over
# these rather than averaged: one compartment expressing it is enough to matter.
STROMAL_IMMUNE = [FIBROBLAST, IMMUNE, ENDOTHELIAL]

ROW_BLOCK = 8192


@dataclass
class Atlas:
    genes: np.ndarray
    groups: list[tuple[str, str]]          # (subset, cell type)
    group_means: np.ndarray                # groups x genes, mean of expm1
    compartments: list[str]
    compartment_means: np.ndarray          # compartments x genes
    compartment_counts: dict[str, int]
    level3_to_level1: dict[str, str]
    ensembl: np.ndarray | None = None
    per_cell_total: np.ndarray | None = None

    def gene_index(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for i, g in enumerate(self.genes):
            out.setdefault(str(g), i)
        return out

    def compartment_value(self, compartment: str, gene_index: int) -> float:
        return float(self.compartment_means[self.compartments.index(compartment)][gene_index])

    def peak_group(self, gene_index: int) -> float:
        """Highest group mean, ignoring groups that hold no cells."""
        column = self.group_means[:, gene_index]
        if np.all(np.isnan(column)):
            return float("nan")
        return float(np.nanmax(column))


class SingleCellSource(DataSource):
    name = f"GEO {SERIES}"
    namespace = "singlecell"

    def cache_entries(self) -> Iterable[CacheEntry]:
        fp = {"series": SERIES, "file": ARCHIVE}
        return [
            CacheEntry(key="archive", filename=ARCHIVE, fingerprint=fp),
            CacheEntry(
                key="matrix",
                filename="totaldata-final-toshare.h5ad",
                fingerprint=fp,
            ),
            CacheEntry(
                key="group_means",
                filename="group_means.npz",
                fingerprint={
                    **fp,
                    "epsilon": DROPOUT_EPSILON,
                    # Bumped when the derivation changes in a way that alters
                    # what is stored. Version 2 decodes the identifier column
                    # through its category table; version 1 stored the raw
                    # codes, which no lookup could ever match.
                    "derived_version": 2,
                },
            ),
        ]

    def _entry(self, key: str) -> CacheEntry:
        return next(e for e in self.cache_entries() if e.key == key)

    def archive_path(self) -> Path:
        entry = self._entry("archive")

        def fetcher(tmp: Path) -> dict:
            print(f"  fetching {ARCHIVE}", flush=True)
            return stream_to_file(URL, tmp, progress_label="archive", timeout=1800)

        return self.cache.ensure(entry, fetcher)

    def matrix_path(self) -> Path:
        entry = self._entry("matrix")
        if self.cache.is_valid(entry):
            return self.cache.path(entry)

        source = self.archive_path()

        def fetcher(tmp: Path) -> dict:
            print("  expanding archive", flush=True)
            import hashlib
            import os

            digest = hashlib.sha256()
            total = 0
            with gzip.open(source, "rb") as fin, open(tmp, "wb") as fout:
                while True:
                    block = fin.read(1 << 22)
                    if not block:
                        break
                    fout.write(block)
                    digest.update(block)
                    total += len(block)
                    if total % (1 << 30) < (1 << 22):
                        print(f"    expanded {total / 1_073_741_824:.1f} GB", flush=True)
                fout.flush()
                os.fsync(fout.fileno())
            print(f"    expanded {total / 1_073_741_824:.2f} GB", flush=True)
            return {"digest": digest.hexdigest(), "extra": {"expanded_bytes": total}}

        return self.cache.ensure(entry, fetcher)

    # -- annotation -------------------------------------------------------

    @staticmethod
    def _decode(values) -> np.ndarray:
        return np.asarray(
            [v.decode() if isinstance(v, bytes) else str(v) for v in values]
        )

    @classmethod
    def _read_categorical(
        cls, obs: h5py.Group, name: str
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return integer codes and their category labels.

        Two encodings are in circulation for this: one stores categories beside
        the codes, the older one keeps them in a shared table.
        """
        node = obs[name]
        if isinstance(node, h5py.Group) and "categories" in node:
            return node["codes"][:], cls._decode(node["categories"][:])
        if "__categories" in obs and name in obs["__categories"]:
            return node[:], cls._decode(obs["__categories"][name][:])
        # Stored plainly, with no table behind it.
        return node[:], None

    @classmethod
    def _read_column(cls, group: h5py.Group, name: str) -> np.ndarray | None:
        """Decode a column that may be stored as codes into a shared table.

        Reading such a column raw yields the integer codes rendered as text,
        which look like perfectly good values and match nothing. That failure is
        entirely silent: every lookup against them simply misses, and a missed
        lookup is indistinguishable from a value the source never carried.
        """
        if name not in group:
            return None
        codes, cats = cls._read_categorical(group, name)
        if cats is None:
            return cls._decode(codes)
        # Built as text rather than as objects: an object array cannot be
        # reloaded from the cache without allowing arbitrary deserialisation,
        # which is not a thing a data cache should ever need.
        return np.asarray(
            [cats[c] if c >= 0 else "" for c in codes], dtype=str
        )

    @classmethod
    def _read_index(cls, group: h5py.Group) -> np.ndarray:
        key = group.attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        return cls._decode(group[key][:])

    # -- aggregation ------------------------------------------------------

    def build_group_means(self) -> Path:
        entry = self._entry("group_means")

        def fetcher(tmp: Path) -> dict:
            with h5py.File(self.matrix_path(), "r") as fh:
                obs = fh["obs"]
                genes = self._read_index(fh["var"])
                ensembl = self._read_column(fh["var"], "ensg")
                if ensembl is None:
                    ensembl = np.asarray([""] * len(genes))
                # A bridge built from values of the wrong kind matches nothing
                # and reports every protein as absent. Checked here so the
                # failure is loud at load rather than invisible at join.
                recognisable = sum(1 for e in ensembl if str(e).startswith("ENSG"))
                if recognisable < 0.5 * len(ensembl):
                    raise ValueError(
                        "the identifier column did not decode to identifiers "
                        f"({recognisable:,} of {len(ensembl):,} recognisable); "
                        f"first few read as {list(ensembl[:3])}"
                    )

                l1_codes, l1_cats = self._read_categorical(obs, LEVEL1)
                l3_codes, l3_cats = self._read_categorical(obs, LEVEL3)
                tr_codes, tr_cats = self._read_categorical(obs, TREATMENT)

                n_cells = l3_codes.shape[0]
                n_genes = len(genes)
                n_types = len(l3_cats)

                # The finer axis must sit entirely inside the coarser one. If a
                # cell type straddled two branches, every compartment figure
                # below would be a blend of two things.
                nesting: dict[str, set[str]] = {c: set() for c in l3_cats}
                for t_idx, t_name in enumerate(l3_cats):
                    parents = np.unique(l1_codes[l3_codes == t_idx])
                    nesting[t_name] = {l1_cats[p] for p in parents}
                straddling = {k: v for k, v in nesting.items() if len(v) > 1}
                if straddling:
                    raise ValueError(
                        "cell types spanning more than one branch: "
                        + "; ".join(f"{k} -> {sorted(v)}" for k, v in straddling.items())
                    )

                print("  cell type to compartment")
                level3_to_level1: dict[str, str] = {}
                for t_name in l3_cats:
                    parent = next(iter(nesting[t_name])) if nesting[t_name] else ""
                    compartment = COMPARTMENT_MAP.get(parent, OTHER)
                    level3_to_level1[str(t_name)] = compartment
                    print(f"    {t_name:38s} {parent:30s} -> {compartment}")

                compartment_of_l1 = np.asarray(
                    [
                        COMPARTMENT_ORDER.index(COMPARTMENT_MAP.get(str(c), OTHER))
                        for c in l1_cats
                    ]
                )
                comp_of_cell = compartment_of_l1[l1_codes]

                untreated_idx = int(np.where(tr_cats == UNTREATED_LABEL)[0][0])
                untreated = tr_codes == untreated_idx

                x = fh["X"]
                indptr = x["indptr"][:]
                data_ds = x["data"]
                indices_ds = x["indices"]

                sums_all = np.zeros(n_types * n_genes, dtype=np.float64)
                sums_unt = np.zeros(n_types * n_genes, dtype=np.float64)
                sums_comp = np.zeros(
                    len(COMPARTMENT_ORDER) * n_genes, dtype=np.float64
                )
                per_cell_total = np.zeros(n_cells, dtype=np.float64)

                print(f"  aggregating {n_cells:,} cells in blocks", flush=True)
                for start in range(0, n_cells, ROW_BLOCK):
                    stop = min(start + ROW_BLOCK, n_cells)
                    lo, hi = int(indptr[start]), int(indptr[stop])
                    if hi == lo:
                        continue
                    values = np.expm1(data_ds[lo:hi].astype(np.float64))
                    cols = indices_ds[lo:hi].astype(np.int64)
                    counts = np.diff(indptr[start : stop + 1]).astype(np.int64)
                    rows = np.repeat(np.arange(start, stop, dtype=np.int64), counts)

                    # Per-cell totals, accumulated blockwise. These are the
                    # scale check: the transform is reversible only if they land
                    # where the stated normalisation says they should.
                    per_cell_total[start:stop] += np.bincount(
                        rows - start, weights=values, minlength=stop - start
                    )

                    flat_all = l3_codes[rows].astype(np.int64) * n_genes + cols
                    sums_all += np.bincount(
                        flat_all, weights=values, minlength=n_types * n_genes
                    )

                    flat_comp = comp_of_cell[rows].astype(np.int64) * n_genes + cols
                    sums_comp += np.bincount(
                        flat_comp,
                        weights=values,
                        minlength=len(COMPARTMENT_ORDER) * n_genes,
                    )

                    mask = untreated[rows]
                    if mask.any():
                        sums_unt += np.bincount(
                            flat_all[mask],
                            weights=values[mask],
                            minlength=n_types * n_genes,
                        )

                    if (start // ROW_BLOCK) % 5 == 0:
                        print(f"    {stop:,}/{n_cells:,} cells", flush=True)

                cells_all = np.bincount(l3_codes, minlength=n_types).astype(np.float64)
                cells_unt = np.bincount(
                    l3_codes[untreated], minlength=n_types
                ).astype(np.float64)
                cells_comp = np.bincount(
                    comp_of_cell, minlength=len(COMPARTMENT_ORDER)
                ).astype(np.float64)

                def normalise(sums, counts, rows_n):
                    """Mean per group, or not-a-number where the group is empty.

                    A group with no cells has no mean. Dividing by a substituted
                    one would report it as measured and silent, which is a
                    different and much more reassuring statement than no cells
                    of that type having been captured.
                    """
                    m = sums.reshape(rows_n, n_genes)
                    out = np.full((rows_n, n_genes), np.nan, dtype=np.float32)
                    present = counts > 0
                    out[present] = (
                        m[present] / counts[present][:, None]
                    ).astype(np.float32)
                    return out

                means_all = normalise(sums_all, cells_all, n_types)
                means_unt = normalise(sums_unt, cells_unt, n_types)
                comp_means = normalise(
                    sums_comp, cells_comp, len(COMPARTMENT_ORDER)
                )

                group_labels = np.asarray(
                    [f"{SUBSET_ALL}|{c}" for c in l3_cats]
                    + [f"{SUBSET_UNTREATED}|{c}" for c in l3_cats]
                )
                group_means = np.vstack([means_all, means_unt])

                with open(tmp, "wb") as out:
                    np.savez_compressed(
                        out,
                        genes=genes,
                        ensembl=ensembl,
                        group_labels=group_labels,
                        group_means=group_means,
                        group_cells=np.concatenate([cells_all, cells_unt]),
                        compartments=np.asarray(COMPARTMENT_ORDER),
                        compartment_means=comp_means,
                        compartment_cells=cells_comp,
                        level3=l3_cats,
                        level3_parent=np.asarray(
                            [level3_to_level1[str(c)] for c in l3_cats]
                        ),
                        per_cell_total=per_cell_total.astype(np.float32),
                    )

            return {
                "observed_rows": int(group_means.shape[0]),
                "extra": {
                    "cells": int(n_cells),
                    "genes": int(n_genes),
                    "cell_types": int(n_types),
                },
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> Atlas:
        path = self.build_group_means()
        with np.load(path, allow_pickle=False) as d:
            groups = [tuple(str(g).split("|", 1)) for g in d["group_labels"]]
            atlas = Atlas(
                genes=d["genes"],
                groups=groups,
                group_means=d["group_means"],
                compartments=[str(c) for c in d["compartments"]],
                compartment_means=d["compartment_means"],
                compartment_counts={
                    str(c): int(n)
                    for c, n in zip(d["compartments"], d["compartment_cells"])
                },
                level3_to_level1={
                    str(a): str(b) for a, b in zip(d["level3"], d["level3_parent"])
                },
            )
            atlas.ensembl = d["ensembl"]
            atlas.per_cell_total = d["per_cell_total"]
        return atlas


JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"


def match_surface(
    atlas: Atlas, surface, atlas_by_accession: dict
) -> dict[str, tuple[int, str]]:
    """Symbol-first join with an identifier bridge for renamed genes.

    The route is returned alongside the column so a bridged row can be counted
    as one. An unrecorded route cannot be audited, and this join reaches the two
    heaviest components in the score.
    """
    by_symbol = atlas.gene_index()
    by_ensembl: dict[str, int] = {}
    if atlas.ensembl is not None:
        for i, e in enumerate(atlas.ensembl):
            key = str(e).split(".")[0]
            if key:
                by_ensembl.setdefault(key, i)

    out: dict[str, tuple[int, str]] = {}
    for rec in surface:
        if rec.gene and rec.gene in by_symbol:
            out[rec.accession] = (by_symbol[rec.gene], JOIN_SYMBOL)
            continue
        entry = atlas_by_accession.get(rec.accession)
        if entry is not None and entry.ensembl in by_ensembl:
            out[rec.accession] = (by_ensembl[entry.ensembl], JOIN_ENSEMBL_BRIDGE)
    return out


def describe(path: Path) -> None:
    """Print the layout of the stored matrix and its annotation columns."""
    with h5py.File(path, "r") as fh:
        print("  root keys:", list(fh.keys()))
        x = fh["X"]
        if isinstance(x, h5py.Group):
            print("  X is sparse:", dict(x.attrs))
            for k in x.keys():
                print(f"    X/{k}: shape={x[k].shape} dtype={x[k].dtype}")
        else:
            print(f"  X dense: shape={x.shape} dtype={x.dtype}")
        print("  obs columns:", list(fh["obs"].keys()))
        print("  var columns:", list(fh["var"].keys()))


if __name__ == "__main__":
    from car_pipeline.data.hpa import HPASource, index as atlas_index
    from car_pipeline.data.uniprot import load_surface

    atlas = SingleCellSource().load()
    idx = atlas.gene_index()

    print("\ncompartments")
    expected_comp = {
        MALIGNANT: 64538,
        FIBROBLAST: 54935,
        EPITHELIAL: 38208,
        IMMUNE: 21461,
        ENDOTHELIAL: 17175,
        OTHER: 28671,
    }
    for name, exp in expected_comp.items():
        got = atlas.compartment_counts[name]
        pct = abs(got - exp) / exp * 100
        print(f"  {name:26s} {got:7,}  expected {exp:7,}  ({pct:.2f}% off)")

    print("\nshape")
    for label, got, exp in [
        ("genes", len(atlas.genes), 22164),
        ("groups", len(atlas.groups), 78),
        ("cells", int(sum(atlas.compartment_counts.values())), 224988),
    ]:
        pct = abs(got - exp) / exp * 100
        print(f"  {label}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    if atlas.per_cell_total is not None:
        med = float(np.median(atlas.per_cell_total))
        print(f"  per-cell reversed total (median): {med:,.0f}   expected ~9,600")

    print("\nmarkers, mean per compartment")
    markers = [
        ("KRT19", MALIGNANT, 14.1, FIBROBLAST, 0.02),
        ("COL1A1", FIBROBLAST, 18.3, None, None),
        ("PTPRC", IMMUNE, 16.7, None, None),
        ("VWF", ENDOTHELIAL, 16.2, None, None),
    ]
    for gene, comp, exp, other, exp_other in markers:
        i = idx.get(gene)
        if i is None:
            print(f"  {gene}: absent")
            continue
        got = atlas.compartment_value(comp, i)
        line = f"  {gene:8s} {comp:12s} {got:8.2f}  expected {exp}"
        if other:
            line += (
                f"   {other} {atlas.compartment_value(other, i):.3f}"
                f"  expected {exp_other}"
            )
        print(line)

    print("\ntop genes by peak across cell type groups")
    # Ranked by the largest value any group reaches, not by a population mean.
    # A mean over every cell is dominated by the nuclear transcripts this assay
    # always captures, which say nothing about any cell type in particular.
    peak_by_gene = np.nanmax(atlas.group_means, axis=0)
    top = np.argsort(peak_by_gene)[::-1][:6]
    print("  " + ", ".join(str(atlas.genes[i]) for i in top))
    for enzyme in ("CTRB1", "CPA1", "PRSS1"):
        i = idx.get(enzyme)
        if i is not None:
            rank = int((peak_by_gene > peak_by_gene[i]).sum()) + 1
            print(f"    {enzyme}: rank {rank} of {len(atlas.genes):,}")
    print("  expected the acinar enzymes at the top")

    print("\npurity: malignant against peak stromal or immune")
    def ratio(gene: str) -> float | None:
        i = idx.get(gene)
        if i is None:
            return None
        mal = atlas.compartment_value(MALIGNANT, i)
        peak = max(atlas.compartment_value(c, i) for c in STROMAL_IMMUNE)
        return mal / peak if peak else float("inf")

    for gene, exp in [
        ("HLA-DRA", 0.07), ("CD74", 0.07),
        ("CEACAM6", 193), ("CLDN18", 43), ("MUC1", 35), ("MSLN", 20),
    ]:
        r = ratio(gene)
        print(f"  {gene:9s} {r:9.2f}x   expected {exp}x" if r is not None
              else f"  {gene}: absent")

    print("\ndropout sensitivity across the surface set")
    surface, _ = load_surface()
    by_acc, _ = atlas_index(HPASource().load())
    joined = match_surface(atlas, surface, by_acc)
    peaks = np.asarray([atlas.peak_group(i) for i, _ in joined.values()])
    bridged = sum(1 for _, p in joined.values() if p == JOIN_ENSEMBL_BRIDGE)
    print(f"  joined through the identifier bridge: {bridged}")
    for eps, exp in [(0.0, 267), (DROPOUT_EPSILON, 357), (0.01, 534)]:
        got = int(np.sum(peaks <= eps))
        pct = abs(got - exp) / exp * 100
        print(f"  at {eps:<6}: {got:,} silent everywhere  expected {exp:,}  ({pct:.2f}% off)")

    zero_everywhere = int(np.sum(peaks <= 0.0))
    counts = {
        "surface matched": (len(joined), 3311),
        "zero everywhere": (zero_everywhere, 267),
        "usable": (len(joined) - zero_everywhere, 3044),
        "no row": (len(surface) - len(joined), 185),
    }
    print()
    for label, (got, exp) in counts.items():
        pct = abs(got - exp) / exp * 100
        print(f"  {label}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    i5 = idx.get("CEACAM5")
    if i5 is not None:
        print(
            f"\n  CEACAM5 peak across groups: {atlas.peak_group(i5):.6f}"
            "   bulk puts it near 300 transcripts and 409x normal"
        )
