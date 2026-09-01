"""Per cell type expression for the indication."""

from __future__ import annotations

import gzip
import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np

from car_pipeline.data.source import (
    CacheEntry, CacheError, DataSource, stream_to_file)


DROPOUT_EPSILON = 0.001

IMMUNE = "immune"
MALIGNANT = "malignant"
FIBROBLAST = "fibroblast"
EPITHELIAL = "epithelial non-malignant"
ENDOTHELIAL = "endothelial"
OTHER = "other"

SUBSET_ALL = "all"
SUBSET_UNTREATED = "untreated"


COMPARTMENT_ORDER = [MALIGNANT, FIBROBLAST, EPITHELIAL, IMMUNE, ENDOTHELIAL, OTHER]


STROMAL_IMMUNE = [FIBROBLAST, IMMUNE, ENDOTHELIAL]

ROW_BLOCK = 8192


@dataclass
class Atlas:
    genes: np.ndarray
    groups: list[tuple[str, str]]
    group_means: np.ndarray
    compartments: list[str]
    compartment_means: np.ndarray
    compartment_counts: dict[str, int]
    level3_to_level1: dict[str, str]
    ensembl: np.ndarray | None = None
    per_cell_total: np.ndarray | None = None

    def gene_index(self) -> dict[str, int]:
        """Column positions for the genes this matrix carries."""
        out: dict[str, int] = {}
        for i, g in enumerate(self.genes):
            out.setdefault(str(g), i)
        return out

    def compartment_value(self, compartment: str, gene_index: int) -> float:
        """One gene's mean in one compartment."""
        return float(self.compartment_means[self.compartments.index(compartment)][gene_index])

    def peak_group(self, gene_index: int) -> float:
        """Highest group mean, ignoring groups that hold no cells."""
        column = self.group_means[:, gene_index]
        if np.all(np.isnan(column)):
            return float("nan")
        return float(np.nanmax(column))


class SingleCellSource(DataSource):
    name = "Single-cell tumour atlas"
    namespace = "singlecell"

    def __init__(self, atlas=None, root=None) -> None:
        """`atlas` is an AtlasSchema; None keeps the reference submission."""
        super().__init__(root=root)
        if atlas is None:
            from car_pipeline.configs.pdac import PDAC_ATLAS
            atlas = PDAC_ATLAS
        self.atlas = atlas

        self.series_name = f"GEO {atlas.series}"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The archive, the expanded matrix and the summaries derived from it."""
        a = self.atlas
        tag = a.slug
        fp = {"series": a.series, "file": a.archive}
        return [
            CacheEntry(key=f"archive__{tag}", filename=a.archive, fingerprint=fp),
            CacheEntry(
                key=f"matrix__{tag}",
                filename=f"matrix__{tag}.h5ad",
                fingerprint=fp,
            ),
            CacheEntry(
                key=f"group_means__{tag}",
                filename=f"group_means__{tag}.npz",
                fingerprint={
                    **fp,
                    "epsilon": DROPOUT_EPSILON,
                    "derived_version": 2,
                },
            ),
        ]

    def _entry(self, kind: str) -> CacheEntry:
        """Look an entry up by unqualified kind; the key carries the accession."""
        want = f"{kind}__{self.atlas.slug}"
        for entry in self.cache_entries():
            if entry.key == want:
                return entry
        raise KeyError(f"no cache entry {want!r} for atlas {self.atlas.series!r}")

    def archive_path(self) -> Path:
        """The downloaded archive, fetching it if absent."""
        entry = self._entry("archive")

        def fetcher(tmp: Path) -> dict:
            """Fetch or derive one artifact into a temporary file."""
            print(f"  fetching {self.atlas.archive}", flush=True)
            return stream_to_file(self.atlas.url, tmp,
                                  progress_label="archive", timeout=1800)

        return self.cache.ensure(entry, fetcher)

    OFFLINE_ENV = "CART_NO_MATRIX_FETCH"

    def matrix_path(self) -> Path:
        """The expanded matrix, refusing to fetch it where that is forbidden."""
        entry = self._entry("matrix")
        if self.cache.is_valid(entry):
            return self.cache.path(entry)

        if os.environ.get(self.OFFLINE_ENV):
            raise CacheError(
                "the single-cell matrix is absent and this deployment cannot "
                "fetch it. The derived summaries under data/singlecell cover "
                "one gene pool; a different indication changes the pool digest "
                "and needs the matrix, which must be materialised ahead of "
                "time rather than downloaded here. See specs/deployment.md."
            )

        if not self.atlas.archive.endswith(".gz"):
            return self.archive_path()

        source = self.archive_path()

        def fetcher(tmp: Path) -> dict:
            """Fetch or derive one artifact into a temporary file."""
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

    @staticmethod
    def _decode(values) -> np.ndarray:
        """Decode a stored string field to text."""
        return np.asarray(
            [v.decode() if isinstance(v, bytes) else str(v) for v in values]
        )

    @classmethod
    def _read_categorical(
        cls, obs: h5py.Group, name: str
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return integer codes and their category labels."""
        node = obs[name]
        if isinstance(node, h5py.Group) and "categories" in node:
            return node["codes"][:], cls._decode(node["categories"][:])
        if "__categories" in obs and name in obs["__categories"]:
            return node[:], cls._decode(obs["__categories"][name][:])

        return node[:], None

    @classmethod
    def _read_column(cls, group: h5py.Group, name: str) -> np.ndarray | None:
        """Decode a column that may be stored as codes into a shared table."""
        if name not in group:
            return None
        codes, cats = cls._read_categorical(group, name)
        if cats is None:
            return cls._decode(codes)

        return np.asarray(
            [cats[c] if c >= 0 else "" for c in codes], dtype=str
        )

    def _read_field(self, group: "h5py.Group", field: str) -> np.ndarray | None:
        """One var/obs field, whether it is the index or a named column."""
        if field == "_index":
            return self._read_index(group)
        return self._read_column(group, field)

    def _counts_matrix(self, fh: "h5py.File"):
        """The raw-counts matrix, wherever this submission put it."""
        node = fh
        for part in self.atlas.counts_path.split("/"):
            if part not in node:
                raise KeyError(
                    f"{self.atlas.series}: no raw counts at "
                    f"{self.atlas.counts_path!r} (looking for {part!r}); "
                    f"available here: {sorted(node.keys())}"
                )
            node = node[part]
        return node

    @classmethod
    def _read_index(cls, group: h5py.Group) -> np.ndarray:
        """Read an index column from the stored annotation."""
        key = group.attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        return cls._decode(group[key][:])

    def build_group_means(self) -> Path:
        """Accumulate per-group means by streaming the matrix in row blocks."""
        entry = self._entry("group_means")

        def fetcher(tmp: Path) -> dict:
            """Fetch or derive one artifact into a temporary file."""
            a = self.atlas
            with h5py.File(self.matrix_path(), "r") as fh:
                obs = fh["obs"]
                genes = self._read_field(fh["var"], a.symbol_field)
                ensembl = self._read_field(fh["var"], a.ensembl_field)
                if ensembl is None:
                    ensembl = np.asarray([""] * len(genes))

                recognisable = sum(1 for e in ensembl if str(e).startswith("ENSG"))
                if recognisable < 0.5 * len(ensembl):
                    raise ValueError(
                        "the identifier column did not decode to identifiers "
                        f"({recognisable:,} of {len(ensembl):,} recognisable); "
                        f"first few read as {list(ensembl[:3])}"
                    )

                a = self.atlas
                l1_codes, l1_cats = self._read_categorical(obs, a.level1_column)
                l3_codes, l3_cats = self._read_categorical(obs, a.level3_column)

                tr_codes, tr_cats = (
                    self._read_categorical(obs, a.treatment_column)
                    if a.treatment_column else (None, None))

                n_cells = l3_codes.shape[0]
                n_genes = len(genes)
                n_types = len(l3_cats)

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
                    compartment = self.atlas.compartment_map.get(parent, OTHER)
                    level3_to_level1[str(t_name)] = compartment
                    print(f"    {t_name:38s} {parent:30s} -> {compartment}")

                compartment_of_l1 = np.asarray(
                    [
                        COMPARTMENT_ORDER.index(
                            self.atlas.compartment_map.get(str(c), OTHER))
                        for c in l1_cats
                    ]
                )
                comp_of_cell = compartment_of_l1[l1_codes]

                if tr_cats is not None and a.untreated_label in list(tr_cats):
                    untreated_idx = list(tr_cats).index(a.untreated_label)
                    untreated = tr_codes == untreated_idx
                else:
                    if a.treatment_column:
                        print(f"    no {a.untreated_label!r} category in "
                              f"{a.treatment_column!r}; untreated subset empty")
                    untreated = np.zeros(l1_codes.shape[0], dtype=bool)

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
                    """Mean per group, or not-a-number where the group is empty."""
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
        """The cached group means and the axes they are indexed by."""
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

    def malignant_entry(self, genes: list[str]) -> CacheEntry:
        """The cache entry for one gene set's malignant-cell matrix."""
        digest = _gene_digest(genes)
        tag = self.atlas.slug
        return CacheEntry(
            key=f"malignant_cells__{tag}_{digest}",
            filename=f"malignant_cells__{tag}_{digest}.npz",
            fingerprint={
                "series": self.atlas.series,
                "file": self.atlas.archive,
                "layer": self.atlas.counts_path,
                "compartment": self.atlas.malignant_label,
                "genes": digest,
                "n_genes": len(genes),
                "derived_version": 1,
            },
        )

    def build_malignant(self, genes: list[str]) -> Path:
        """Stream raw counts for one gene set over malignant cells only."""
        entry = self.malignant_entry(genes)

        def fetcher(tmp: Path) -> dict:
            """Fetch or derive one artifact into a temporary file."""
            with h5py.File(self.matrix_path(), "r") as fh:
                var_names = list(
                    self._read_field(fh["var"], self.atlas.symbol_field))
                column_of = {g: i for i, g in enumerate(var_names)}
                present = [g for g in genes if g in column_of]
                missing = [g for g in genes if g not in column_of]
                cols = np.asarray([column_of[g] for g in present], dtype=np.int64)

                lut = np.full(len(var_names), -1, dtype=np.int32)
                lut[cols] = np.arange(len(cols), dtype=np.int32)

                obs = fh["obs"]
                a = self.atlas
                l1_codes, l1_cats = self._read_categorical(obs, a.level1_column)
                if l1_cats is None or a.malignant_label not in list(l1_cats):
                    raise KeyError(
                        f"the annotation has no {a.malignant_label!r} branch; "
                        f"found {sorted(set(map(str, l1_cats or [])))}"
                    )
                malignant = l1_codes == list(l1_cats).index(a.malignant_label)
                n_cells = int(malignant.sum())
                if n_cells == 0:
                    raise ValueError("no malignant cells selected")

                if not a.patient_column:
                    raise CacheError(
                        f"{a.series}: no patient column declared, so per-patient "
                        "prevalence cannot be measured for this atlas"
                    )
                pid = self._read_field(obs, a.patient_column)
                if pid is None:
                    raise CacheError(
                        f"{a.series}: obs has no {a.patient_column!r} column; "
                        f"available: {sorted(obs.keys())[:12]}..."
                    )

                if a.treatment_column:
                    tr_codes, tr_cats = self._read_categorical(obs, a.treatment_column)
                    if tr_cats is not None and a.untreated_label in list(tr_cats):
                        untreated_all = tr_codes == list(tr_cats).index(
                            a.untreated_label)
                    else:
                        untreated_all = np.zeros(n_cells_total, dtype=bool)
                else:
                    untreated_all = np.zeros(l1_codes.shape[0], dtype=bool)

                layer = self._counts_matrix(fh)
                data, indices = layer["data"], layer["indices"]
                iptr = layer["indptr"][:]

                counts = np.zeros((n_cells, len(present)), dtype=np.uint16)
                depth = np.zeros(n_cells, dtype=np.int64)
                written = 0
                for start in range(0, len(iptr) - 1, ROW_BLOCK):
                    stop = min(start + ROW_BLOCK, len(iptr) - 1)
                    keep_rows = np.nonzero(malignant[start:stop])[0]
                    if keep_rows.size == 0:
                        continue
                    lo, hi = int(iptr[start]), int(iptr[stop])
                    blk_i = indices[lo:hi]
                    blk_d = data[lo:hi]
                    off = iptr[start:stop + 1] - lo
                    rows_here = stop - start
                    row_of = np.repeat(np.arange(rows_here), np.diff(off))

                    out_row = np.full(rows_here, -1, dtype=np.int64)
                    out_row[keep_rows] = np.arange(
                        written, written + keep_rows.size
                    )
                    target = out_row[row_of]
                    wanted = target >= 0

                    depth += np.bincount(
                        target[wanted], weights=blk_d[wanted], minlength=n_cells
                    ).astype(np.int64)

                    take = wanted & (lut[blk_i] >= 0)
                    if take.any():
                        vals = blk_d[take]
                        if vals.max() > np.iinfo(np.uint16).max:
                            raise ValueError("a count exceeds the stored width")
                        counts[target[take], lut[blk_i][take]] = vals
                    written += keep_rows.size

                if written != n_cells:
                    raise ValueError(
                        f"wrote {written} rows for {n_cells} malignant cells"
                    )

                with open(tmp, "wb") as out:
                    np.savez_compressed(
                        out,
                        genes=np.asarray(present, dtype=str),
                        missing=np.asarray(missing, dtype=str),
                        counts=counts,
                        patient=np.asarray(
                            [str(p) for p in pid[malignant]], dtype=str
                        ),
                        untreated=untreated_all[malignant],
                        depth=depth,
                    )

            return {
                "observed_rows": n_cells,
                "extra": {
                    "cells": n_cells,
                    "genes": len(present),
                    "missing": len(missing),
                },
            }

        return self.cache.ensure(entry, fetcher)

    def load_malignant(self, genes: list[str]) -> MalignantCells:
        """Per-cell counts over malignant cells for one gene set."""
        path = self.build_malignant(genes)
        with np.load(path, allow_pickle=False) as d:
            return MalignantCells(
                genes=[str(g) for g in d["genes"]],
                counts=d["counts"],
                patient=d["patient"],
                untreated=d["untreated"],
                depth=d["depth"],
                missing=[str(g) for g in d["missing"]],
            )


JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"


@dataclass
class MalignantCells:
    """Per-cell counts for a fixed gene set, malignant compartment only."""

    genes: list[str]
    counts: np.ndarray
    patient: np.ndarray
    untreated: np.ndarray
    depth: np.ndarray
    missing: list[str]

    def positive(self, threshold: int = 1) -> np.ndarray:
        """Cells at or above the detection threshold for this gene."""
        return self.counts >= threshold

    def evaluable_patients(self, minimum: int = 100) -> list[str]:
        """Patients carrying enough cells to be counted."""
        labels, counts = np.unique(self.patient, return_counts=True)
        return sorted(str(l) for l, c in zip(labels, counts) if c >= minimum)

    def subset(self, genes: list[str]) -> "MalignantCells":
        """Narrow to a subset of the columns already held."""
        column = {g: i for i, g in enumerate(self.genes)}

        present = [g for g in genes if g in column]
        absent = [g for g in genes if g not in column]
        return MalignantCells(
            genes=present,
            counts=self.counts[:, [column[g] for g in present]],
            patient=self.patient,
            untreated=self.untreated,
            depth=self.depth,
            missing=sorted(set(self.missing) | set(absent)),
        )


def _gene_digest(genes: Iterable[str]) -> str:
    """A stable digest of the gene set an artifact was built for."""
    canonical = "\n".join(genes)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def match_surface(
    atlas: Atlas, surface, atlas_by_accession: dict
) -> dict[str, tuple[int, str]]:
    """Symbol-first join with an identifier bridge for renamed genes."""
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
        """One gene's malignant mean against its peak stromal compartment."""
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
