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
import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np

from car_pipeline.data.source import (
    CacheEntry, CacheError, DataSource, stream_to_file)

#: No dataset-specific constants live here any more. Which series, which
#: columns, which category values and which compartment map all belong to the
#: atlas being read, and are declared on an AtlasSchema in the indication
#: config. A module global naming one submission is exactly what made this
#: loader single-indication.

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


# The authors' top-level branches, collapsed only where two of them describe the
# same compartment. Everything not named here is reported as other rather than
# forced into one of the five.
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
    name = "Single-cell tumour atlas"
    namespace = "singlecell"

    def __init__(self, atlas=None, root=None) -> None:
        """`atlas` is an AtlasSchema; None keeps the reference submission.

        Every reference to a column name, a category value or a filename now
        goes through this object, so a second atlas is a declaration rather than
        an edit to library code.
        """
        super().__init__(root=root)
        if atlas is None:
            from car_pipeline.configs.pdac import PDAC_ATLAS
            atlas = PDAC_ATLAS
        self.atlas = atlas
        #: The registry key stays the class name; the accession is detail,
        #: reported per run rather than baked into the dataset identity.
        self.series_name = f"GEO {atlas.series}"

    def cache_entries(self) -> Iterable[CacheEntry]:
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
                    # Bumped when the derivation changes in a way that alters
                    # what is stored. Version 2 decodes the identifier column
                    # through its category table; version 1 stored the raw
                    # codes, which no lookup could ever match.
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
        entry = self._entry("archive")

        def fetcher(tmp: Path) -> dict:
            print(f"  fetching {self.atlas.archive}", flush=True)
            return stream_to_file(self.atlas.url, tmp,
                                  progress_label="archive", timeout=1800)

        return self.cache.ensure(entry, fetcher)

    #: Set where recovering the matrix is not an option — a container that ships
    #: only the derived summaries. Fetching the archive there would download
    #: 2.6 GB and expand it to 8.3 GB onto an in-memory filesystem, so the
    #: instance is killed mid-job instead of reporting anything.
    OFFLINE_ENV = "CART_NO_MATRIX_FETCH"

    def matrix_path(self) -> Path:
        entry = self._entry("matrix")
        if self.cache.is_valid(entry):
            return self.cache.path(entry)

        # The offline guard comes first. It used to sit below the non-gzip
        # branch, so a deployment that forbids fetching the matrix would happily
        # download an uncompressed one instead -- 844 MB on an in-memory
        # filesystem, which is the exact outcome the variable exists to stop.
        if os.environ.get(self.OFFLINE_ENV):
            raise CacheError(
                "the single-cell matrix is absent and this deployment cannot "
                "fetch it. The derived summaries under data/singlecell cover "
                "one gene pool; a different indication changes the pool digest "
                "and needs the matrix, which must be materialised ahead of "
                "time rather than downloaded here. See specs/deployment.md."
            )

        # Not every atlas arrives compressed. The reference submission is a
        # gzipped h5ad that has to be expanded; a CELLxGENE export is already an
        # h5ad, and running it through the expander would fail on the first read
        # rather than having nothing to do.
        if not self.atlas.archive.endswith(".gz"):
            return self.archive_path()

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

    def _read_field(self, group: "h5py.Group", field: str) -> np.ndarray | None:
        """One var/obs field, whether it is the index or a named column.

        The two submissions are mirror images: the reference indexes by gene
        symbol and carries Ensembl in a column, the CELLxGENE export indexes by
        Ensembl and carries the symbol in `feature_name`. Resolving both here
        means the readers below never have to know which is which.
        """
        if field == "_index":
            return self._read_index(group)
        return self._read_column(group, field)

    def _counts_matrix(self, fh: "h5py.File"):
        """The raw-counts matrix, wherever this submission put it.

        Reading the wrong matrix is the failure that does not announce itself:
        normalised values would be consumed as though they were counts and every
        downstream count would be wrong but plausible. So a missing path raises
        by name rather than falling back to X.
        """
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
        key = group.attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        return cls._decode(group[key][:])

    # -- aggregation ------------------------------------------------------

    def build_group_means(self) -> Path:
        entry = self._entry("group_means")

        def fetcher(tmp: Path) -> dict:
            a = self.atlas
            with h5py.File(self.matrix_path(), "r") as fh:
                obs = fh["obs"]
                genes = self._read_field(fh["var"], a.symbol_field)
                ensembl = self._read_field(fh["var"], a.ensembl_field)
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

                a = self.atlas
                l1_codes, l1_cats = self._read_categorical(obs, a.level1_column)
                l3_codes, l3_cats = self._read_categorical(obs, a.level3_column)
                # An atlas without a treatment split is read without one
                # rather than refused: the untreated subset is then absent,
                # which is a third state, not a zero.
                tr_codes, tr_cats = (
                    self._read_categorical(obs, a.treatment_column)
                    if a.treatment_column else (None, None))

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

                # An atlas with no treatment split, or one whose split does not
                # carry the declared label, yields an empty untreated subset
                # rather than raising. The group table then holds the "all" rows
                # with the untreated rows present and empty, which is what a
                # reader must see: not measured, rather than measured as zero.
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

    # -- per cell, malignant compartment ----------------------------------

    def malignant_entry(self, genes: list[str]) -> CacheEntry:
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
                # The gene set is part of what this artifact *is*. Without it a
                # changed pool would silently reuse the wrong columns.
                "genes": digest,
                "n_genes": len(genes),
                "derived_version": 1,
            },
        )

    def build_malignant(self, genes: list[str]) -> Path:
        """Stream raw counts for one gene set over malignant cells only.

        Kept as its own cache entry rather than folded into the group means.
        That artifact consumes the cell axis inside its accumulation loop and
        stores 78 x 22,164 group means, so no conjunction over cells can be
        recovered from it; and leaving it untouched means the ranking stage is
        not invalidated by anything done here.
        """
        entry = self.malignant_entry(genes)

        def fetcher(tmp: Path) -> dict:
            with h5py.File(self.matrix_path(), "r") as fh:
                var_names = list(
                    self._read_field(fh["var"], self.atlas.symbol_field))
                column_of = {g: i for i, g in enumerate(var_names)}
                present = [g for g in genes if g in column_of]
                missing = [g for g in genes if g not in column_of]
                cols = np.asarray([column_of[g] for g in present], dtype=np.int64)

                # A full width lookup rather than a search over each row. The
                # column indices in this file are stored in descending order
                # within a row, so anything relying on the ascending order this
                # format usually carries matches nothing and reports every gene
                # as absent, silently.
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

                # The donor identifier, whatever this submission calls it.
                # Reading a hardcoded name gave a bare TypeError three frames
                # later when the column was absent, which named nothing.
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
                # Guarded the same way the group-means reader is. It was not,
                # so an atlas whose treatment column is stored plainly, or whose
                # declared label is absent, built its group means fine and then
                # died here -- the same shape as the donor-column bug this
                # replaced.
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

                # Written through a handle: passing the path would have numpy
                # append its own suffix, and the cache would then commit a file
                # that is not the one it just wrote.
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


#: The layer holding raw integer counts. `X` is log1p(CP10K) and is normalised
#: per cell, so a fixed threshold on it is not a fixed count threshold: raw depth
#: across malignant cells runs from 96 to 9,642. Detection is defined on counts,
#: so counts is what gets read.


JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"


@dataclass
class MalignantCells:
    """Per-cell counts for a fixed gene set, malignant compartment only."""

    genes: list[str]              # requested order, one column each
    counts: np.ndarray            # cells x genes, uint16
    patient: np.ndarray           # cells, patient label
    untreated: np.ndarray         # cells, bool
    depth: np.ndarray             # cells, total counts across all genes
    missing: list[str]            # requested genes with no column in the matrix

    def positive(self, threshold: int = 1) -> np.ndarray:
        return self.counts >= threshold

    def evaluable_patients(self, minimum: int = 100) -> list[str]:
        labels, counts = np.unique(self.patient, return_counts=True)
        return sorted(str(l) for l, c in zip(labels, counts) if c >= minimum)

    def subset(self, genes: list[str]) -> "MalignantCells":
        """Narrow to a subset of the columns already held.

        A smaller gene set is a slice of this one, not a different derivation,
        so taking it here avoids streaming 8.3 GB again to answer a question the
        loaded matrix already contains.
        """
        column = {g: i for i, g in enumerate(self.genes)}
        # Genes with no column are recorded rather than rejected, exactly as the
        # builder records them. Some requested genes have no row in the matrix
        # at all, so a subset that refused them would behave differently from
        # the derivation it is a subset of.
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
    canonical = "\n".join(genes)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


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
