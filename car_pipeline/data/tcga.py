"""Bulk tumour expression for the indication cohort.

Files are requested in small batches. The whole cohort as a single request is
roughly three quarters of a gigabyte that either arrives complete or fails
complete, with no way to resume; batching turns one all-or-nothing transfer
into many recoverable ones.

The gene axis of every file is compared against the first before any values are
stacked. After stacking, columns are positional — two files disagreeing on gene
order would misalign the matrix silently, and nothing downstream could detect
it.

The single metastatic sample is kept as its own category. Folding it into
either group would misstate that group's median, and the cohort is small enough
that one sample moves it.
"""

from __future__ import annotations

import io
import json
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from car_pipeline.data.source import (
    CacheEntry,
    DataSource,
    IntegrityError,
    _open,
    post_json,
)

PROJECT = "TCGA-PAAD"
WORKFLOW = "STAR - Counts"
MEASURE = "tpm_unstranded"
BATCH_SIZE = 30

FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

PRIMARY_TUMOUR = "primary_tumour"
SOLID_NORMAL = "solid_normal"
METASTATIC = "metastatic"

SAMPLE_TYPES = {
    "Primary Tumor": PRIMARY_TUMOUR,
    "Solid Tissue Normal": SOLID_NORMAL,
    "Metastatic": METASTATIC,
}

# Leading rows of each file are alignment tallies, not genes.
_TALLY_PREFIX = "N_"


@dataclass
class Cohort:
    genes: np.ndarray
    gene_symbols: np.ndarray
    samples: np.ndarray
    sample_types: np.ndarray
    values: np.ndarray  # samples x genes, float32

    def group(self, category: str) -> np.ndarray:
        return self.values[self.sample_types == category]

    def median_for(self, gene_index: int, category: str) -> float:
        block = self.group(category)
        if block.size == 0:
            return float("nan")
        return float(np.median(block[:, gene_index]))


class TCGASource(DataSource):
    name = "GDC TCGA"
    namespace = "tcga"

    def cache_entries(self) -> Iterable[CacheEntry]:
        fp = {
            "project": PROJECT,
            "workflow": WORKFLOW,
            "measure": MEASURE,
        }
        return [
            CacheEntry(key="file_index", filename="file_index.json", fingerprint=fp),
            CacheEntry(key="cohort", filename="cohort.npz", fingerprint=fp),
        ]

    # -- file index -------------------------------------------------------

    def _query_index(self) -> list[dict]:
        payload = {
            "filters": {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "cases.project.project_id",
                            "value": [PROJECT],
                        },
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "data_type",
                            "value": ["Gene Expression Quantification"],
                        },
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "analysis.workflow_type",
                            "value": [WORKFLOW],
                        },
                    },
                ],
            },
            "fields": ",".join(
                [
                    "file_id",
                    "file_name",
                    "cases.samples.sample_type",
                    "cases.samples.submitter_id",
                ]
            ),
            "format": "JSON",
            "size": "10000",
        }
        raw = post_json(FILES_ENDPOINT, payload)
        body = json.loads(raw)
        hits = body["data"]["hits"]
        declared = body["data"]["pagination"]["total"]
        if declared != len(hits):
            raise IntegrityError(
                f"file index: source declared {declared} files, received {len(hits)}"
            )

        rows = []
        for hit in hits:
            sample = hit["cases"][0]["samples"][0]
            rows.append(
                {
                    "file_id": hit["file_id"],
                    "file_name": hit["file_name"],
                    "submitter_id": sample.get("submitter_id", ""),
                    "sample_type": sample["sample_type"],
                }
            )
        rows.sort(key=lambda r: r["file_id"])
        return rows

    def file_index(self) -> list[dict]:
        entry = next(e for e in self.cache_entries() if e.key == "file_index")

        def fetcher(tmp: Path) -> dict:
            print("  querying cohort file index", flush=True)
            rows = self._query_index()
            tmp.write_text(json.dumps(rows, indent=2), encoding="utf-8")
            return {"observed_rows": len(rows), "declared_rows": len(rows)}

        path = self.cache.ensure(entry, fetcher)
        return json.loads(path.read_text(encoding="utf-8"))

    # -- payload ----------------------------------------------------------

    def _download_batch(self, ids: list[str]) -> dict[str, bytes]:
        body = json.dumps({"ids": ids}).encode("utf-8")
        with _open(
            DATA_ENDPOINT,
            method="POST",
            data=body,
            headers={"Content-Type": "application/json"},
            timeout=600,
        ) as resp:
            payload = resp.read()

        out: dict[str, bytes] = {}
        # A single-file request returns the file itself; several return an archive.
        if len(ids) == 1:
            out[ids[0]] = payload
            return out
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                parts = member.name.split("/")
                if len(parts) < 2:
                    continue
                handle = tar.extractfile(member)
                if handle is not None:
                    out[parts[0]] = handle.read()
        return out

    @staticmethod
    def _parse_file(blob: bytes) -> tuple[list[str], list[str], np.ndarray]:
        text = blob.decode("utf-8")
        gene_ids: list[str] = []
        symbols: list[str] = []
        values: list[float] = []

        header: list[str] | None = None
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            row = line.split("\t")
            if header is None:
                header = row
                continue
            if row[0].startswith(_TALLY_PREFIX):
                continue
            if len(row) <= header.index(MEASURE):
                continue
            gene_ids.append(row[0])
            symbols.append(row[1])
            values.append(float(row[header.index(MEASURE)]))

        return gene_ids, symbols, np.asarray(values, dtype=np.float32)

    def build_cohort(self) -> Path:
        entry = next(e for e in self.cache_entries() if e.key == "cohort")
        rows = self.file_index()

        def fetcher(tmp: Path) -> dict:
            axis: list[str] | None = None
            symbols: list[str] | None = None
            samples: list[str] = []
            types: list[str] = []
            blocks: list[np.ndarray] = []

            for start in range(0, len(rows), BATCH_SIZE):
                batch = rows[start : start + BATCH_SIZE]
                ids = [r["file_id"] for r in batch]
                print(
                    f"  batch {start // BATCH_SIZE + 1}"
                    f"/{(len(rows) + BATCH_SIZE - 1) // BATCH_SIZE}"
                    f" ({len(ids)} files)",
                    flush=True,
                )
                payload = self._download_batch(ids)
                if len(payload) != len(ids):
                    raise IntegrityError(
                        f"batch requested {len(ids)} files, archive held {len(payload)}"
                    )
                for meta in batch:
                    blob = payload[meta["file_id"]]
                    gene_ids, syms, values = self._parse_file(blob)
                    if axis is None:
                        axis, symbols = gene_ids, syms
                    elif gene_ids != axis:
                        raise IntegrityError(
                            f"{meta['file_name']}: gene axis differs from the first "
                            "file; values cannot be stacked positionally"
                        )
                    samples.append(meta["submitter_id"] or meta["file_id"])
                    # An unrecognised category keeps its own label rather than
                    # being folded into a catch-all. Folded, it would vanish
                    # from every per-category count while still inflating the
                    # sample total, and the two would stop adding up silently.
                    types.append(
                        SAMPLE_TYPES.get(meta["sample_type"], meta["sample_type"])
                    )
                    blocks.append(values)

            matrix = np.vstack(blocks).astype(np.float32)
            with open(tmp, "wb") as fh:
                np.savez_compressed(
                    fh,
                    genes=np.asarray(axis),
                    gene_symbols=np.asarray(symbols),
                    samples=np.asarray(samples),
                    sample_types=np.asarray(types),
                    values=matrix,
                )
            return {
                "observed_rows": matrix.shape[0],
                "declared_rows": len(rows),
                "extra": {"genes": int(matrix.shape[1])},
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> Cohort:
        path = self.build_cohort()
        with np.load(path, allow_pickle=False) as data:
            return Cohort(
                genes=data["genes"],
                gene_symbols=data["gene_symbols"],
                samples=data["samples"],
                sample_types=data["sample_types"],
                values=data["values"],
            )


JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"


def symbol_index(cohort: Cohort) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, sym in enumerate(cohort.gene_symbols):
        out.setdefault(str(sym), i)
    return out


def ensembl_index(cohort: Cohort) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, gid in enumerate(cohort.genes):
        out.setdefault(str(gid).split(".")[0], i)
    return out


def match_surface(
    cohort: Cohort, surface, atlas_by_accession: dict
) -> dict[str, tuple[int, str]]:
    """Map surface accessions onto column positions, recording the route taken.

    Symbol first; renamed genes are reached through the identifier bridge, which
    this source can support because it carries stable gene identifiers of its
    own. A protein with no column is simply absent — never a column of zeros.
    """
    by_symbol = symbol_index(cohort)
    by_ensembl = ensembl_index(cohort)

    out: dict[str, tuple[int, str]] = {}
    for rec in surface:
        if rec.gene and rec.gene in by_symbol:
            out[rec.accession] = (by_symbol[rec.gene], JOIN_SYMBOL)
            continue
        atlas = atlas_by_accession.get(rec.accession)
        if atlas is not None and atlas.ensembl in by_ensembl:
            out[rec.accession] = (by_ensembl[atlas.ensembl], JOIN_ENSEMBL_BRIDGE)
    return out


if __name__ == "__main__":
    from car_pipeline.data.hpa import HPASource, index as atlas_index
    from car_pipeline.data.uniprot import load_surface

    cohort = TCGASource().load()
    idx = symbol_index(cohort)
    surface, _ = load_surface()
    by_acc, _ = atlas_index(HPASource().load())
    joined = match_surface(cohort, surface, by_acc)
    matched = len(joined)
    bridged = sum(1 for _, p in joined.values() if p == JOIN_ENSEMBL_BRIDGE)

    counts = {
        "samples": (len(cohort.samples), 183),
        "primary tumour": (int((cohort.sample_types == PRIMARY_TUMOUR).sum()), 178),
        "solid normal": (int((cohort.sample_types == SOLID_NORMAL).sum()), 4),
        "metastatic": (int((cohort.sample_types == METASTATIC).sum()), 1),
        "genes": (len(cohort.genes), 60660),
        "surface matched": (matched, 3435),
    }
    for label, (got, exp) in counts.items():
        pct = abs(got - exp) / exp * 100 if exp else 0.0
        print(f"  {label}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    print(f"  joined through the identifier bridge: {bridged}")

    # Every sample must sit in exactly one named category, or the per-category
    # counts and the sample total stop describing the same set of files.
    known = {PRIMARY_TUMOUR, SOLID_NORMAL, METASTATIC}
    stray = sorted({str(t) for t in cohort.sample_types} - known)
    print(f"  categories outside the expected three: {stray or 'none'}")

    print("\n  tumour vs normal median TPM")
    checks = [
        ("CEACAM5", 299.3, 0.7, 409),
        ("MSLN", 298.1, 4.3, 69),
        ("CLDN18", 151.7, 2.3, 65),
        ("CEACAM6", 1194.5, 31.9, 37),
        ("PRSS1", 105.0, 799.0, 0.1),
        ("GPC3", 14.8, 28.5, 0.5),
    ]
    for gene, exp_t, exp_n, exp_ratio in checks:
        i = idx.get(gene)
        if i is None:
            print(f"    {gene}: absent from the gene axis")
            continue
        t = cohort.median_for(i, PRIMARY_TUMOUR)
        n = cohort.median_for(i, SOLID_NORMAL)
        ratio = t / n if n else float("inf")
        print(
            f"    {gene:9s} {t:9.1f} / {n:7.1f} ({ratio:6.1f}x)"
            f"   expected {exp_t:.1f} / {exp_n:.1f} ({exp_ratio}x)"
        )
