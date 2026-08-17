"""Normal tissue transcript baseline, release 10.

Joined on symbol first. A minority of genes were renamed between the proteome
annotation and this one and cannot be reached that way at all; those go through
an identifier bridge built from the tissue atlas gene table. The bridge is
recorded per gene at join time rather than inferred afterwards, because after
the fact a renamed gene and a directly matched one look identical.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from car_pipeline.data.source import CacheEntry, DataSource, stream_to_file

RELEASE_PIN = "v10"
URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
    "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
)

JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"


@dataclass
class TissueProfile:
    ensembl: str
    symbol: str
    values: np.ndarray
    join_path: str

    def max_tpm(self) -> float:
        return float(self.values.max())

    def silent_everywhere(self, threshold: float = 1.0) -> bool:
        """Measured in every tissue and below the threshold in all of them."""
        return bool((self.values < threshold).all())


class GTExSource(DataSource):
    name = "GTEx"
    namespace = "gtex"

    def cache_entries(self) -> Iterable[CacheEntry]:
        return [
            CacheEntry(
                key="gene_median_tpm",
                filename="gene_median_tpm.gct.gz",
                fingerprint={"release": RELEASE_PIN, "measure": "gene_median_tpm"},
            )
        ]

    def fetch(self) -> Path:
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            print("  fetching normal tissue medians", flush=True)
            return stream_to_file(URL, tmp)

        return self.cache.ensure(entry, fetcher)

    def _path(self) -> Path:
        entry = next(iter(self.cache_entries()))
        if not self.cache.is_valid(entry):
            return self.fetch()
        return self.cache.path(entry)

    def tissues(self) -> list[str]:
        with gzip.open(self._path(), "rt", encoding="utf-8") as fh:
            fh.readline()
            fh.readline()
            header = fh.readline().rstrip("\r\n").split("\t")
        return header[2:]

    def match_surface(
        self,
        surface,
        atlas_by_accession: dict,
    ) -> tuple[dict[str, TissueProfile], list[str], int]:
        """Return profiles keyed by accession, the tissue axis, and the gene total.

        A protein with no row here is absent from the result entirely. It is
        never given a row of zeros: not measured and measured at zero are
        different findings, and only one of them is reassuring.
        """
        by_symbol = {}
        by_ensembl = {}
        for rec in surface:
            if rec.gene:
                by_symbol.setdefault(rec.gene, rec)
            atlas = atlas_by_accession.get(rec.accession)
            if atlas is not None and atlas.ensembl:
                by_ensembl.setdefault(atlas.ensembl, rec)

        profiles: dict[str, TissueProfile] = {}
        gene_total = 0

        with gzip.open(self._path(), "rt", encoding="utf-8") as fh:
            fh.readline()
            fh.readline()
            header = fh.readline().rstrip("\r\n").split("\t")
            tissues = header[2:]

            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) < 3:
                    continue
                gene_total += 1
                ensembl = row[0].split(".")[0]
                symbol = row[1]

                rec = by_symbol.get(symbol)
                path = JOIN_SYMBOL
                if rec is None:
                    rec = by_ensembl.get(ensembl)
                    path = JOIN_ENSEMBL_BRIDGE
                if rec is None or rec.accession in profiles:
                    continue

                values = np.asarray(row[2:], dtype=np.float32)
                profiles[rec.accession] = TissueProfile(
                    ensembl=ensembl,
                    symbol=symbol,
                    values=values,
                    join_path=path,
                )

        return profiles, tissues, gene_total


if __name__ == "__main__":
    from car_pipeline.data.hpa import HPASource, index
    from car_pipeline.data.uniprot import load_surface

    surface, _ = load_surface()
    by_acc, _ = index(HPASource().load())

    src = GTExSource()
    profiles, tissues, gene_total = src.match_surface(surface, by_acc)

    bridged = [p for p in profiles.values() if p.join_path == JOIN_ENSEMBL_BRIDGE]

    results = {
        "genes": (gene_total, 59033),
        "tissues": (len(tissues), 68),
        "surface matched": (len(profiles), 3432),
    }
    for label, (got, exp) in results.items():
        pct = abs(got - exp) / exp * 100
        print(f"  {label}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    print(f"\n  joined only through the identifier bridge: {len(bridged)}  expected 29")
    for p in sorted(bridged, key=lambda x: x.symbol)[:12]:
        print(f"    {p.symbol} ({p.ensembl})")
    print(f"  unmatched surface proteins: {len(surface) - len(profiles):,}")
