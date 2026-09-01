"""Genomic span per gene, from a pinned annotation release."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, stream_to_file

RELEASE_PIN = "47"
URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    f"release_{RELEASE_PIN}/gencode.v{RELEASE_PIN}.basic.annotation.gtf.gz"
)

_GENE_NAME = re.compile(r'gene_name "([^"]+)"')


class GeneSpanSource(DataSource):
    name = "GENCODE"
    namespace = "genespan"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The pinned annotation release this source reads."""
        return [
            CacheEntry(
                key="annotation",
                filename=f"gencode.v{RELEASE_PIN}.basic.annotation.gtf.gz",
                fingerprint={"release": RELEASE_PIN, "measure": "gene_span"},
            )
        ]

    def fetch(self) -> Path:
        """Download the annotation if it is absent."""
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            """Stream the annotation into a temporary file."""
            print("  fetching gene annotation", flush=True)
            return stream_to_file(URL, tmp)

        return self.cache.ensure(entry, fetcher)

    def load(self) -> dict[str, int]:
        """Symbol to genomic span in bases."""
        entry = next(iter(self.cache_entries()))
        path = self.cache.path(entry)
        if not self.cache.is_valid(entry):
            path = self.fetch()

        spans: dict[str, int] = {}
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                match = _GENE_NAME.search(parts[8])
                if match is None:
                    continue
                length = int(parts[4]) - int(parts[3]) + 1
                if length > spans.get(match.group(1), 0):
                    spans[match.group(1)] = length
        return spans


if __name__ == "__main__":
    spans = GeneSpanSource().load()
    print(f"  genes with a span: {len(spans):,}")
    for gene in ("NRG3", "MSLN", "CLDN18", "CEACAM5", "MUC1"):
        value = spans.get(gene)
        print(f"    {gene:9s} {value / 1000:9.1f} kb" if value else f"    {gene}: absent")
