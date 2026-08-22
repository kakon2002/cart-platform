"""Genomic span per gene, from a pinned annotation release.

Carried for one reason: the per-cell detection rate this project measures tracks
how long a gene is more strongly than it tracks how much of it is expressed.
Measured over the pool, the rank correlation between detection rate and genomic
span is +0.68 against +0.20 for bulk tumour expression, and the span effect holds
inside every quartile of expression. The cause is known — the cell atlas was
quantified against a pre-mRNA reference, so intronic reads are counted and
intronic content scales with span.

This source does not correct that. It measures it, so a co-expression figure can
be reported beside the span it is confounded with and a reader can tell an
absolutely high overlap from one that is only high for genes of that length.

**Not a blocking dataset and not declared as one in Stage 1.** Nothing gates on
it; it annotates. If a later stage ever gates on a span-derived quantity, it needs
a Stage 1 row and a connector registration at that point, not before.
"""

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
        return [
            CacheEntry(
                key="annotation",
                filename=f"gencode.v{RELEASE_PIN}.basic.annotation.gtf.gz",
                fingerprint={"release": RELEASE_PIN, "measure": "gene_span"},
            )
        ]

    def fetch(self) -> Path:
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            print("  fetching gene annotation", flush=True)
            return stream_to_file(URL, tmp)

        return self.cache.ensure(entry, fetcher)

    def load(self) -> dict[str, int]:
        """Symbol to genomic span in bases.

        The longest annotated span wins where a symbol appears more than once.
        A symbol on more than one contig is a real thing in this annotation and
        the maximum is the conservative reading: it is the figure the capture
        artefact would scale with.
        """
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
