"""Tissue atlas connector, pinned to release 23.

The pin is not conservatism. Per-tissue staining and pathology tables were
withdrawn from the per-dataset downloads after this release; later ones ship
only the consolidated gene table, with no per-tissue calls at all. All four
files are taken from the same release so the gene set and the tissue vocabulary
stay consistent with each other.

The level column carries values outside the documented four, and two of them
change the answer:

* Placeholder rows mark genes that were never stained. Counting them as
  coverage inflates the surface-protein figure by roughly a fifth. They are
  dropped at parse, but the gene's symbol, accession and localisation call are
  kept — absence of staining is not absence of the gene.
* Gradient calls on gut enterocytes mean the protein *is* present. Sorting them
  off the end of the scale would place them below "not detected", which
  understates expression — the dangerous direction for a dataset whose whole
  job is flagging risk. They rank with the lowest positive level.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import (
    CacheEntry,
    DataSource,
    extract_from_zip,
    stream_to_file,
)

RELEASE_PIN = "v23"
BASE = f"https://{RELEASE_PIN}.proteinatlas.org/download/"

FILES = {
    "normal_tissue": "normal_tissue.tsv",
    "pathology": "pathology.tsv",
    "subcellular_location": "subcellular_location.tsv",
    "proteinatlas": "proteinatlas.tsv",
}

# Ordinal scale. The gradient calls sit with the lowest positive level: they
# report presence, and ranking them below "not detected" would understate it.
LEVEL_RANK = {
    "Not detected": 0,
    "Low": 1,
    "Ascending": 1,
    "Descending": 1,
    "Medium": 2,
    "High": 3,
}

# Never stained. A placeholder, not a measurement.
PLACEHOLDER_LEVELS = {"N/A", ""}

# A withdrawn call: unscored, and not a detection.
WITHDRAWN_LEVELS = {"Not representative"}

PLASMA_MEMBRANE = "Plasma membrane"


@dataclass
class AtlasGene:
    ensembl: str
    symbol: str = ""
    accession: str = ""
    staining: list[tuple[str, str, int]] = field(default_factory=list)
    main_location: list[str] = field(default_factory=list)
    additional_location: list[str] = field(default_factory=list)

    @property
    def has_staining(self) -> bool:
        return bool(self.staining)

    @property
    def has_subcellular_call(self) -> bool:
        return bool(self.main_location or self.additional_location)

    @property
    def at_plasma_membrane(self) -> bool:
        return PLASMA_MEMBRANE in self.main_location or (
            PLASMA_MEMBRANE in self.additional_location
        )

    def peak_level(self) -> int | None:
        if not self.staining:
            return None
        return max(level for _, _, level in self.staining)


def _split_locations(cell: str) -> list[str]:
    return [p.strip() for p in cell.split(";") if p.strip()]


class HPASource(DataSource):
    name = "Human Protein Atlas"
    namespace = "hpa"

    def cache_entries(self) -> Iterable[CacheEntry]:
        return [
            CacheEntry(
                key=key,
                filename=name,
                fingerprint={"release": RELEASE_PIN, "file": name},
            )
            for key, name in FILES.items()
        ]

    def fetch(self) -> None:
        for entry in self.cache_entries():
            if self.cache.is_valid(entry):
                continue

            def fetcher(tmp: Path, entry: CacheEntry = entry) -> dict:
                archive = tmp.with_suffix(".zip")
                print(f"  fetching {entry.filename}", flush=True)
                stream_to_file(f"{BASE}{entry.filename}.zip", archive)
                meta = extract_from_zip(archive, entry.filename, tmp)
                archive.unlink(missing_ok=True)
                return meta

            self.cache.ensure(entry, fetcher)

    def path_for(self, key: str) -> Path:
        entry = next(e for e in self.cache_entries() if e.key == key)
        if not self.cache.is_valid(entry):
            self.fetch()
        return self.cache.path(entry)

    # -- parsing ----------------------------------------------------------

    def load(self) -> dict[str, AtlasGene]:
        genes: dict[str, AtlasGene] = {}

        def get(ensembl: str) -> AtlasGene:
            g = genes.get(ensembl)
            if g is None:
                g = AtlasGene(ensembl=ensembl)
                genes[ensembl] = g
            return g

        # gene table: symbols and accessions
        with open(self.path_for("proteinatlas"), encoding="utf-8", newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            i_ens = idx.get("Ensembl")
            i_sym = idx.get("Gene")
            i_acc = idx.get("Uniprot")
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if i_ens is None or len(row) <= i_ens:
                    continue
                ensembl = row[i_ens].strip()
                if not ensembl:
                    continue
                g = get(ensembl)
                if i_sym is not None and len(row) > i_sym:
                    g.symbol = row[i_sym].strip()
                if i_acc is not None and len(row) > i_acc:
                    # The column can carry several accessions; the first is the
                    # one the atlas treats as canonical.
                    acc = row[i_acc].strip()
                    g.accession = acc.split(",")[0].strip() if acc else ""

        # normal tissue staining
        dropped_placeholder = 0
        dropped_withdrawn = 0
        with open(self.path_for("normal_tissue"), encoding="utf-8", newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            i_ens = idx.get("Gene")
            i_sym = idx.get("Gene name")
            i_tis = idx.get("Tissue")
            i_cell = idx.get("Cell type")
            i_lvl = idx.get("Level")
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) <= max(i_ens or 0, i_lvl or 0):
                    continue
                level_text = row[i_lvl].strip()
                ensembl = row[i_ens].strip()
                if not ensembl:
                    continue
                g = get(ensembl)
                if not g.symbol and i_sym is not None and len(row) > i_sym:
                    g.symbol = row[i_sym].strip()
                if level_text in PLACEHOLDER_LEVELS:
                    dropped_placeholder += 1
                    continue
                if level_text in WITHDRAWN_LEVELS:
                    dropped_withdrawn += 1
                    continue
                rank = LEVEL_RANK.get(level_text)
                if rank is None:
                    dropped_withdrawn += 1
                    continue
                g.staining.append(
                    (sys.intern(row[i_tis]), sys.intern(row[i_cell]), rank)
                )

        # subcellular localisation
        with open(
            self.path_for("subcellular_location"), encoding="utf-8", newline=""
        ) as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            i_ens = idx.get("Gene")
            i_main = idx.get("Main location")
            i_add = idx.get("Additional location")
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if i_ens is None or len(row) <= i_ens:
                    continue
                ensembl = row[i_ens].strip()
                if not ensembl:
                    continue
                g = get(ensembl)
                if i_main is not None and len(row) > i_main:
                    g.main_location = _split_locations(row[i_main])
                if i_add is not None and len(row) > i_add:
                    g.additional_location = _split_locations(row[i_add])

        self.dropped_placeholder_rows = dropped_placeholder
        self.dropped_withdrawn_rows = dropped_withdrawn
        return genes


def summarise(genes: dict[str, AtlasGene]) -> dict:
    values = list(genes.values())
    return {
        "genes": len(values),
        "with_staining": sum(1 for g in values if g.has_staining),
        "with_subcellular_call": sum(1 for g in values if g.has_subcellular_call),
        "with_accession": sum(1 for g in values if g.accession),
    }


def index(genes: dict[str, AtlasGene]) -> tuple[dict[str, AtlasGene], dict[str, AtlasGene]]:
    """Lookup by accession and by symbol."""
    by_accession: dict[str, AtlasGene] = {}
    by_symbol: dict[str, AtlasGene] = {}
    for g in genes.values():
        if g.accession and g.accession not in by_accession:
            by_accession[g.accession] = g
        if g.symbol and g.symbol not in by_symbol:
            by_symbol[g.symbol] = g
    return by_accession, by_symbol


if __name__ == "__main__":
    from car_pipeline.data.uniprot import load_surface

    src = HPASource()
    atlas = src.load()
    stats = summarise(atlas)

    by_acc, by_sym = index(atlas)
    surface, _ = load_surface()

    matched = []
    for rec in surface:
        g = by_acc.get(rec.accession) or (by_sym.get(rec.gene) if rec.gene else None)
        if g is not None:
            matched.append(g)

    surface_stats = {
        "surface matched": len(matched),
        "with staining": sum(1 for g in matched if g.has_staining),
        "with subcellular call": sum(1 for g in matched if g.has_subcellular_call),
        "at plasma membrane": sum(1 for g in matched if g.at_plasma_membrane),
        "no atlas entry": len(surface) - len(matched),
    }

    expected = {
        "genes": 20162,
        "with_staining": 13468,
        "with_subcellular_call": 13147,
        "with_accession": 19300,
    }
    expected_surface = {
        "surface matched": 3402,
        "with staining": 1945,
        "with subcellular call": 1586,
        "at plasma membrane": 682,
        "no atlas entry": 94,
    }

    print("\natlas")
    for k, v in stats.items():
        exp = expected[k]
        pct = abs(v - exp) / exp * 100
        print(f"  {k}: {v:,}  expected {exp:,}  ({pct:.2f}% off)")
    print("\nsurface intersection")
    for k, v in surface_stats.items():
        exp = expected_surface[k]
        pct = abs(v - exp) / exp * 100 if exp else 0.0
        print(f"  {k}: {v:,}  expected {exp:,}  ({pct:.2f}% off)")
    print(
        f"\nplaceholder rows dropped: {src.dropped_placeholder_rows:,}"
        f"   withdrawn/unscored rows dropped: {src.dropped_withdrawn_rows:,}"
    )
