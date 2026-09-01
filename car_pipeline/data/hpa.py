"""Tissue atlas connector, pinned to release 23."""

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


LEVEL_RANK = {
    "Not detected": 0,
    "Low": 1,
    "Ascending": 1,
    "Descending": 1,
    "Medium": 2,
    "High": 3,
}


PLACEHOLDER_LEVELS = {"N/A", ""}


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
        """Whether any normal tissue was stained for this entry."""
        return bool(self.staining)

    @property
    def has_subcellular_call(self) -> bool:
        """Whether a subcellular location was recorded."""
        return bool(self.main_location or self.additional_location)

    @property
    def at_plasma_membrane(self) -> bool:
        """Whether a location statement places it at the plasma membrane."""
        return PLASMA_MEMBRANE in self.main_location or (
            PLASMA_MEMBRANE in self.additional_location
        )

    def peak_level(self) -> int | None:
        """The highest staining level recorded across tissues."""
        if not self.staining:
            return None
        return max(level for _, _, level in self.staining)


def _split_locations(cell: str) -> list[str]:
    """Split a location field into its separate statements."""
    return [p.strip() for p in cell.split(";") if p.strip()]


def _columns(header: list[str], required: list[str], filename: str) -> dict[str, int]:
    """Resolve required column positions, refusing to proceed without them."""
    idx = {name: i for i, name in enumerate(header)}
    missing = [c for c in required if c not in idx]
    if missing:
        raise KeyError(f"{filename} is missing columns: {', '.join(missing)}")
    return idx


class HPASource(DataSource):
    name = "Human Protein Atlas"
    namespace = "hpa"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The pinned atlas release this source reads."""
        return [
            CacheEntry(
                key=key,
                filename=name,
                fingerprint={"release": RELEASE_PIN, "file": name},
            )
            for key, name in FILES.items()
        ]

    def fetch(self) -> None:
        """Download the atlas if it is absent."""
        for entry in self.cache_entries():
            if self.cache.is_valid(entry):
                continue

            def fetcher(tmp: Path, entry: CacheEntry = entry) -> dict:
                """Stream the atlas into a temporary file."""
                archive = tmp.with_suffix(".zip")
                print(f"  fetching {entry.filename}", flush=True)
                stream_to_file(f"{BASE}{entry.filename}.zip", archive)
                meta = extract_from_zip(archive, entry.filename, tmp)
                archive.unlink(missing_ok=True)
                return meta

            self.cache.ensure(entry, fetcher)

    def path_for(self, key: str) -> Path:
        """The cached path for one key, fetching first if it is absent."""
        entry = next(e for e in self.cache_entries() if e.key == key)
        if not self.cache.is_valid(entry):
            self.fetch()
        return self.cache.path(entry)

    def load(self) -> dict[str, AtlasGene]:
        """Every atlas entry, parsed."""
        genes: dict[str, AtlasGene] = {}

        def get(ensembl: str) -> AtlasGene:
            """One atlas entry by accession."""
            g = genes.get(ensembl)
            if g is None:
                g = AtlasGene(ensembl=ensembl)
                genes[ensembl] = g
            return g

        with open(self.path_for("proteinatlas"), encoding="utf-8", newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = _columns(header, ["Ensembl", "Gene", "Uniprot"], "gene table")
            i_ens = idx["Ensembl"]
            i_sym = idx["Gene"]
            i_acc = idx["Uniprot"]
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) <= i_ens:
                    continue
                ensembl = row[i_ens].strip()
                if not ensembl:
                    continue
                g = get(ensembl)
                if i_sym is not None and len(row) > i_sym:
                    g.symbol = row[i_sym].strip()
                if i_acc is not None and len(row) > i_acc:
                    acc = row[i_acc].strip()
                    g.accession = acc.split(",")[0].strip() if acc else ""

        dropped_placeholder = 0
        dropped_withdrawn = 0
        with open(self.path_for("normal_tissue"), encoding="utf-8", newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = _columns(
                header,
                ["Gene", "Gene name", "Tissue", "Cell type", "Level"],
                "normal tissue",
            )
            i_ens = idx["Gene"]
            i_sym = idx["Gene name"]
            i_tis = idx["Tissue"]
            i_cell = idx["Cell type"]
            i_lvl = idx["Level"]
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) <= max(i_ens, i_lvl, i_tis, i_cell):
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

        with open(
            self.path_for("subcellular_location"), encoding="utf-8", newline=""
        ) as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            idx = _columns(
                header,
                ["Gene", "Main location", "Additional location"],
                "subcellular location",
            )
            i_ens = idx["Gene"]
            i_main = idx["Main location"]
            i_add = idx["Additional location"]
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) <= i_ens:
                    continue
                ensembl = row[i_ens].strip()
                if not ensembl:
                    continue
                g = get(ensembl)
                if len(row) > i_main:
                    g.main_location = _split_locations(row[i_main])
                if len(row) > i_add:
                    g.additional_location = _split_locations(row[i_add])

        self.dropped_placeholder_rows = dropped_placeholder
        self.dropped_withdrawn_rows = dropped_withdrawn
        return genes


def summarise(genes: dict[str, AtlasGene]) -> dict:
    """Counts of what the atlas records across the entries read."""
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
