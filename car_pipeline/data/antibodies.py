"""Antibody annotation over deposited structures, and named therapeutics."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, stream_to_file

SUMMARY_URL = "https://sabdab.opig.stats.ox.ac.uk/api/download/all-summary"
THERAPEUTIC_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
    "TheraSAbDab_SeqStruc_OnlineDownload.csv"
)
RELEASE_PIN = "2.0.10"


@dataclass
class AntibodyStructure:
    """One antibody instance in one deposited entry."""

    pdb: str
    extended_id: str
    heavy_chain: str
    light_chain: str
    antigen_chain: str
    antigen_type: str
    antigen_name: str
    method: str
    resolution: str


@dataclass
class Therapeutic:
    name: str
    fmt: str
    targets: list[str]
    highest_trial: str
    status: str
    heavy_sequence: str
    light_sequence: str
    structures: str = ""
    conditions: str = ""

    @property
    def has_sequence(self) -> bool:
        """Whether both variable regions are present."""
        return bool(self.heavy_sequence and self.light_sequence)


_EXTENDED_PREFIX = "pdb_0000"


def short_entry_id(identifier: str) -> str:
    """The four-character entry code, from either identifier format."""
    value = identifier.strip().lower()
    if value.startswith(_EXTENDED_PREFIX):
        value = value[len(_EXTENDED_PREFIX):]
    return value


def split_targets(field_value: str) -> list[str]:
    """Tokens of a compound target field, exact and de-duplicated."""
    tokens: list[str] = []
    for antigen in field_value.split(";"):
        for synonym in antigen.split("/"):
            token = synonym.strip()
            if token and token not in tokens:
                tokens.append(token)
    return tokens


class AntibodySource(DataSource):
    name = "SAbDab"
    namespace = "antibodies"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The two tables this source caches, and the release that defines them."""
        return [
            CacheEntry(
                key="structure_summary",
                filename="sabdab_summary_all.csv",
                fingerprint={"release": RELEASE_PIN, "measure": "structure_summary"},
            ),
            CacheEntry(
                key="therapeutics",
                filename="therasabdab_seqstruc.csv",
                fingerprint={"release": RELEASE_PIN, "measure": "therapeutics"},
            ),
        ]

    def _entry(self, key: str) -> CacheEntry:
        """One cache entry by key."""
        return next(e for e in self.cache_entries() if e.key == key)

    def fetch(self) -> Path:
        """Download whichever of the two tables is missing."""
        for key, url, label in (
            ("structure_summary", SUMMARY_URL, "antibody structure summary"),
            ("therapeutics", THERAPEUTIC_URL, "therapeutic antibodies"),
        ):
            entry = self._entry(key)
            if self.cache.is_valid(entry):
                continue

            def fetcher(tmp: Path, _url=url, _label=label) -> dict:
                """Stream one table into a temporary file."""
                print(f"  fetching {_label}", flush=True)
                return stream_to_file(_url, tmp)

            self.cache.ensure(entry, fetcher)
        return self.cache.path(self._entry("structure_summary"))

    def _rows(self, key: str) -> list[dict]:
        """One cached table as dict rows, fetching first if it is absent."""
        entry = self._entry(key)
        if not self.cache.is_valid(entry):
            self.fetch()
        text = self.cache.path(entry).read_text(encoding="utf-8-sig", errors="replace")
        return list(csv.DictReader(io.StringIO(text)))

    def structures(self) -> dict[str, list[AntibodyStructure]]:
        """Antibody instances keyed by the four-character entry code, lower case."""
        out: dict[str, list[AntibodyStructure]] = {}
        for row in self._rows("structure_summary"):
            raw = (row.get("PDB") or "").strip()
            if not raw:
                continue
            pdb = short_entry_id(raw)
            out.setdefault(pdb, []).append(
                AntibodyStructure(
                    pdb=pdb,
                    extended_id=raw,
                    heavy_chain=(row.get("Hchain") or "").strip(),
                    light_chain=(row.get("Lchain") or "").strip(),
                    antigen_chain=(row.get("antigen_chain") or "").strip(),
                    antigen_type=(row.get("antigen_type") or "").strip(),
                    antigen_name=(row.get("antigen_name") or "").strip(),
                    method=(row.get("method") or "").strip(),
                    resolution=(row.get("resolution") or "").strip(),
                )
            )
        return out

    def therapeutics(self) -> list[Therapeutic]:
        """The therapeutic antibodies, one record per named entry."""
        out: list[Therapeutic] = []
        for row in self._rows("therapeutics"):
            name = (row.get("Therapeutic") or "").strip()
            if not name:
                continue
            trial = ""
            for key in row:
                if key.startswith("Highest_Clin_Trial"):
                    trial = (row[key] or "").strip()
                    break
            out.append(
                Therapeutic(
                    name=name,
                    fmt=(row.get("Format") or "").strip(),
                    targets=split_targets(row.get("Target") or ""),
                    highest_trial=trial,
                    status=(row.get("Est. Status") or "").strip(),
                    heavy_sequence=(row.get("HeavySequence") or "").strip(),
                    light_sequence=(row.get("LightSequence") or "").strip(),
                    structures=(row.get("100% SI Structure") or "").strip(),
                    conditions=(row.get("Conditions Approved") or "").strip(),
                )
            )
        return out

    def therapeutics_by_target(self) -> dict[str, list[Therapeutic]]:
        """Therapeutics indexed by every target token they name."""
        index: dict[str, list[Therapeutic]] = {}
        for therapeutic in self.therapeutics():
            for token in therapeutic.targets:
                index.setdefault(token, []).append(therapeutic)
        return index


if __name__ == "__main__":
    source = AntibodySource()
    structures = source.structures()
    therapeutics = source.therapeutics()
    index = source.therapeutics_by_target()
    print(f"  entries with an antibody instance   {len(structures):,}")
    for probe in ("4f3f", "9v32", "8bw0"):
        print(f"    {probe}: {len(structures.get(probe, []))} instance(s)")
    print(f"  named therapeutics                  {len(therapeutics):,}")
    print(f"  distinct target tokens              {len(index):,}")
    for gene in ("MSLN", "CLDN18", "CEACAM5", "CEACAM6", "MUC1"):
        hits = index.get(gene, [])
        print(f"    {gene:9s} {len(hits):2d} therapeutic(s): "
              f"{', '.join(t.name for t in hits[:4])}")
