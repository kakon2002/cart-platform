"""Antibody annotation over deposited structures, and named therapeutics.

Two sources, cached under the usual discipline.

**The curated structure summary** says which chain of a complex is heavy, which
is light and which is the antigen. Without it those have to be guessed from
entity description text, which the specification called the weak link. It is a
bulk CSV; the front end that serves it is a JavaScript application, and an
earlier probe read that shell as the source being unreachable. It is not — the
application has a REST interface and this is one of its routes.

**The therapeutic table** carries named antibodies with their clinical stage and,
the part that matters most, their actual variable-region sequences. A binder does
not have to be a solved structure to be usable; a sequence is what a construct is
built from.

Target matching on the therapeutic table is by **token**, never by substring. The
field is compound and synonym-laden — measured values include `CEACAM5/CD66e`,
`MUC1/PEM/EMA`, `CLDN18;CD3E` and `IAP/CD47;CLDN18` — so the field is split on
`;` to separate the antigens of a bispecific and on `/` to separate synonyms,
then matched exactly. Substring matching would put MUC16 and MUC18 in the MUC1
bucket, and would match `CLDN1` inside `CLDN18`, which this table contains as
separate targets.
"""

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
        return bool(self.heavy_sequence and self.light_sequence)


#: The two identifier formats these sources use for the same entry. The curated
#: summary keys on the extended form, `pdb_00004f3f`; the structure search returns
#: the short form, `4F3F`. Joining one directly against the other matches nothing
#: — and matches nothing *silently*, producing an empty candidate list for every
#: target that reads exactly like "no binder exists". It did, for all 200, until a
#: positive known answer was added to catch it.
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
        return next(e for e in self.cache_entries() if e.key == key)

    def fetch(self) -> Path:
        for key, url, label in (
            ("structure_summary", SUMMARY_URL, "antibody structure summary"),
            ("therapeutics", THERAPEUTIC_URL, "therapeutic antibodies"),
        ):
            entry = self._entry(key)
            if self.cache.is_valid(entry):
                continue

            def fetcher(tmp: Path, _url=url, _label=label) -> dict:
                print(f"  fetching {_label}", flush=True)
                return stream_to_file(_url, tmp)

            self.cache.ensure(entry, fetcher)
        return self.cache.path(self._entry("structure_summary"))

    def _rows(self, key: str) -> list[dict]:
        entry = self._entry(key)
        if not self.cache.is_valid(entry):
            self.fetch()
        text = self.cache.path(entry).read_text(encoding="utf-8-sig", errors="replace")
        return list(csv.DictReader(io.StringIO(text)))

    def structures(self) -> dict[str, list[AntibodyStructure]]:
        """Antibody instances keyed by the four-character entry code, lower case.

    Keyed on the short form because that is what the structure search returns;
    the extended form is retained on each instance so the row can be traced back.
    """
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
