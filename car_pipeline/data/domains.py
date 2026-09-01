"""Construct parts, fetched by accession and located by feature annotation."""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, _write_json_atomic
from car_pipeline.data.uniprot import RELEASE_PIN

BASE = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "car-platform/stage6"

PROTEOME = "proteome"
SYNTHETIC = "synthetic"
STRUCTURE = "structure"


PARTS = {
    "CD8A": "P01732",
    "TNFRSF9": "Q07011",
    "CD247": "P20963",
    "FKBP1A": "P62942",
    "CASP9": "P55211",
}


HINGE_RESIDUES = 45


@dataclass(frozen=True)
class Part:
    name: str
    provenance: str
    sequence: str
    accession: str = ""
    feature: str = ""
    start: int | None = None
    end: int | None = None

    declared_residues: int | None = None

    @property
    def supplied(self) -> bool:
        """Whether this part carries an actual sequence."""
        return bool(self.sequence)

    @property
    def residues(self) -> int:
        """The part's length in residues, declared or actual."""
        if not self.sequence and self.declared_residues is not None:
            return self.declared_residues
        return len(self.sequence)

    @property
    def bases(self) -> int:
        """The part's coding length in bases."""
        return self.residues * 3

    @property
    def described(self) -> bool:
        """A part must say where it came from. See criterion K4."""
        if self.provenance == SYNTHETIC:
            return bool(self.name)
        if self.provenance == STRUCTURE:
            return bool(self.accession)
        return bool(self.accession and self.start and self.end)


class DomainSource(DataSource):
    name = "UniProt parts"
    namespace = "domains"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The assembled part table and the accessions behind it."""
        return [
            CacheEntry(
                key="parts",
                filename="parts.json",
                fingerprint={
                    "release": RELEASE_PIN,
                    "accessions": sorted(PARTS.values()),
                    "measure": "sequence_and_features",
                },
            )
        ]

    def fetch(self) -> Path:
        """Build the part table from the proteome entries it names."""
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            """Assemble the part table into a temporary file."""
            print("  fetching construct parts", flush=True)
            payload = {}
            for gene, accession in sorted(PARTS.items()):
                url = f"{BASE}/{accession}.json"
                request = urllib.request.Request(
                    url, headers={"User-Agent": USER_AGENT}
                )
                with urllib.request.urlopen(request, timeout=60) as response:
                    data = json.loads(response.read())
                payload[accession] = {
                    "gene": gene,
                    "sequence": data["sequence"]["value"],
                    "features": [
                        {
                            "type": f.get("type", ""),
                            "description": f.get("description", ""),
                            "start": f["location"]["start"].get("value"),
                            "end": f["location"]["end"].get("value"),
                        }
                        for f in data.get("features", [])
                    ],
                }
            _write_json_atomic(tmp, payload)
            blob = json.dumps(payload, sort_keys=True).encode("utf-8")
            import hashlib

            return {
                "digest": hashlib.sha256(blob).hexdigest(),
                "declared_rows": len(PARTS),
                "observed_rows": len(payload),
                "extra": {"source": BASE},
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> dict[str, dict]:
        """The cached construct parts, keyed by role."""
        entry = next(iter(self.cache_entries()))
        if not self.cache.is_valid(entry):
            self.fetch()
        return json.loads(self.cache.path(entry).read_text(encoding="utf-8"))


def _feature(record: dict, kind: str, contains: str = "") -> dict:
    """One named feature's residue range from a proteome entry."""
    for f in record["features"]:
        if f["type"] != kind:
            continue
        if contains and contains.lower() not in (f["description"] or "").lower():
            continue
        if f["start"] and f["end"]:
            return f
    raise LookupError(
        f"{record['gene']}: no {kind} feature"
        + (f" describing {contains!r}" if contains else "")
    )


def _slice(record: dict, start: int, end: int) -> str:
    """Residues start..end inclusive, one-based as the annotation states them."""
    return record["sequence"][start - 1:end]


def build_parts(store: dict[str, dict] | None = None) -> dict[str, Part]:
    """Every proteome-derived part, located from its own entry's features."""
    store = store or DomainSource().load()
    cd8a = store[PARTS["CD8A"]]
    bb = store[PARTS["TNFRSF9"]]
    zeta = store[PARTS["CD247"]]
    fkbp = store[PARTS["FKBP1A"]]
    casp = store[PARTS["CASP9"]]

    signal = _feature(cd8a, "Signal")
    transmem = _feature(cd8a, "Transmembrane")

    hinge_end = transmem["start"] - 1
    hinge_start = hinge_end - HINGE_RESIDUES + 1
    if hinge_start < 1:
        raise LookupError(
            f"CD8A: a {HINGE_RESIDUES}-residue stalk does not fit before the "
            f"transmembrane segment at {transmem['start']}"
        )

    bb_tail = _feature(bb, "Topological domain", "cytoplasmic")
    zeta_tail = _feature(zeta, "Topological domain", "cytoplasmic")
    card = _feature(casp, "Domain", "card")
    fkbp_chain = _feature(fkbp, "Chain")

    out = {
        "leader": Part(
            "CD8A leader", PROTEOME, _slice(cd8a, signal["start"], signal["end"]),
            PARTS["CD8A"], "Signal", signal["start"], signal["end"]),
        "hinge": Part(
            "CD8A hinge", PROTEOME, _slice(cd8a, hinge_start, hinge_end),
            PARTS["CD8A"], "stalk before Transmembrane", hinge_start, hinge_end),
        "transmembrane": Part(
            "CD8A transmembrane", PROTEOME,
            _slice(cd8a, transmem["start"], transmem["end"]),
            PARTS["CD8A"], "Transmembrane", transmem["start"], transmem["end"]),
        "costimulatory": Part(
            "4-1BB cytoplasmic", PROTEOME,
            _slice(bb, bb_tail["start"], bb_tail["end"]),
            PARTS["TNFRSF9"], "Topological domain", bb_tail["start"], bb_tail["end"]),
        "activation": Part(
            "CD3zeta cytoplasmic", PROTEOME,
            _slice(zeta, zeta_tail["start"], zeta_tail["end"]),
            PARTS["CD247"], "Topological domain", zeta_tail["start"], zeta_tail["end"]),
        "switch_fkbp": Part(
            "FKBP12", PROTEOME,
            _slice(fkbp, fkbp_chain["start"], fkbp_chain["end"]),
            PARTS["FKBP1A"], "Chain", fkbp_chain["start"], fkbp_chain["end"]),
        "switch_caspase": Part(
            "caspase-9 without CARD", PROTEOME,
            _slice(casp, card["end"] + 1, len(casp["sequence"])),
            PARTS["CASP9"], "after Domain CARD", card["end"] + 1,
            len(casp["sequence"])),
    }
    return out


SYNTHETIC_PARTS = {
    "linker": Part("(G4S)x3 linker", SYNTHETIC, "GGGGSGGGGSGGGGS"),
    "switch_linker": Part("SGGGS linker", SYNTHETIC, "SGGGS"),
    "skip": Part("T2A skip peptide", SYNTHETIC, "EGRGSLLTCGDVEENPGP"),
    "adaptor_binder": Part(
        "anti-tag binder (sequence not supplied)", SYNTHETIC, "",
        declared_residues=240,
    ),
}


if __name__ == "__main__":
    parts = build_parts()
    print(f"  {'part':24s} {'src':10s} {'accession':10s} {'range':>12s} {'aa':>5s}")
    for key, part in parts.items():
        span = f"{part.start}-{part.end}" if part.start else "-"
        print(f"  {part.name:24s} {part.provenance:10s} {part.accession:10s} "
              f"{span:>12s} {part.residues:5d}")
    for key, part in SYNTHETIC_PARTS.items():
        print(f"  {part.name:24s} {part.provenance:10s} {'':10s} {'-':>12s} "
              f"{part.residues:5d}")


def anti_tag_binder() -> Part:
    """The adaptor binder: retrieved if the structure is cached, declared if not."""
    from car_pipeline.data.antitag import (
        ANTIGEN_ENTITY, AntiTagError, AntiTagSource, BINDER_ENTITIES, ENTRY_ID)

    try:
        source = AntiTagSource()
        payload = source.load()
        sequence = source.sequence()
    except (AntiTagError, OSError, ValueError, KeyError):
        return SYNTHETIC_PARTS["adaptor_binder"]
    if not sequence:
        return SYNTHETIC_PARTS["adaptor_binder"]
    return Part(
        name=(f"anti-tag binder, {payload['tag_system']} "
              f"(PDB {ENTRY_ID} entities {'+'.join(BINDER_ENTITIES)}, "
              f"antigen entity {ANTIGEN_ENTITY} excluded)"),
        provenance=STRUCTURE,
        sequence=sequence,
        accession=f"{ENTRY_ID}_{'+'.join(BINDER_ENTITIES)}",
        feature=f"deposited revision {payload['revision']}",
    )
