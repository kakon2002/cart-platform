"""Construct parts, fetched by accession and located by feature annotation.

No coordinate in this module is transcribed from memory. Each part names the
annotated feature it is taken from — a signal peptide, a transmembrane segment, a
cytoplasmic topological domain — and the range is read from the entry. A part
whose feature is absent raises rather than falling back to a remembered number,
because a remembered number that is close enough to look right is the failure
this project keeps finding.

Two provenance classes, kept apart and never blurred:

* **proteome** — an accession, a residue range and the release pin. Re-derivable.
* **synthetic** — a designed sequence with no database entry: the flexible linker
  and the ribosomal skip peptide. Recorded as a literal with its name.

Fetched per accession rather than by widening the proteome query, which would
invalidate that cache and re-fetch twenty thousand entries for eight sequences.
"""

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

#: Accessions this stage draws from, and what each supplies.
PARTS = {
    "CD8A": "P01732",
    "TNFRSF9": "Q07011",
    "CD247": "P20963",
    "FKBP1A": "P62942",
    "CASP9": "P55211",
}

#: The membrane-proximal stalk used as a hinge. There is no "hinge" feature to
#: read, so its length is a stated design choice and its position is derived:
#: the segment immediately preceding the annotated transmembrane start.
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

    @property
    def residues(self) -> int:
        return len(self.sequence)

    @property
    def bases(self) -> int:
        return len(self.sequence) * 3

    @property
    def described(self) -> bool:
        """A part must say where it came from. See criterion K4."""
        if self.provenance == SYNTHETIC:
            return bool(self.name)
        return bool(self.accession and self.start and self.end)


class DomainSource(DataSource):
    name = "UniProt parts"
    namespace = "domains"

    def cache_entries(self) -> Iterable[CacheEntry]:
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
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
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
        entry = next(iter(self.cache_entries()))
        if not self.cache.is_valid(entry):
            self.fetch()
        return json.loads(self.cache.path(entry).read_text(encoding="utf-8"))


def _feature(record: dict, kind: str, contains: str = "") -> dict:
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
    # The stalk immediately before the membrane. Length is the design choice in
    # HINGE_RESIDUES; the position is read from the annotated segment.
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
        # The switch is FKBP12 fused to caspase-9 with its CARD removed, which is
        # what makes it dimeriser-inducible rather than constitutively active.
        # The CARD boundary is read from the entry, not assumed.
        # Read from the annotated chain rather than assumed to be the whole
        # entry: the mature protein starts after the initiator methionine, and
        # this was the one part bypassing the feature lookup every other uses.
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


#: Designed sequences with no database entry. Named, marked, and never given a
#: fabricated accession.
SYNTHETIC_PARTS = {
    "linker": Part("(G4S)x3 linker", SYNTHETIC, "GGGGSGGGGSGGGGS"),
    "switch_linker": Part("SGGGS linker", SYNTHETIC, "SGGGS"),
    "skip": Part("T2A skip peptide", SYNTHETIC, "EGRGSLLTCGDVEENPGP"),
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
