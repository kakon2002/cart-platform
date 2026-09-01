"""Reviewed human proteome and the surface filter built on it."""

from __future__ import annotations

import re
import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import (
    CacheEntry,
    CacheError,
    DataSource,
    stream_paginated_to_file,
)

QUERY = "(reviewed:true) AND (organism_id:9606)"
FIELDS = [
    "accession",
    "gene_primary",
    "protein_name",
    "cc_subcellular_location",
    "ft_transmem",
    "ft_topo_dom",
    "ft_lipid",
    "ft_chain",
]


RELEASE_PIN = "2026_02"
RELEASE_HEADER = "X-UniProt-Release"
RELEASE_DATE_HEADER = "X-UniProt-Release-Date"
PAGE_SIZE = 500
BASE = "https://rest.uniprot.org/uniprotkb/search"

_NOTE = re.compile(r'/note="([^"]*)"')
_TRANSMEM = re.compile(r"\bTRANSMEM\b")
_TOPO_SEGMENT = re.compile(r"TOPO_DOM\s+(\S+?)\.\.(\S+?);((?:(?!TOPO_DOM).)*)", re.S)
_CHAIN_SEGMENT = re.compile(r"CHAIN\s+(\S+?)\.\.(\S+?);((?:(?!CHAIN).)*)", re.S)
_GPI = re.compile(r"GPI-anchor", re.IGNORECASE)


_PLASMA_MEMBRANE_TERMS = ("cell membrane", "cell surface")


_NOTE_FIELD = "Note="
_LOCATION_BLOCK = "SUBCELLULAR LOCATION:"


def location_statements(subcellular: str) -> str:
    """The localisation text with every free-text note removed."""
    kept = []
    for block in subcellular.split(_LOCATION_BLOCK):
        if not block.strip():
            continue
        cut = block.find(_NOTE_FIELD)
        kept.append(block if cut == -1 else block[:cut])
    return " ".join(kept)


_OUTWARD_NOTE = "Extracellular"


_INTRACELLULAR_COMPARTMENTS = (
    "endoplasmic reticulum",
    "golgi",
    "mitochondri",
    "lysosom",
    "peroxisom",
    "nucleus",
    "nuclear",
    "endosom",
    "vesicle",
    "cytoplasm",
    "cytosol",
    "sarcoplasmic reticulum",
    "microsom",
    "melanosom",
    "phagosom",
    "autophagosom",
    "vacuol",
    "secretory",
    "synaptic",
    "acrosom",
    "centrosom",
    "membrane raft",
    "lipid droplet",
    "midbody",
    "spindle",
    "endomembrane",
    "photoreceptor",
    "recycling",
    "secreted",
    "dendrite",
)


class MembraneClass:
    SINGLE_PASS = "single_pass"
    MULTI_PASS = "multi_pass"
    GPI_ANCHORED = "gpi_anchored"


@dataclass
class ProteinRecord:
    accession: str
    gene: str
    protein_name: str
    subcellular: str
    transmem_count: int
    gpi_anchored: bool
    topo_notes: list[str] = field(default_factory=list)
    extracellular_residues: int | None = None

    chains: list["Chain"] = field(default_factory=list)

    attached: bool = False
    outward: bool = False
    membrane_class: str | None = None

    outward_note_only: bool = False

    @property
    def is_surface(self) -> bool:
        """Whether both surface gates admit this entry."""
        return self.attached and self.outward


def _count_extracellular_residues(topo_field: str) -> int | None:
    """Total annotated outward-facing residues, or None when never annotated."""
    total = 0
    measured = False
    for start, end, tail in _TOPO_SEGMENT.findall(topo_field):
        note_match = _NOTE.search(tail)
        if not note_match or note_match.group(1) != _OUTWARD_NOTE:
            continue
        try:
            lo = int(start.lstrip("<>?"))
            hi = int(end.lstrip("<>?"))
        except ValueError:
            continue
        if hi >= lo:
            total += hi - lo + 1
            measured = True

    return total if measured else None


_UNCERTAIN = "<>?"
_CHAIN_ID = re.compile(r'/id="([^"]*)"')


def _bound(text: str) -> tuple[int | None, bool]:
    """A residue position and whether the annotation hedged it."""
    stripped = text.lstrip(_UNCERTAIN)
    uncertain = stripped != text
    try:
        return int(stripped), uncertain
    except ValueError:
        return None, True


@dataclass(frozen=True)
class Chain:
    """One mature chain carved out of a precursor."""

    start: int | None
    end: int | None
    note: str
    chain_id: str
    start_uncertain: bool = False
    end_uncertain: bool = False

    @property
    def exact(self) -> bool:
        """True only when both boundaries are numbers the annotation stated flatly."""
        return (
            self.start is not None
            and self.end is not None
            and not self.start_uncertain
            and not self.end_uncertain
        )

    def contains(self, position: int) -> bool:
        """Whether a residue falls inside, refusing to guess on a missing bound."""
        if self.start is None or self.end is None:
            return False
        return self.start <= position <= self.end

    @property
    def length(self) -> int | None:
        """The chain's length, or nothing where a bound is missing."""
        if self.start is None or self.end is None:
            return None
        return self.end - self.start + 1


def parse_chains(chain_field: str) -> list[Chain]:
    """Mature chains carved out of the precursor, in annotation order."""
    chains: list[Chain] = []
    for start, end, tail in _CHAIN_SEGMENT.findall(chain_field):
        note_match = _NOTE.search(tail)
        id_match = _CHAIN_ID.search(tail)
        lo, lo_unc = _bound(start)
        hi, hi_unc = _bound(end)
        chains.append(
            Chain(
                start=lo,
                end=hi,
                note=note_match.group(1) if note_match else "",
                chain_id=id_match.group(1) if id_match else "",
                start_uncertain=lo_unc,
                end_uncertain=hi_unc,
            )
        )
    return chains


def names_compartment(subcellular: str) -> bool:
    """True when the localisation text places the protein somewhere specific."""
    text = subcellular.lower()
    return any(c in text for c in _INTRACELLULAR_COMPARTMENTS)


def parse_row(row: list[str]) -> ProteinRecord:
    """One proteome row into a record."""
    accession, gene, protein_name, subcellular, transmem, topo, lipid, chain = (
        (row + [""] * len(FIELDS))[: len(FIELDS)]
    )

    transmem_count = len(_TRANSMEM.findall(transmem))
    gpi = bool(_GPI.search(lipid))
    notes = _NOTE.findall(topo)

    rec = ProteinRecord(
        accession=accession,
        gene=gene,
        protein_name=protein_name,
        subcellular=subcellular,
        transmem_count=transmem_count,
        gpi_anchored=gpi,
        topo_notes=notes,
        extracellular_residues=_count_extracellular_residues(topo),
        chains=parse_chains(chain),
    )

    rec.attached = transmem_count > 0 or gpi

    statements = location_statements(subcellular).lower()
    rec.outward = (
        gpi
        or _OUTWARD_NOTE in notes
        or any(term in statements for term in _PLASMA_MEMBRANE_TERMS)
    )
    if not rec.outward and rec.attached:
        full = subcellular.lower()
        rec.outward_note_only = any(t in full for t in _PLASMA_MEMBRANE_TERMS)

    if rec.is_surface:
        if gpi:
            rec.membrane_class = MembraneClass.GPI_ANCHORED
        elif transmem_count == 1:
            rec.membrane_class = MembraneClass.SINGLE_PASS
        else:
            rec.membrane_class = MembraneClass.MULTI_PASS

    return rec


def check_release(meta: dict) -> str:
    """Confirm the service served the release this project is pinned to."""
    served = meta.get("extra", {}).get(RELEASE_HEADER.lower())
    if served is None:
        raise CacheError(
            "the proteome service did not state which release it served; "
            f"refusing to cache an unlabelled fetch under {RELEASE_PIN}"
        )
    if served != RELEASE_PIN:
        raise CacheError(
            f"the proteome service is serving release {served} and this project "
            f"is pinned to {RELEASE_PIN}. Every count downstream is measured "
            "against the pinned release, so the fetch stops here rather than "
            "replacing them with a different proteome. Bump RELEASE_PIN "
            "deliberately and re-run every verifier."
        )
    return served


class UniProtSource(DataSource):
    name = "UniProt"
    namespace = "uniprot"

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The pinned proteome release and the fields requested from it."""
        return [
            CacheEntry(
                key="human_reviewed",
                filename="human_reviewed.tsv",
                fingerprint={
                    "query": QUERY,
                    "fields": FIELDS,
                    "release": RELEASE_PIN,
                    "format": "tsv",
                },
            )
        ]

    def _url(self) -> str:
        """The query URL for the pinned field list."""
        params = {
            "query": QUERY,
            "fields": ",".join(FIELDS),
            "format": "tsv",
            "size": PAGE_SIZE,
        }
        return f"{BASE}?{urllib.parse.urlencode(params)}"

    def fetch(self) -> Path:
        """Download the proteome if it is absent."""
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            """Stream the proteome into a temporary file."""
            print("  fetching reviewed human proteome", flush=True)
            meta = stream_paginated_to_file(
                self._url(),
                tmp,
                progress_label="proteome",
                capture_headers=(RELEASE_HEADER, RELEASE_DATE_HEADER),
            )
            served = check_release(meta)
            print(f"    proteome release {served} confirmed against the pin", flush=True)
            return meta

        return self.cache.ensure(entry, fetcher)

    def load(self) -> list[ProteinRecord]:
        """Every reviewed entry, parsed."""
        entry = next(iter(self.cache_entries()))
        path = self.cache.path(entry)
        if not self.cache.is_valid(entry):
            path = self.fetch()

        records: list[ProteinRecord] = []
        with open(path, "r", encoding="utf-8", newline="") as fh:
            header = fh.readline()
            if not header:
                return records
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                records.append(parse_row(line.split("\t")))
        return records


def summarise(records: list[ProteinRecord]) -> dict:
    """Counts of what the surface filter admitted and withheld."""
    surface = [r for r in records if r.is_surface]
    note_only = [r for r in records if r.outward_note_only]
    attached = [r for r in records if r.attached]
    withheld = [r for r in attached if not r.outward]

    internal = [r for r in withheld if names_compartment(r.subcellular)]

    unresolved = [r for r in withheld if not names_compartment(r.subcellular)]

    return {
        "entries": len(records),
        "surface": len(surface),
        "single_pass": sum(
            1 for r in surface if r.membrane_class == MembraneClass.SINGLE_PASS
        ),
        "multi_pass": sum(
            1 for r in surface if r.membrane_class == MembraneClass.MULTI_PASS
        ),
        "gpi_anchored": sum(
            1 for r in surface if r.membrane_class == MembraneClass.GPI_ANCHORED
        ),
        "internal_anchored": len(internal),
        "compartment_unresolved": len(unresolved),
        "outward_note_only": len(note_only),
    }


def load_surface() -> tuple[list[ProteinRecord], dict]:
    """The surface-accessible proteome and the counts behind it."""
    records = UniProtSource().load()
    return [r for r in records if r.is_surface], summarise(records)


if __name__ == "__main__":
    recs = UniProtSource().load()
    stats = summarise(recs)

    expected = {
        "entries": 20431,
        "surface": 3466,
        "single_pass": 1446,
        "multi_pass": 1884,
        "gpi_anchored": 136,
        "internal_anchored": 1362,
        "compartment_unresolved": 534,
        "outward_note_only": 14,
    }
    for k, v in stats.items():
        exp = expected[k]
        print(f"  {'ok  ' if v == exp else 'DIFF'}  {k}: {v:,}  expected {exp:,}")
