"""Reviewed human proteome and the surface filter built on it.

The seventh requested field, lipidation, is not optional. Anchors of that kind
leave no transmembrane segment behind and are not consistently spelled out in
the localisation text, so without it that entire class of protein is invisible
to the attachment gate — which would silently drop several of the best known
targets.

The eighth, chain boundaries, is requested because a mature protein is not
always one molecule. Where a precursor is cleaved, some of the resulting chains
are released rather than held at the surface, and a binder raised against a
released chain meets its antigen in plasma rather than on a cell. The surface
filter itself does not consult this field — attachment and orientation are
decided by the membrane evidence alone — so adding it left the surface set
identical, 3,480 before and after. It is not free even so: the field list is
inside the cache fingerprint, so adding it invalidated the proteome cache and
re-fetched all 20,431 entries. The field is carried so that a later stage can
tell the anchored chain from the shed one instead of treating the precursor as
uniformly reachable.
"""

from __future__ import annotations

import re
import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, stream_paginated_to_file

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
PAGE_SIZE = 500
BASE = "https://rest.uniprot.org/uniprotkb/search"

_NOTE = re.compile(r'/note="([^"]*)"')
_TRANSMEM = re.compile(r"\bTRANSMEM\b")
_TOPO_SEGMENT = re.compile(r"TOPO_DOM\s+(\S+?)\.\.(\S+?);((?:(?!TOPO_DOM).)*)", re.S)
_CHAIN_SEGMENT = re.compile(r"CHAIN\s+(\S+?)\.\.(\S+?);((?:(?!CHAIN).)*)", re.S)
_GPI = re.compile(r"GPI-anchor", re.IGNORECASE)

# Localisation phrases that place a protein at the outward face of the cell
# itself. Matched case-insensitively as substrings of the localisation text so
# that polarised variants ("Apical cell membrane") are covered.
_PLASMA_MEMBRANE_TERMS = ("cell membrane", "cell surface")

# The one topological note that means the outward face of the cell. Matched
# exactly: "Lumenal" and "Perinuclear space" sit on the same side of the bilayer
# but inside an organelle, and would be wrongly admitted by a substring test.
_OUTWARD_NOTE = "Extracellular"

# Compartments that place an anchored protein inside the cell. Used only to
# separate "measured and placed internally" from "nothing on record", which are
# reported apart rather than pooled: the second group is held out for want of
# evidence, not judged against.
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

    #: Mature chains carved out of the precursor, as (start, end, note). A bound
    #: recorded as uncertain is kept as None rather than dropping the chain: the
    #: chain exists either way, and losing the row would understate how many
    #: pieces a precursor is cut into. An empty list means the entry carries no
    #: chain annotation at all, which is not the same as being a single chain.
    chains: list[tuple[int | None, int | None, str]] = field(default_factory=list)

    # filter outcome
    attached: bool = False
    outward: bool = False
    membrane_class: str | None = None

    @property
    def is_surface(self) -> bool:
        return self.attached and self.outward


def _count_extracellular_residues(topo_field: str) -> int | None:
    """Total annotated outward-facing residues, or None when never annotated.

    None is a third state on purpose. A protein nobody annotated must not be
    scored as though it were measured and found small.
    """
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
            # Bounds recorded as uncertain. The segment exists but its length
            # does not, so it contributes nothing and cannot make the total
            # meaningful on its own.
            continue
        if hi >= lo:
            total += hi - lo + 1
            measured = True
    # A protein annotated as outward-facing but with no segment whose length
    # could be read is not a protein measured at zero residues. Returning zero
    # here would let it be scored as measured and tiny, which is the imputation
    # the accessibility component is required to avoid.
    return total if measured else None


def _bound(text: str) -> int | None:
    """A residue position, or None where the annotation records it as uncertain."""
    try:
        return int(text.lstrip("<>?"))
    except ValueError:
        return None


def parse_chains(chain_field: str) -> list[tuple[int | None, int | None, str]]:
    """Mature chains carved out of the precursor, in annotation order.

    A chain whose bounds cannot be read is kept with None bounds rather than
    skipped. Dropping it would make a cleaved precursor look like an uncleaved
    one, which is the direction that matters: it is the difference between a
    protein held at the surface and one released into plasma.
    """
    chains: list[tuple[int | None, int | None, str]] = []
    for start, end, tail in _CHAIN_SEGMENT.findall(chain_field):
        note_match = _NOTE.search(tail)
        chains.append(
            (_bound(start), _bound(end), note_match.group(1) if note_match else "")
        )
    return chains


def names_compartment(subcellular: str) -> bool:
    """True when the localisation text places the protein somewhere specific."""
    text = subcellular.lower()
    return any(c in text for c in _INTRACELLULAR_COMPARTMENTS)


def parse_row(row: list[str]) -> ProteinRecord:
    # Width tracks FIELDS. A row padded or truncated to the old width would drop
    # the new column silently, and the drop would only surface as an empty chain
    # list — indistinguishable from a protein that genuinely has no chain
    # annotation, and only after the re-fetch had already been paid for.
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

    # Gate 1 — is it held in a membrane at all. Secreted proteins fail here.
    rec.attached = transmem_count > 0 or gpi

    # Gate 2 — positive evidence that it faces outward. Stated as evidence to
    # look for rather than compartments to exclude: that is what keeps out the
    # multi-pass proteins of internal compartments, which are topologically
    # outward facing but unreachable from outside the cell.
    localisation = subcellular.lower()
    rec.outward = (
        gpi
        or _OUTWARD_NOTE in notes
        or any(term in localisation for term in _PLASMA_MEMBRANE_TERMS)
    )

    if rec.is_surface:
        if gpi:
            rec.membrane_class = MembraneClass.GPI_ANCHORED
        elif transmem_count == 1:
            rec.membrane_class = MembraneClass.SINGLE_PASS
        else:
            rec.membrane_class = MembraneClass.MULTI_PASS

    return rec


class UniProtSource(DataSource):
    name = "UniProt"
    namespace = "uniprot"

    def cache_entries(self) -> Iterable[CacheEntry]:
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
        params = {
            "query": QUERY,
            "fields": ",".join(FIELDS),
            "format": "tsv",
            "size": PAGE_SIZE,
        }
        return f"{BASE}?{urllib.parse.urlencode(params)}"

    def fetch(self) -> Path:
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            print("  fetching reviewed human proteome", flush=True)
            return stream_paginated_to_file(
                self._url(), tmp, progress_label="proteome"
            )

        return self.cache.ensure(entry, fetcher)

    def load(self) -> list[ProteinRecord]:
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
    surface = [r for r in records if r.is_surface]
    attached = [r for r in records if r.attached]
    withheld = [r for r in attached if not r.outward]

    # Anchored, not shown to face outward, and placed in a named compartment.
    internal = [r for r in withheld if names_compartment(r.subcellular)]
    # Anchored, but nothing on record says where. Held out for want of evidence,
    # reported separately rather than discarded or counted against them.
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
    }


def load_surface() -> tuple[list[ProteinRecord], dict]:
    records = UniProtSource().load()
    return [r for r in records if r.is_surface], summarise(records)


if __name__ == "__main__":
    recs = UniProtSource().load()
    stats = summarise(recs)
    expected = {
        "entries": 20431,
        "surface": 3496,
        "single_pass": 1464,
        "multi_pass": 1894,
        "gpi_anchored": 138,
        "internal_anchored": 1322,
        "compartment_unresolved": 550,
    }
    for k, v in stats.items():
        exp = expected[k]
        print(f"  {'ok  ' if v == exp else 'DIFF'}  {k}: {v:,}  expected {exp:,}")
