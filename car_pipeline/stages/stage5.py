"""Stage 5 — binder discovery.

Implements `specs/stage5-binder-discovery.md`. Retrieval only: this stage finds
and characterises binders that already exist. It does not design one, does not
predict affinity and does not dock anything, because a prediction printed beside
a measurement is read as one.

**Two routes, reported apart and never summed.**

* the **structure route** — a deposited complex, found by accession and confirmed
  as antibody-containing by curated chain annotation rather than by reading entity
  description text. It carries an epitope.
* the **sequence route** — a named therapeutic with its variable-region
  sequences, a clinical stage and a target annotation. It carries no epitope, and
  it carries the thing a construct is actually built from.

A target with a sequence-route binder and no structure is **not** `NO_BINDER`.
Collapsing the two would discard the sequences, which are the more useful half.

**Affinity is `NOT_CONNECTED` for every candidate**, and the reason is measured
rather than assumed: the curated structure summary carried affinity in a previous
release and does not in this one, and the bioactivity database returns zero
records for both validation molecules and both validation targets. "We did not
rank on affinity" and "we could not" are different statements and only the second
is true.

**This stage does not re-rank targets.** Binder availability tracks how much
attention a protein has had. Ordering anything by it would let the literature
choose the targets, which is the thing the discovery stages exist to avoid.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from car_pipeline.data.antibodies import AntibodySource, Therapeutic
from car_pipeline.data.structures import entries_for, entry_summary

NO_BINDER = "NO_BINDER"
STRUCTURE_AND_SEQUENCE = "STRUCTURE_AND_SEQUENCE"
BINDER_STRUCTURE_ONLY = "BINDER_STRUCTURE_ONLY"
BINDER_SEQUENCE_ONLY = "BINDER_SEQUENCE_ONLY"

NOT_CONNECTED = "NOT_CONNECTED"
ISOFORM_UNRESOLVED = "ISOFORM_UNRESOLVED"

#: Coding length of the antigen-recognition domain once converted to a
#: single-chain format, in base pairs: heavy variable, a linker, light variable.
#: Three times the residue count, stated so the arithmetic is visible rather than
#: hidden inside a number.
SCFV_RESIDUES = 250
SCFV_BP = SCFV_RESIDUES * 3


@dataclass
class Candidate:
    route: str
    identifier: str
    name: str = ""
    fmt: str = ""
    clinical_stage: str = ""
    status: str = ""
    heavy_sequence: str = ""
    light_sequence: str = ""
    antigen_chain: str = ""
    antigen_name: str = ""
    method: str = ""
    #: Never a number. See the module docstring.
    affinity: str = NOT_CONNECTED
    #: Neither route can say which isoform a binder engages: the therapeutic
    #: table records the gene, and the structure route would need the deposited
    #: antigen sequence aligned against each isoform. Reported, never guessed.
    isoform: str = ISOFORM_UNRESOLVED

    @property
    def car_bp(self) -> int | None:
        """Size of the CAR-converted binder, where a sequence is available."""
        if self.heavy_sequence and self.light_sequence:
            return (len(self.heavy_sequence) + len(self.light_sequence) + 15) * 3
        if self.route == "structure":
            return SCFV_BP
        return None


@dataclass
class TargetBinders:
    gene: str
    accession: str
    pool_index: int
    outcome: str
    partner: str | None
    entries: list[str] = field(default_factory=list)
    structure: list[Candidate] = field(default_factory=list)
    sequence: list[Candidate] = field(default_factory=list)
    #: Entries that exist for this accession but carry no antibody instance.
    #: Reported because "the protein has structures" and "the protein has a
    #: binder" are different claims and the gap between them is the point.
    entries_without_antibody: int = 0
    #: Entries dropped because their coordinates are computed rather than
    #: measured. Kept apart from the row above: "no antibody in it" and "not an
    #: experiment" are different reasons to discard an entry.
    entries_excluded_as_model: int = 0

    @property
    def verdict(self) -> str:
        if self.structure and self.sequence:
            return STRUCTURE_AND_SEQUENCE
        if self.structure:
            return BINDER_STRUCTURE_ONLY
        if self.sequence:
            return BINDER_SEQUENCE_ONLY
        return NO_BINDER

    @property
    def structure_verdict(self) -> str:
        return BINDER_STRUCTURE_ONLY if self.structure else NO_BINDER

    @property
    def sequence_verdict(self) -> str:
        return BINDER_SEQUENCE_ONLY if self.sequence else NO_BINDER


def _sequence_candidates(therapeutics: list[Therapeutic]) -> list[Candidate]:
    out = []
    for t in sorted(therapeutics, key=lambda x: x.name):
        out.append(
            Candidate(
                route="sequence",
                identifier=t.name,
                name=t.name,
                fmt=t.fmt,
                clinical_stage=t.highest_trial,
                status=t.status,
                heavy_sequence=t.heavy_sequence,
                light_sequence=t.light_sequence,
            )
        )
    return out


def retrieve(
    decisions: list[dict],
    source: AntibodySource | None = None,
    progress: bool = True,
) -> list[TargetBinders]:
    """One record per pool member, in the order Stage 4 emitted them.

    The order is carried unchanged and nothing here reorders it.
    """
    source = source or AntibodySource()
    by_pdb = source.structures()
    by_target = source.therapeutics_by_target()

    out: list[TargetBinders] = []
    for n, row in enumerate(decisions, 1):
        gene = row["gene"]
        record = TargetBinders(
            gene=gene,
            accession=row["accession"],
            pool_index=row["pool_index"],
            outcome=row["outcome"],
            partner=row.get("partner"),
        )
        record.entries = entries_for(record.accession)

        for entry_id in record.entries:
            instances = by_pdb.get(entry_id.lower())
            if not instances:
                record.entries_without_antibody += 1
                continue
            summary = entry_summary(entry_id)
            if summary["is_model"]:
                # A computed model is not retrieved evidence. Counted, so the
                # structures-versus-binders gap reported downstream stays exact.
                record.entries_excluded_as_model += 1
                continue
            for inst in instances:
                record.structure.append(
                    Candidate(
                        route="structure",
                        identifier=f"{entry_id}:{inst.heavy_chain}{inst.light_chain}",
                        name=summary["title"][:90],
                        fmt="Fab" if inst.light_chain else "single domain",
                        antigen_chain=inst.antigen_chain,
                        antigen_name=inst.antigen_name,
                        method=inst.method,
                    )
                )

        record.sequence = _sequence_candidates(by_target.get(gene, []))
        out.append(record)
        if progress and n % 25 == 0:
            print(f"    retrieved {n}/{len(decisions)}", flush=True)
    return out


def configuration_hash(stage4_hash: str, genes: list[str]) -> str:
    payload = {
        "stage4": stage4_hash,
        "genes": genes,
        "routes": ["structure", "sequence"],
        "scfv_bp": SCFV_BP,
        "affinity": NOT_CONNECTED,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
