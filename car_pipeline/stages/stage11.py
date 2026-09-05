"""Stage 11 — multi-objective ranking."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from car_pipeline.stages import scoring

NO_DESIGN_REACHES_THE_END = "NO_DESIGN_REACHES_THE_END"
RANKED = "RANKED"


RANKED_AWAITING_BINDER = "RANKED_AWAITING_BINDER"


GATES = (
    "blocked on normal tissue risk",
    "no design recommended",
    "no binder retrieved",
    "no construct assembled",
    "construct over budget",
)


RECOMMENDED = ("SINGLE", "DUAL", "ADAPTOR")


PASSED_ALL_GATES = "PASSED_ALL_GATES"

# A machine token per gate, because the GATES strings are prose written to be
# read in an attrition table and a column a frontend renders should not be one.
GATE_STATUS = {
    GATES[0]: "BLOCKED_ON_NORMAL_TISSUE_RISK",
    GATES[1]: "NO_DESIGN_RECOMMENDED",
    GATES[2]: "NO_BINDER_RETRIEVED",
    GATES[3]: "NO_CONSTRUCT_ASSEMBLED",
    GATES[4]: "OVER_PAYLOAD_BUDGET",
}

ADVANCE = "ADVANCE"
BACKUP = "BACKUP"
VALIDATE = "VALIDATE"
REQUIRES_EVIDENCE = "REQUIRES_EVIDENCE"
EXCLUDED = "EXCLUDED"

DECISIONS = (ADVANCE, BACKUP, VALIDATE, REQUIRES_EVIDENCE, EXCLUDED)

# Which failures a measurement could later clear, and which it could not. A
# target over the safety ceiling is not waiting on evidence; a target with no
# binder retrieved is waiting on exactly that.
GATE_DECISION = {
    GATES[0]: EXCLUDED,
    GATES[1]: EXCLUDED,
    GATES[2]: REQUIRES_EVIDENCE,
    GATES[3]: REQUIRES_EVIDENCE,
    GATES[4]: EXCLUDED,
}

# Position is Stage 4's composite ordering, inherited and carried through
# unchanged. It is recorded here so a reader is told where the order came from
# rather than inferring that Stage 11 produced it. Nothing decisional hangs off
# it: the decision column reads front membership, which Stage 11 does compute.
POSITION_BASIS = "stage4 composite order, inherited; not a ranking Stage 11 computed"


@dataclass
class Ranked:
    gene: str
    accession: str
    pool_index: int
    survived: bool
    failed_at: str | None
    attractiveness: float | None = None
    safety_margin: float | None = None
    binder_count: int = 0

    binder_supplied: bool = True

    cleanliness: int = 0
    on_front: bool = False

    position: int | None = None
    candidate_id: str | None = None
    gate_status: str = ""
    decision: str = ""

    # Level B. None on every candidate that failed a gate, because scoring is
    # only ever reached by survivors -- criterion W4.
    scorecard: object | None = None
    overall: float | None = None

    @property
    def objectives(self) -> tuple[float, float, float, float] | None:
        """The objectives this candidate is compared on, or None if it did not survive."""
        if not self.survived:
            return None
        return (
            self.attractiveness or 0.0,
            self.safety_margin or 0.0,
            float(self.binder_count),
            float(self.cleanliness),
        )


def dominates(a: tuple, b: tuple) -> bool:
    """True when a is at least as good everywhere and better somewhere."""
    return all(x >= y for x, y in zip(a, b)) and any(x > y for x, y in zip(a, b))


def pareto_front(points: list[tuple]) -> list[int]:
    """Indices of the non-dominated points."""
    front = []
    for i, p in enumerate(points):
        if not any(dominates(q, p) for j, q in enumerate(points) if j != i):
            front.append(i)
    return front


def decision_for(entry: Ranked, scored: bool | None) -> str:
    """What should happen to this candidate, from its gate status and the front.

    `scored` is tri-state on purpose. True means a score was emitted, False
    means the candidate cleared every gate but too little of it was measured to
    score, and None means no scoring stage has run at all. Only False produces
    VALIDATE, so that branch is unreachable from the pipeline until scoring
    exists; it is exercised directly rather than left untested until then.
    """
    if not entry.survived:
        return GATE_DECISION[entry.failed_at]
    if scored is False:
        return VALIDATE
    return ADVANCE if entry.on_front else BACKUP


def candidate_id(indication_key: str, position: int) -> str:
    """A within-run label. The durable identity is the gene and the hash chain."""
    return f"CAR-{indication_key.upper()}-{position:03d}"


def rank(
    decisions: list[dict],
    binders: dict,
    constructs: dict,
    gated: dict,
    liabilities: dict,
    composites: dict,
    ceiling: float,
    *,
    indication_key: str,
    stage3_rows: dict,
    budget_bp: int,
) -> tuple[list[Ranked], dict[str, int], str]:
    """Attribute every pool member to its first failed gate, then rank survivors."""
    rows: list[Ranked] = []
    attrition: dict[str, int] = {g: 0 for g in GATES}

    for row in decisions:
        gene = row["gene"]
        safety = gated.get(gene)
        binder = binders.get(gene)
        construct = constructs.get(gene)
        entry = Ranked(
            gene=gene, accession=row["accession"], pool_index=row["pool_index"],
            survived=False, failed_at=None,
            attractiveness=composites.get(gene),
            binder_count=(len(binder.sequence) + len(binder.structure)) if binder else 0,
            cleanliness=-min(
                (l.flag_count for l in liabilities.get(gene, [])), default=0),
        )
        if safety is not None and safety.risk is not None:
            entry.safety_margin = round(ceiling - safety.risk, 4)

        if safety is None or safety.verdict == "BLOCKED":
            entry.failed_at = GATES[0]
        elif row["outcome"] not in RECOMMENDED:
            entry.failed_at = GATES[1]
        elif (row["outcome"] != "ADAPTOR"
              and (not binder or not (binder.sequence or binder.structure))):
            entry.failed_at = GATES[2]
        elif construct is None or not construct.segments:
            entry.failed_at = GATES[3]
        elif construct.verdict != "BUILDABLE":
            entry.failed_at = GATES[4]
        else:
            entry.survived = True
            entry.binder_supplied = bool(construct.amino_acid_sequence)

        if entry.failed_at:
            attrition[entry.failed_at] += 1
        rows.append(entry)

    survivors = [r for r in rows if r.survived]
    if survivors:
        points = [r.objectives for r in survivors]
        for i in pareto_front(points):
            survivors[i].on_front = True

        for position, entry in enumerate(survivors, 1):
            entry.position = position
            entry.candidate_id = candidate_id(indication_key, position)

        # Level B, survivors only. A candidate that failed a gate never reaches
        # this loop, which is what makes the weighted sum safe: no weight can
        # rescue a gate failure because no gate failure is scored.
        for entry in survivors:
            entry.scorecard = scoring.score(
                entry,
                stage3_rows.get(entry.gene),
                gated.get(entry.gene),
                constructs.get(entry.gene),
                binders.get(entry.gene),
                budget_bp,
            )
            entry.overall = entry.scorecard.overall

        status = (RANKED if any(r.binder_supplied for r in survivors)
                  else RANKED_AWAITING_BINDER)
    else:
        status = NO_DESIGN_REACHES_THE_END

    # After the front and the scoring, because the decision reads both.
    for entry in rows:
        entry.gate_status = (PASSED_ALL_GATES if entry.survived
                             else GATE_STATUS[entry.failed_at])
        entry.decision = decision_for(
            entry, entry.scorecard.scored if entry.scorecard else None)

    return rows, attrition, status


def configuration_hash(stage9_hash: str, genes: list[str]) -> str:
    """Fingerprint the ranking configuration, gates and recommendations included."""
    payload = {"stage9": stage9_hash, "genes": genes, "gates": list(GATES),
               "recommended": list(RECOMMENDED),
               "decisions": list(DECISIONS),
               "gate_status": [GATE_STATUS[g] for g in GATES],
               "scoring": scoring.configuration_hash()}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
