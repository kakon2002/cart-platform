"""Stage 11 — multi-objective ranking.

Implements `specs/stage11-ranking.md`. Combines numbers the earlier stages made;
measures nothing new.

**No weighted sum.** The objectives are not commensurable and are not made so.
Weights would have to be invented, and a weighted total would let a design with
an unmanageable safety margin be rescued by a strong tumour score — the failure
Stage 3 refuses by keeping its three numbers apart.

The attrition chain is the primary output. For this indication it is the only
part with content.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

NO_DESIGN_REACHES_THE_END = "NO_DESIGN_REACHES_THE_END"
RANKED = "RANKED"

#: Gates in pipeline order. Every pool member is attributed to the first it fails,
#: so the counts partition the pool rather than overlapping.
GATES = (
    "blocked on normal tissue risk",
    "no design recommended",
    "no binder retrieved",
    "no construct assembled",
    "construct over budget",
)


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
    #: Negated so that every objective points the same way: higher is better.
    cleanliness: int = 0
    on_front: bool = False

    @property
    def objectives(self) -> tuple[float, float, float, float] | None:
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


def rank(
    decisions: list[dict],
    binders: dict,
    constructs: dict,
    gated: dict,
    liabilities: dict,
    composites: dict,
    ceiling: float,
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

        # First gate failed, in pipeline order.
        if safety is None or safety.verdict == "BLOCKED":
            entry.failed_at = GATES[0]
        elif row["outcome"] not in ("SINGLE", "DUAL"):
            entry.failed_at = GATES[1]
        elif not binder or not (binder.sequence or binder.structure):
            entry.failed_at = GATES[2]
        elif construct is None or not construct.amino_acid_sequence:
            entry.failed_at = GATES[3]
        elif construct.verdict != "BUILDABLE":
            entry.failed_at = GATES[4]
        else:
            entry.survived = True

        if entry.failed_at:
            attrition[entry.failed_at] += 1
        rows.append(entry)

    survivors = [r for r in rows if r.survived]
    if survivors:
        points = [r.objectives for r in survivors]
        for i in pareto_front(points):
            survivors[i].on_front = True
        status = RANKED
    else:
        status = NO_DESIGN_REACHES_THE_END
    return rows, attrition, status


def configuration_hash(stage9_hash: str, genes: list[str]) -> str:
    payload = {"stage9": stage9_hash, "genes": genes, "gates": list(GATES)}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
