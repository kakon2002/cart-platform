"""Stage 9 — safety gate.

Implements `specs/stage9-safety-gate.md`. Aggregates what earlier stages measured
and adds two readings of its own. It does not re-decide anything Stage 3 decided:
off-tumour risk is carried, never recomputed, because a second implementation of
the tissue mapping would be a second place for its bugs to live.

Passing this gate is not a safety claim. It means three specific questions failed
to show a problem, which is why the verdict is named `PASSES_STATED_CHECKS` and
not `SAFE`.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field

BLOCKED = "BLOCKED"
FLAGGED = "FLAGGED"
NO_GATE = "NO_GATE"
PASSES = "PASSES_STATED_CHECKS"

NOT_CONNECTED = "NOT_CONNECTED"
ORIGIN_UNKNOWN = "ORIGIN_UNKNOWN"

#: Source-species stems in international non-proprietary names, longest first so
#: that `-xizu-` and `-xi-` cannot be confused. A CONVENTION, not a measurement:
#: a molecule re-engineered after naming keeps its original stem, and a binder
#: with no INN has no stem at all.
ORIGIN_STEMS = [
    ("xizu", "chimeric/humanised"),
    ("xi", "chimeric"),
    ("zu", "humanised"),
    ("o", "murine"),
    ("u", "human"),
]
#: Origins that carry an anti-CAR immunogenicity risk worth flagging.
FOREIGN_ORIGINS = {"chimeric", "murine", "chimeric/humanised"}

_SUFFIX = re.compile(r"(.*)(mab|tug|tamig|tamab|tocel)$", re.IGNORECASE)


def binder_origin(name: str) -> str:
    """Source species implied by the name's stem, or ORIGIN_UNKNOWN.

    Read from the two or three letters before the terminal `-mab`. Never guessed
    from sequence, and never asserted where no stem is recognisable.
    """
    if not name:
        return ORIGIN_UNKNOWN
    match = _SUFFIX.match(name.strip())
    if not match:
        return ORIGIN_UNKNOWN
    body = match.group(1).lower()
    for stem, origin in ORIGIN_STEMS:
        if body.endswith(stem):
            return origin
    return ORIGIN_UNKNOWN


@dataclass
class SafetyRecord:
    gene: str
    accession: str
    pool_index: int
    verdict: str
    risk: float | None = None
    risk_organ: str | None = None
    ceiling: float = 0.15
    binder_name: str = ""
    binder_origin: str = ORIGIN_UNKNOWN
    #: Every distinct origin among this target's named binders, because a target
    #: with both a human and a murine binder is not described by either alone.
    binder_origins: list[str] = field(default_factory=list)
    #: Not a lookup. See the specification: answering it needs a k-mer scan of the
    #: variable region against the bulk epitope table, which this stage does not do.
    epitope_immunogenicity: str = NOT_CONNECTED
    trials_total: int = 0
    trials_stopped: int = 0
    trials_stopped_ids: list[str] = field(default_factory=list)
    #: True when the tallies cover fewer studies than the registry holds.
    trials_truncated: bool = False
    reasons: list[str] = field(default_factory=list)


def gate(
    decisions: list[dict],
    binders: dict,
    risks: dict[str, tuple[float | None, str | None]],
    trials: dict,
    ceiling: float = 0.15,
) -> list[SafetyRecord]:
    """One record per pool member, in the order Stage 4 emitted them."""
    out: list[SafetyRecord] = []
    for row in decisions:
        gene = row["gene"]
        risk, organ = risks.get(gene, (None, None))
        record = binders.get(gene)
        # Both routes. Stage 5 defines a structure-route binder as a binder, and
        # counting only sequences here would contradict it — a target with a
        # solved complex and no named therapeutic would read as ungateable.
        usable = ([c for c in record.sequence if c.name]
                  + [c for c in record.structure if c.identifier]) if record else []
        named = [c for c in usable if getattr(c, "route", "") == "sequence"]
        summary = trials.get(gene)

        entry = SafetyRecord(
            gene=gene,
            accession=row["accession"],
            pool_index=row["pool_index"],
            verdict=NO_GATE,
            risk=risk,
            risk_organ=organ,
            ceiling=ceiling,
            trials_total=summary.total if summary else 0,
            trials_stopped=summary.stopped if summary else 0,
            trials_stopped_ids=list(summary.stopped_ids) if summary else [],
            trials_truncated=bool(summary.truncated) if summary else False,
        )

        # Every named binder, not the alphabetically first. Picking one would
        # have read CLDN18 as ORIGIN_UNKNOWN from Ciletatug while Zolbetuximab,
        # the spec's own example, sits in the same set and is chimeric.
        origins = sorted({binder_origin(c.name) for c in named})
        entry.binder_origins = origins
        foreign = [o for o in origins if o in FOREIGN_ORIGINS]
        if named:
            entry.binder_name = ", ".join(sorted(c.name for c in named)[:3])
            # The conservative reading: if any binder for this target is foreign,
            # the target carries that risk, because the design may use it.
            entry.binder_origin = foreign[0] if foreign else (
                origins[0] if origins else ORIGIN_UNKNOWN)

        # Carried from Stage 3, never recomputed. Undefined risk blocks: Stage 3
        # §"Undefined risk is not low risk" makes that explicit, and treating a
        # missing measurement as a pass would invert it.
        if risk is None:
            entry.verdict = BLOCKED
            entry.reasons.append(
                "Stage 3 records no risk for this target; undefined risk is not "
                "low risk"
            )
            out.append(entry)
            continue
        # The ceiling this target was actually routed against, not the
        # persistent one. Stage 4a selects an architecture from the risk profile
        # and the ceiling follows from it; gating every target on the persistent
        # ceiling here would re-apply the blind gate that routing exists to
        # replace, and would block an adaptor design for the very risk its
        # architecture was chosen to carry.
        applied = row.get("route_ceiling") or ceiling
        exposure = row.get("route_exposure") or "persistent"
        entry.ceiling = applied
        if risk > applied:
            entry.verdict = BLOCKED
            entry.reasons.append(
                f"Stage 3 risk {risk:.4f} exceeds the {exposure} ceiling "
                f"{applied} (peak organ {organ})"
            )
            out.append(entry)
            continue
        if exposure == "terminable":
            # Recorded, not silent. This target is admitted on a tolerance for
            # an exposure that can be stopped, which is a different claim from
            # clearing the persistent ceiling, and a reader must see which.
            entry.reasons.append(
                f"admitted against the terminable ceiling {applied}, not the "
                f"persistent {ceiling}: activation requires a separately dosed "
                "adaptor, so the exposure is stoppable"
            )

        if not usable:
            entry.verdict = NO_GATE
            entry.reasons.append("no binder, so there is nothing to gate")
            out.append(entry)
            continue

        if foreign:
            entry.reasons.append(
                f"binder set includes {', '.join(foreign)} by name stem, which is "
                "a convention and not a measurement"
            )
        if entry.trials_stopped:
            floor = " (a floor: tallied over one page)" if entry.trials_truncated else ""
            entry.reasons.append(
                f"{entry.trials_stopped} trial(s) mentioning this symbol were "
                f"terminated, withdrawn or suspended{floor}: "
                f"{', '.join(entry.trials_stopped_ids[:3])}"
            )
        entry.verdict = FLAGGED if entry.reasons else PASSES
        out.append(entry)
    return out


def configuration_hash(
    stage6_hash: str, genes: list[str], ceiling: float
) -> str:
    payload = {
        "stage6": stage6_hash,
        "genes": genes,
        "origin_stems": ORIGIN_STEMS,
        # The ceiling decides every BLOCKED verdict, so a run at a different
        # tolerance must not hash the same as this one.
        "ceiling": ceiling,
        "foreign_origins": sorted(FOREIGN_ORIGINS),
        "epitope": NOT_CONNECTED,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
