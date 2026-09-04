"""Stage 9 — safety gate."""

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


ORIGIN_STEMS = [
    ("xizu", "chimeric/humanised"),
    ("xi", "chimeric"),
    ("zu", "humanised"),
    ("o", "murine"),
    ("u", "human"),
]

FOREIGN_ORIGINS = {"chimeric", "murine", "chimeric/humanised"}

_SUFFIX = re.compile(r"(.*)(mab|tug|tamig|tamab|tocel)$", re.IGNORECASE)


def binder_origin(name: str) -> str:
    """Source species implied by the name's stem, or ORIGIN_UNKNOWN."""
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
    binder_structure_accession: str = ""
    binder_source_organism: str = ""

    binder_origins: list[str] = field(default_factory=list)

    epitope_immunogenicity: str = NOT_CONNECTED
    trials_total: int = 0
    trials_stopped: int = 0
    trials_stopped_ids: list[str] = field(default_factory=list)

    trials_truncated: bool = False
    reasons: list[str] = field(default_factory=list)

    construct_safety: dict | None = None


def structure_binders(constructs) -> dict[str, tuple[str, str, str]]:
    """Per gene, the structure-derived binder a construct carries, if any."""
    from car_pipeline.data.domains import STRUCTURE

    out: dict[str, tuple[str, str, str]] = {}
    for construct in constructs or []:
        for segment in getattr(construct, "segments", []):
            if segment.provenance != STRUCTURE:
                continue
            try:
                from car_pipeline.data.antitag import origin
                organism, verdict = origin()
            except Exception:
                organism, verdict = "not established", "non-human"
            out[construct.gene] = (segment.name, segment.accession or "",
                                   f"{organism}|{verdict}")
            break
    return out


def gate(
    decisions: list[dict],
    binders: dict,
    risks: dict[str, tuple[float | None, str | None]],
    trials: dict,
    ceiling: float = 0.15,
    constructs=None,
) -> list[SafetyRecord]:
    """One record per pool member, in the order Stage 4 emitted them."""
    from car_pipeline.stages import construct_safety as cs

    out: list[SafetyRecord] = []
    from_structure = structure_binders(constructs)
    assembled = {c.gene: c for c in (constructs or []) if c.amino_acid_sequence}
    for row in decisions:
        gene = row["gene"]
        risk, organ = risks.get(gene, (None, None))
        record = binders.get(gene)

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

        built = assembled.get(gene)
        if built is not None:
            entry.construct_safety = cs.analyse(
                built.amino_acid_sequence, built.dna, built.segments)

        origins = sorted({binder_origin(c.name) for c in named})
        entry.binder_origins = origins
        foreign = [o for o in origins if o in FOREIGN_ORIGINS]

        structural = from_structure.get(gene)
        if structural:
            part_name, accession, detail = structural
            organism, verdict = detail.split("|", 1)
            entry.binder_name = part_name
            entry.binder_structure_accession = accession
            entry.binder_source_organism = organism
            entry.binder_origin = verdict
            entry.binder_origins = sorted(set(origins) | {verdict})
        if named:
            entry.binder_name = ", ".join(sorted(c.name for c in named)[:3])

            entry.binder_origin = foreign[0] if foreign else (
                origins[0] if origins else ORIGIN_UNKNOWN)

        if risk is None:
            entry.verdict = BLOCKED
            entry.reasons.append(
                "Stage 3 records no risk for this target; undefined risk is not "
                "low risk"
            )
            out.append(entry)
            continue

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
            entry.reasons.append(
                f"admitted against the terminable ceiling {applied}, not the "
                f"persistent {ceiling}: activation requires a separately dosed "
                "adaptor, so the exposure is stoppable"
            )

        if not usable and not structural:
            entry.verdict = NO_GATE
            entry.reasons.append("no binder, so there is nothing to gate")
            out.append(entry)
            continue

        if structural:
            entry.reasons.append(
                f"the receptor carries a structure-derived binder, "
                f"{entry.binder_source_organism} as deposited in "
                f"{entry.binder_structure_accession}, treated as "
                f"{entry.binder_origin} because no humanised sequence is "
                "established for it. This is read from the deposition, not "
                "from a name stem, which a structure-derived binder does not "
                "carry"
            )
            entry.reasons.append(
                "epitope-level immunogenicity remains NOT_CONNECTED for it: "
                "no epitope source is connected, so the species is known and "
                "the immunogenicity is not"
            )

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
    """Fingerprint the safety-gate configuration, ceiling included."""
    payload = {
        "stage6": stage6_hash,
        "genes": genes,
        "origin_stems": ORIGIN_STEMS,
        "ceiling": ceiling,
        "foreign_origins": sorted(FOREIGN_ORIGINS),
        "epitope": NOT_CONNECTED,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
