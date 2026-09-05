"""Stage 11 Level B — the declared scoring frame.

Level A is the gate chain in `stage11`. Nothing here can rescue a candidate that
failed one: `score` is only ever called on survivors, and criterion W4 asserts
that no gate failure carries a score. What this differentiates is residual
margin among candidates that have already passed, which is the thing the
no-weighted-sum rule was never trying to prevent.

Every coefficient is declared before any run and versioned. Nothing here is
fitted to an observation.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

MEASURED = "MEASURED"
UNKNOWN = "UNKNOWN"
NOT_APPLICABLE = "NOT_APPLICABLE"

STATES = (MEASURED, UNKNOWN, NOT_APPLICABLE)


# Bump when any coefficient below changes. It enters the Stage 11 configuration
# hash, so a run under different weights cannot compare equal to one under
# these.
WEIGHT_VERSION = "wm-scoring-1"


# The nine weighted components, in the order the ordering rule puts them:
# measured outranks modelled, and decisiveness outranks granularity.
WEIGHTS: dict[str, float] = {
    "tumour_coverage": 0.18,
    "malignant_specificity": 0.16,
    "normal_tissue_safety": 0.16,
    "binder_quality": 0.12,
    "manufacturability": 0.10,
    "developability": 0.10,
    "structural_feasibility": 0.08,
    "functional_prediction": 0.06,
    "pairing_robustness": 0.04,
}

# The two components that parameterise the adjustment rather than the sum.
# Keeping them outside is what stops evidence confidence being combined with
# normal-tissue risk, which is a standing rule this specification does not
# overturn.
CONFIDENCE_EXPONENT = 1.0
UNCERTAINTY_PENALTY = 0.5

ADJUSTMENT_COEFFICIENTS = {
    "evidence_confidence_exponent": CONFIDENCE_EXPONENT,
    "prediction_uncertainty_penalty": UNCERTAINTY_PENALTY,
}

# All eleven the reference document names. Nine carry weights, two are the
# adjustment. None is dropped.
COMPONENTS = tuple(WEIGHTS) + ("evidence_confidence", "prediction_uncertainty")

MINIMUM_SCORED_FRACTION = 0.50


@dataclass(frozen=True)
class Component:
    """One component on one candidate, in exactly one of the three states."""

    key: str
    state: str
    value: float | None
    source: str

    @property
    def weight(self) -> float:
        """The declared weight, which exists whatever state the component is in."""
        return WEIGHTS[self.key]

    def as_payload(self) -> dict:
        return {"component": self.key, "state": self.state, "value": self.value,
                "weight": self.weight, "source": self.source}


def measured(key: str, value: float, source: str) -> Component:
    """A component with a value in [0, 1] and the stage that produced it named."""
    return Component(key, MEASURED, _clamp(value), source)


def unknown(key: str, reason: str) -> Component:
    """No measurement exists. The reason is carried so a reader sees which."""
    return Component(key, UNKNOWN, None, reason)


def not_applicable(key: str, reason: str) -> Component:
    """The question does not arise for this design. Not the same as missing."""
    return Component(key, NOT_APPLICABLE, None, reason)


def _clamp(x: float) -> float:
    """Components are normalised to [0, 1]; nothing outside it is a score."""
    return 0.0 if x < 0.0 else 1.0 if x > 1.0 else float(x)


@dataclass
class Scorecard:
    """Everything needed to recompute the score, carried on the candidate.

    W3 requires the record to reconstruct its own score, so the component
    values, the weights, the adjustment and its inputs are all here rather than
    only the total.
    """

    gene: str
    components: dict[str, Component]
    applicable: float
    measured_weight: float
    fraction: float
    weighted_sum: float | None
    evidence_confidence: float | None
    prediction_uncertainty: float | None
    confidence_adjustment: float | None
    overall: float | None
    scored: bool
    weight_version: str = WEIGHT_VERSION
    reasons: list[str] = field(default_factory=list)

    @property
    def unknown_components(self) -> list[str]:
        return [k for k, c in self.components.items() if c.state == UNKNOWN]

    @property
    def not_applicable_components(self) -> list[str]:
        return [k for k, c in self.components.items() if c.state == NOT_APPLICABLE]

    def as_payload(self) -> dict:
        """The scorecard as the package and the API carry it."""
        return {
            "weight_version": self.weight_version,
            "components": [c.as_payload() for c in self.components.values()],
            "applicable_weight": round(self.applicable, 6),
            "measured_weight": round(self.measured_weight, 6),
            "scored_fraction": round(self.fraction, 6),
            "minimum_scored_fraction": MINIMUM_SCORED_FRACTION,
            "weighted_sum": (None if self.weighted_sum is None
                             else round(self.weighted_sum, 6)),
            "evidence_confidence": self.evidence_confidence,
            "prediction_uncertainty": self.prediction_uncertainty,
            "confidence_adjustment": (None if self.confidence_adjustment is None
                                      else round(self.confidence_adjustment, 6)),
            "overall": None if self.overall is None else round(self.overall, 6),
            "scored": self.scored,
            "unknown_components": self.unknown_components,
            "not_applicable_components": self.not_applicable_components,
            "reasons": self.reasons,
        }


def confidence_adjustment(
    evidence_confidence: float | None,
    prediction_uncertainty: float | None,
) -> tuple[float | None, str]:
    """The multiplier, and a statement of which of its two terms was applied.

    Where prediction uncertainty is UNKNOWN the penalty term is not applied and
    its absence is said so. Treating it as zero would be the favourable
    imputation the document forbids: it would make an unmeasured candidate
    indistinguishable from one measured to have no uncertainty.
    """
    if evidence_confidence is None:
        return None, "no evidence confidence, so no adjustment could be formed"
    base = float(evidence_confidence) ** CONFIDENCE_EXPONENT
    if prediction_uncertainty is None:
        return base, (
            "evidence confidence applied at exponent "
            f"{CONFIDENCE_EXPONENT}; the prediction-uncertainty penalty is NOT "
            "applied because no uncertainty is measured. It is not treated as "
            "zero uncertainty, which would flatter an unmeasured candidate.")
    penalty = 1.0 - UNCERTAINTY_PENALTY * float(prediction_uncertainty)
    return base * penalty, (
        f"evidence confidence at exponent {CONFIDENCE_EXPONENT}, times a "
        f"{UNCERTAINTY_PENALTY} uncertainty penalty on "
        f"{prediction_uncertainty}")


def combine(gene: str, components: list[Component],
            evidence_confidence: float | None,
            prediction_uncertainty: float | None) -> Scorecard:
    """Normalise over the measured subset and apply the adjustment.

        applicable = sum of weights not NOT_APPLICABLE
        measured   = sum of weights MEASURED
        fraction   = measured / applicable
        overall    = (sum(w*c over MEASURED) / measured) * adjustment

    The denominator is the measured weight and never 1. Dividing by 1 would
    hand the measured components the weight belonging to all nine, and a
    candidate scored on four would be directly comparable to one scored on
    nine. It would not be.
    """
    by_key = {c.key: c for c in components}
    missing = sorted(set(WEIGHTS) - set(by_key))
    if missing:
        raise ValueError(
            f"{gene}: no component supplied for {missing}. Every weighted "
            "component must be present in one of the three states; a component "
            "left out is an imputation by omission.")

    applicable = sum(c.weight for c in components if c.state != NOT_APPLICABLE)
    measured_weight = sum(c.weight for c in components if c.state == MEASURED)
    fraction = (measured_weight / applicable) if applicable else 0.0

    adjustment, adjustment_note = confidence_adjustment(
        evidence_confidence, prediction_uncertainty)

    reasons = [adjustment_note]
    scored = fraction >= MINIMUM_SCORED_FRACTION and adjustment is not None

    if not scored:
        weighted_sum = None
        overall = None
        if fraction < MINIMUM_SCORED_FRACTION:
            reasons.append(
                f"No overall score is emitted: {fraction:.4f} of the "
                f"applicable frame is measured, below the "
                f"{MINIMUM_SCORED_FRACTION} floor. A score resting on less "
                "than half the frame summarises whichever half happened to be "
                "measured, not the frame. The unmeasured components are named "
                "rather than imputed.")
    else:
        weighted_sum = sum(c.weight * c.value
                           for c in components if c.state == MEASURED)
        overall = (weighted_sum / measured_weight) * adjustment
        reasons.append(
            f"Scored on {measured_weight:.4f} of {applicable:.4f} applicable "
            f"weight ({fraction:.4f}). The remaining components are named, not "
            "imputed, and the score is normalised over what was measured.")

    return Scorecard(
        gene=gene,
        components={k: by_key[k] for k in WEIGHTS},
        applicable=applicable,
        measured_weight=measured_weight,
        fraction=fraction,
        weighted_sum=weighted_sum,
        evidence_confidence=evidence_confidence,
        prediction_uncertainty=prediction_uncertainty,
        confidence_adjustment=adjustment,
        overall=overall,
        scored=scored,
        reasons=reasons,
    )


# ---------------------------------------------------------------------------
# Where each component comes from. Every branch names its source or its reason.
# ---------------------------------------------------------------------------

# The reference document's §14.1 names "tumor coverage / target prevalence",
# so coverage reads Stage 3's prevalence component rather than a new quantity.
COVERAGE_SOURCE = "patient_prevalence"
SPECIFICITY_SOURCE = "malignant_vs_stroma"


def _from_stage3(key: str, ranked, source_component: str) -> Component:
    """One Stage 3 component, or UNKNOWN naming why Stage 3 could not measure it."""
    if ranked is None:
        return unknown(key, "no Stage 3 record for this gene")
    entry = ranked.components.get(source_component)
    if entry is None:
        return unknown(key, f"Stage 3 emitted no {source_component} component")
    if entry.value is None:
        return unknown(
            key, f"Stage 3 {source_component} not measured: "
                 f"{getattr(entry, 'reason', None) or 'no reason recorded'}")
    return measured(key, entry.value, f"stage3.{source_component}")


def _safety(safety) -> Component:
    """Residual margin below the ceiling this candidate was actually judged against.

    The applied ceiling is used, not the persistent one. An adaptor design is
    admitted against the terminable ceiling because its exposure is
    terminable, and measuring its margin against the persistent ceiling would
    report a negative margin for a candidate that passed its gate.
    """
    key = "normal_tissue_safety"
    if safety is None:
        return unknown(key, "no safety record for this gene")
    if safety.risk is None:
        return unknown(key, "normal-tissue risk is not measured for this target")
    ceiling = getattr(safety, "ceiling", None)
    if not ceiling:
        return unknown(key, "no ceiling recorded on the safety verdict")
    return measured(key, (ceiling - safety.risk) / ceiling,
                    f"stage9 residual margin below the applied ceiling {ceiling}")


def _manufacturability(construct, budget_bp: int) -> Component:
    """Payload headroom after the hard packaging gate, which is degree not kind."""
    key = "manufacturability"
    if construct is None or not construct.total_bp:
        return unknown(key, "no assembled construct to measure")
    return measured(key, (budget_bp - construct.total_bp) / budget_bp,
                    f"stage6 headroom against the {budget_bp} bp payload budget")


def _binder_quality(construct, binder) -> Component:
    """Binder quality, which no connected source measures for either route."""
    key = "binder_quality"
    if construct is not None and construct.outcome == "ADAPTOR":
        return unknown(
            key, "this is an adaptor design: the receptor binds a tag, no "
                 "anti-tag binder is retrieved, and there is no binder to score")
    if not binder or not (binder.sequence or binder.structure):
        return unknown(key, "no binder was retrieved for this target")
    return unknown(
        key, "a binder is retrieved, but predicted affinity is NOT_CONNECTED "
             "and no other quality measurement exists. Provenance is recorded "
             "and is not a quality score")


def _developability() -> Component:
    """UNKNOWN by standing decision, not for want of data."""
    return unknown(
        "developability",
        "Stage 10 counts sequence liabilities and refuses to sum them into a "
        "score, because one flag fires on every binder in the pool and a sum "
        "would hide it. That decision stands; this specification does not "
        "overturn it, so the component has no value to carry")


def _structural() -> Component:
    return unknown("structural_feasibility",
                   "Stage 7 does not exist. It is buildable and unbuilt")


def _functional() -> Component:
    return unknown("functional_prediction",
                   "Stage 8 does not exist and is not buildable from what is "
                   "connected: the required training inputs are not available")


def _pairing(construct) -> Component:
    """Only multi-target designs have a partner whose stability can be robust."""
    key = "pairing_robustness"
    outcome = getattr(construct, "outcome", None)
    if outcome == "DUAL":
        return unknown(
            key, "this is a multi-target design, so pairing robustness applies, "
                 "and no robustness measurement is connected")
    article = "an" if (outcome or "s")[0].upper() in "AEIOU" else "a"
    return not_applicable(
        key, f"{article} {outcome or 'single-target'} design names no partner, so "
             "pairing robustness is a question that does not arise. This is "
             "not missing evidence")


def build_components(entry, ranked, safety, construct, binder,
                     budget_bp: int) -> list[Component]:
    """Every weighted component for one candidate, each in exactly one state."""
    return [
        _from_stage3("tumour_coverage", ranked, COVERAGE_SOURCE),
        _from_stage3("malignant_specificity", ranked, SPECIFICITY_SOURCE),
        _safety(safety),
        _binder_quality(construct, binder),
        _manufacturability(construct, budget_bp),
        _developability(),
        _structural(),
        _functional(),
        _pairing(construct),
    ]


def score(entry, ranked, safety, construct, binder, budget_bp: int) -> Scorecard:
    """The Level B scorecard for one gate-passing candidate."""
    components = build_components(
        entry, ranked, safety, construct, binder, budget_bp)
    return combine(
        entry.gene, components,
        evidence_confidence=getattr(ranked, "confidence", None),
        prediction_uncertainty=None,
    )


def configuration_hash() -> str:
    """Fingerprint the declared frame, so a run under other weights differs."""
    payload = {
        "version": WEIGHT_VERSION,
        "weights": {k: WEIGHTS[k] for k in sorted(WEIGHTS)},
        "adjustment": ADJUSTMENT_COEFFICIENTS,
        "floor": MINIMUM_SCORED_FRACTION,
        "coverage_source": COVERAGE_SOURCE,
        "specificity_source": SPECIFICITY_SOURCE,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
