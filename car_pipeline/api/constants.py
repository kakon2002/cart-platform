"""Three endpoints that return a refusal, and the reason for it.

Structural evaluation, functional prediction and model calibration do not
exist. The client's own specification gives each of them a refusal as its
example response, so an endpoint returning exactly that is complete rather than
stubbed: it says the true thing, and it says it in the shape the caller expects.

Each carries its reason beside its status. A reader shown a bare UNKNOWN learns
that the platform is broken; a reader shown UNKNOWN beside *no coordinate
source is connected* learns what would have to change. The second is the whole
point of returning it at all.

These are constants by design. Criterion D8 asserts each returns its declared
status, so the day one of them returns a computed value, something deliberate
has happened rather than something accidental.
"""

from __future__ import annotations

STRUCTURE_UNKNOWN = "UNKNOWN"
FUNCTION_INSUFFICIENT = "INSUFFICIENT_VALIDATED_DATA"
CALIBRATION_BLOCKED = "BLOCKED_INSUFFICIENT_DATA"

# The status each endpoint returns, and nothing else may. Read by D8.
DECLARED = {
    "structure/evaluate": STRUCTURE_UNKNOWN,
    "function/predict": FUNCTION_INSUFFICIENT,
    "learning/calibrate": CALIBRATION_BLOCKED,
}


def structure_evaluate(construct_id: str | None = None) -> dict:
    """Structural evaluation, which has no coordinates to evaluate."""
    return {
        "construct_id": construct_id,
        "structure_status": STRUCTURE_UNKNOWN,
        "structure_score": None,
        "binder_antigen_geometry": None,
        "structural_confidence": None,
        "reasons": [
            "No coordinate file exists anywhere in the connected cache. The "
            "structure source returns entry identifiers, chain labels and "
            "experimental methods; it never returns atoms.",
            "The cached anti-tag binder records its antigen entity as "
            "excluded, so what is held is an antibody with no antigen. There "
            "is no binder-antigen pair to compute geometry between.",
            "Retrieved geometry over deposited complexes is buildable from one "
            "further endpoint and needs no model. Predicted geometry for a "
            "designed receptor is not: no folding model fits the available "
            "hardware.",
            "UNKNOWN is returned rather than 0.00. A zero would rank a "
            "candidate as though its structure had been examined and found "
            "poor, which is a different statement from not having looked.",
        ],
        "what_would_change_this": [
            "A shipping design carrying a target-specific binder with a "
            "deposited complex. Retrieval then becomes worthwhile immediately.",
            "For predicted geometry rather than retrieved: external compute. "
            "That is a purchasing decision and should follow a binder rather "
            "than precede one.",
        ],
    }


def function_predict(construct_id: str | None = None) -> dict:
    """Functional prediction, which has no validated model behind it."""
    return {
        "construct_id": construct_id,
        "status": FUNCTION_INSUFFICIENT,
        "activation_prediction": None,
        "cytotoxic_potential": None,
        "cytokine_profile": None,
        "persistence": None,
        "exhaustion": None,
        "tonic_signalling": None,
        "applicability_domain": None,
        "reasons": [
            "The required training inputs are partner-generated experimental "
            "measurements, and none is connected. This is the one absent stage "
            "that a decision alone cannot close.",
            "No functional quantity is estimated anywhere in the platform, so "
            "there is no proxy to report either.",
            "Every field is null rather than 0.00, because a predicted "
            "activation of zero and an unmeasured activation are different "
            "claims and only one of them is true here.",
        ],
        "what_would_change_this": [
            "Candidate-linked experimental outcomes ingested through the "
            "results endpoint, in enough quantity to train and to state an "
            "applicability domain.",
        ],
    }


def learning_calibrate(dataset_version: str | None = None,
                       model_family: str | None = None) -> dict:
    """Calibration, which has nothing to calibrate against."""
    return {
        "dataset_version": dataset_version,
        "model_family": model_family,
        "status": CALIBRATION_BLOCKED,
        "required_action": "ADD_EXPERIMENTAL_RESULTS",
        "experimental_results_ingested": 0,
        "calibrated_model_version": None,
        "reasons": [
            "No experimental result has ever been ingested, so there is "
            "nothing to calibrate against and no error to measure.",
            "Calibrating against the platform's own predictions would fit the "
            "model to itself. The learning loop closes on measured outcomes or "
            "it does not close.",
        ],
        "what_would_change_this": [
            "Experimental outcomes for candidates this platform proposed, "
            "which is the feedback path the reference document's stage 13 and "
            "14 describe and which begins with the results endpoint.",
        ],
    }
