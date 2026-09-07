"""B10 — two binder rankings and a front, over the same candidates.

The reference document requires binding rank and CAR-suitability rank as
separate outputs and is explicit that highest affinity must not equal best CAR
binder. They are two different *weighted* questions, so both reuse the Level B
machinery unchanged: normalisation over the measured subset, a floor below
which no score is emitted, three states where missing is never a favourable
zero, and a confidence adjustment kept outside the sum.

The front is kept as a third output because it is the strongest guarantee of
the document's own rule. Front membership is computed from component values and
never from weights, so no choice of weights can move it: where the two rankings
disagree, the front says which binders are not beaten on every measured axis at
once.

Level A comes first here as it does everywhere else. A binder whose recorded
antigen names a different protein is excluded rather than scored, so no weight
can rescue it.
"""

from __future__ import annotations

from car_pipeline.stages import binder_check, scoring, stage11

BINDING_VERSION = "binder_binding_v1"
SUITABILITY_VERSION = "binder_suitability_v1"

# Which binder appears to bind the target most strongly? Affinity is the direct
# answer, so it carries the most weight; the other two are evidence that
# binding occurs at all and evidence that the first two can be read.
BINDING_WEIGHTS = {
    "affinity": 0.55,
    "structural_verification": 0.30,
    "evidence_completeness": 0.15,
}

# Which binder makes the best CAR? Affinity is deliberately not the largest
# weight: putting it on top would contradict the document's own thesis in the
# act of implementing it.
SUITABILITY_WEIGHTS = {
    "epitope_accessibility": 0.20,
    "membrane_proximity": 0.15,
    "shedding_suitability": 0.15,
    "specificity": 0.15,
    "affinity": 0.15,
    "developability": 0.12,
    "humanness": 0.08,
}

EXCLUDED = "EXCLUDED"
RANKED = "RANKED"
NOT_SCORED = "NOT_SCORED"

# The metadata a deposited record is expected to carry, so completeness is a
# fraction of something declared rather than of whatever happened to be there.
EXPECTED_METADATA = ("method", "resolution", "antigen_chain", "antigen_type")


def _component(key: str, weights: dict, state: str, value=None, source=""):
    """One component under a named weight set."""
    return {"key": key, "weight": weights[key], "state": state,
            "value": value, "source": source}


def binding_components(row: dict) -> list[dict]:
    """The binding-rank components for one binder."""
    evidence = row.get("structural_evidence") or {}
    verified = evidence.get("structural_status") == "VERIFIED_EXPERIMENTAL"

    present = [f for f in EXPECTED_METADATA if evidence.get(f) is not None]
    completeness = len(present) / len(EXPECTED_METADATA)

    return [
        _component("affinity", BINDING_WEIGHTS, scoring.UNKNOWN, None,
                   "no connected evidence release carries an affinity, KD or "
                   "free-energy column, so no binder has a retrieved value and "
                   "none is predicted"),
        _component(
            "structural_verification", BINDING_WEIGHTS,
            scoring.MEASURED if verified else scoring.UNKNOWN,
            1.0 if verified else None,
            f"a deposited experimental complex whose recorded antigen names "
            f"the target ({evidence.get('entry')})" if verified else
            "no verified deposited complex for this binder"),
        _component(
            "evidence_completeness", BINDING_WEIGHTS,
            scoring.MEASURED if verified else scoring.UNKNOWN,
            completeness if verified else None,
            f"{len(present)} of {len(EXPECTED_METADATA)} expected metadata "
            f"fields recorded" if verified else
            "no deposited record to read metadata from"),
    ]


def suitability_components(row: dict) -> list[dict]:
    """The CAR-suitability components for one binder."""
    reasons = {
        "epitope_accessibility":
            "no epitope is established for this binder; locating one needs "
            "coordinates, and none is connected",
        "membrane_proximity":
            "membrane proximity is a property of the epitope, which is not "
            "established",
        "shedding_suitability":
            "the target's cleavage is known from its annotated chains, but "
            "whether this binder engages a released fragment needs an epitope",
        "specificity":
            "no homolog screen exists: the connected proteome cache carries no "
            "sequence column and no alignment method is installed",
        "affinity":
            "no connected evidence release carries an affinity column",
        "developability":
            "liability flags are counted and listed and never summed into a "
            "score, by a standing decision this ranking does not overturn",
        "humanness":
            "no germline reference is connected; the naming-convention signal "
            "held elsewhere is a convention rather than a measurement",
    }
    return [_component(key, SUITABILITY_WEIGHTS, scoring.UNKNOWN, None, reason)
            for key, reason in reasons.items()]


def _combine(components: list[dict], weights: dict) -> dict:
    """Normalise over the measured subset, using the declared weight set."""
    applicable = sum(c["weight"] for c in components
                     if c["state"] != scoring.NOT_APPLICABLE)
    measured = sum(c["weight"] for c in components
                   if c["state"] == scoring.MEASURED)
    fraction = (measured / applicable) if applicable else 0.0
    scored = fraction >= scoring.MINIMUM_SCORED_FRACTION

    overall = None
    if scored:
        overall = sum(c["weight"] * c["value"] for c in components
                      if c["state"] == scoring.MEASURED) / measured

    return {
        "components": components,
        "applicable_weight": round(applicable, 6),
        "measured_weight": round(measured, 6),
        "scored_fraction": round(fraction, 6),
        "minimum_scored_fraction": scoring.MINIMUM_SCORED_FRACTION,
        "overall": None if overall is None else round(overall, 6),
        "unknown_components": [c["key"] for c in components
                               if c["state"] == scoring.UNKNOWN],
        "reason": (
            f"Scored on {measured:.4f} of {applicable:.4f} applicable weight."
            if scored else
            f"No score emitted: {fraction:.4f} of the applicable frame is "
            f"measured, below the {scoring.MINIMUM_SCORED_FRACTION} floor. The "
            f"unmeasured components are named rather than imputed."),
    }


def rank(rows: list[dict]) -> dict:
    """Both rankings and the front, over binders that clear the hard gate."""
    excluded = [r for r in rows
                if r.get("target_match") == binder_check.FAIL]
    eligible = [r for r in rows
                if r.get("target_match") != binder_check.FAIL]

    scored = []
    for row in eligible:
        binding = _combine(binding_components(row), BINDING_WEIGHTS)
        suitability = _combine(suitability_components(row), SUITABILITY_WEIGHTS)
        scored.append({
            "binder_id": row.get("binder_id"),
            "target_id": row.get("target_id"),
            "identifier": row.get("identifier"),
            "route": row.get("route"),
            "target_match": row.get("target_match"),
            "binding": binding,
            "suitability": suitability,
            "binding_rank": None,
            "car_suitability_rank": None,
        })

    for key, field in (("binding", "binding_rank"),
                       ("suitability", "car_suitability_rank")):
        have = [s for s in scored if s[key]["overall"] is not None]
        for position, entry in enumerate(
                sorted(have, key=lambda s: -s[key]["overall"]), 1):
            entry[field] = position

    # The front, over every component either ranking measured. Computed from
    # values and never from weights, so no weight set can move its membership.
    axes: list[str] = []
    for entry in scored:
        for block in ("binding", "suitability"):
            for c in entry[block]["components"]:
                if c["state"] == scoring.MEASURED and c["key"] not in axes:
                    axes.append(c["key"])

    front: list[str] = []
    def measured_values(entry) -> dict:
        """Every component this binder actually measured, by key."""
        values = {}
        for block in ("binding", "suitability"):
            for c in entry[block]["components"]:
                if c["state"] == scoring.MEASURED:
                    values[c["key"]] = c["value"]
        return values

    incomparable = []
    if axes:
        # A binder is comparable only if it measured every axis the front is
        # computed over. Substituting 0.0 for an axis it did not measure would
        # impute missing evidence as the worst possible value and make it
        # dominated by anything that did measure it -- which is the imputation
        # this platform refuses everywhere else, in the one place that is
        # supposed to be free of it.
        comparable = []
        for entry in scored:
            values = measured_values(entry)
            if all(a in values for a in axes):
                comparable.append((entry, tuple(values[a] for a in axes)))
            else:
                missing = [a for a in axes if a not in values]
                incomparable.append({"binder_id": entry["binder_id"],
                                     "unmeasured_axes": missing})

        points = [pt for _e, pt in comparable]
        front = [comparable[i][0]["binder_id"]
                 for i in stage11.pareto_front(points)] if points else []
        for entry in scored:
            entry["on_front"] = entry["binder_id"] in front
            entry["front_comparable"] = not any(
                x["binder_id"] == entry["binder_id"] for x in incomparable)
        # A front whose members are all the same point has not discriminated
        # between them: it is the shape of an answer without being one, and
        # reporting membership alone would read as a result.
        on_front_points = {tuple(measured_values(e)[a] for a in axes)
                           for e in scored if e["on_front"]}
        discriminates = len(on_front_points) > 1
    else:
        for entry in scored:
            entry["on_front"] = False
            entry["front_comparable"] = False
        discriminates = False

    return {
        "status": RANKED if any(s["binding_rank"] or s["car_suitability_rank"]
                                for s in scored) else NOT_SCORED,
        "weight_versions": {"binding": BINDING_VERSION,
                            "suitability": SUITABILITY_VERSION},
        "binders": scored,
        "excluded_binders": [
            {"binder_id": r.get("binder_id"), "identifier": r.get("identifier"),
             "target_match": r.get("target_match"),
             "recorded_antigens": r.get("recorded_antigens"),
             "decision": EXCLUDED,
             "reason": "the recorded antigen names a different protein, so "
                       "this binder is excluded before scoring and no weight "
                       "can rescue it"}
            for r in excluded
        ],
        "front_axes": axes,
        "front": front,
        "front_discriminates": discriminates,
        "front_incomparable": incomparable,
        "front_note": (
            f"{len(front)} binder(s) are on the front over {len(axes)} measured "
            f"axis/axes, and they are all the same point. Every one is tied at "
            f"the same value on every measured axis, so none dominates any "
            f"other and all are non-dominated by default. This front is "
            f"populated and carries no discrimination: it says nothing "
            f"distinguishes these binders on what was measured, which is the "
            f"same statement an empty front would make."
            if front and not discriminates else
            f"{len(front)} binder(s) are non-dominated over {len(axes)} "
            f"measured axis/axes, and they differ from one another."
            if front else
            "No component carries a value, so there is nothing for the front "
            "to compare and it is empty."),
        "reasons": [
            "Two rankings, two declared weight sets, and a front computed from "
            "component values rather than weights. Where the rankings disagree "
            "the front says which binders are not beaten on every measured "
            "axis at once, and no choice of weights can move it.",
            "Binding rank asks which binder appears to bind most strongly, so "
            "affinity carries the most weight in it. CAR-suitability rank "
            "deliberately does not put affinity on top, because the document's "
            "own rule is that highest affinity must not equal best CAR binder.",
            "A binder whose recorded antigen names a different protein is "
            "excluded before scoring rather than scored badly. Level A first.",
            "A binder that did not measure every axis the front is computed "
            "over is reported as incomparable rather than placed on the front "
            "with zeros standing in for what it did not measure. Substituting "
            "a value for missing evidence would make it dominated by anything "
            "that measured that axis, which is imputation deciding the "
            "comparison.",
        ],
    }
