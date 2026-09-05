"""Runs the final ranking and tests it against its criteria."""

from __future__ import annotations

import dataclasses
import sys

from car_pipeline.configs.pdac import PDAC, PDAC_PROJECT
from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
from car_pipeline.data.trials import TrialSource
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import (
    scoring, stage3, stage4, stage5, stage6, stage9, stage10, stage11,
    validation,
)
from car_pipeline.stages.stage1 import build_spec


def main() -> int:
    """Run the final-ranking criteria."""
    print("loading every upstream stage", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    records = stage5.load_or_retrieve(
        decisions, AntibodySource(), manifest["stage4_hash"])
    binders = {r.gene: r for r in records}
    constructs = {c.gene: c for c in stage6.build(decisions, binders)}

    surface, _ = load_surface()
    atlas = HPASource().load()
    by_acc, by_sym = atlas_index(atlas)
    gtex_profiles, gtex_tissues, _ = GTExSource().match_surface(surface, by_acc)
    cohort = TCGASource().load()
    cohort_join = tcga_match(cohort, surface, by_acc)
    cells = SingleCellSource().load()
    cell_index = cell_match(cells, surface, by_acc)
    dependency, _ = DepMapSource().load()
    dep_index = gene_index(dependency)
    coverage_rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, cohort_join)
    spec = build_spec(PDAC_PROJECT)
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    overrides = {o: ov.tier
                 for o, ov in spec.inputs.tissue_criticality_overrides.items()}
    calibration = stage3.calibrate_atlas_levels(
        surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
        stage3.RiskModel(overrides=overrides))
    ranked3, _m, _s = stage3.rank(
        coverage_rows, {r.accession: r for r in surface}, by_acc, by_sym,
        cells, cell_index, gtex_profiles, gtex_tissues, cohort, cohort_join,
        dependency, dep_index, overrides, ceiling, calibration)
    composites = {r.gene: r.composite for r in ranked3 if r.gene}
    risks = {r.gene: (r.risk, r.risk_organ) for r in ranked3 if r.gene}

    genes = [d["gene"] for d in decisions]
    trials = TrialSource(antigens=genes).load()
    gated = {g.gene: g for g in stage9.gate(
        decisions, binders, risks, trials, ceiling,
        constructs=list(constructs.values()))}

    stage5_hash = stage5.configuration_hash(
        manifest["stage4_hash"], [r.gene for r in records])
    stage6_hash = stage6.configuration_hash(
        stage5_hash, [c.gene for c in constructs.values()])
    stage9_hash = stage9.configuration_hash(
        stage6_hash, [d["gene"] for d in decisions], ceiling)

    dev_rows, _status = stage10.assess(binders)
    liabilities: dict[str, list] = {}
    for row in dev_rows:
        liabilities.setdefault(row.gene, []).append(row)

    stage3_by_gene = {r.gene: r for r in ranked3 if r.gene}
    rows, attrition, status = stage11.rank(
        decisions, binders, constructs, gated, liabilities, composites,
        ceiling, indication_key=PDAC.key, stage3_rows=stage3_by_gene,
        budget_bp=stage6.BUDGET_BP)

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    checked: list[str] = []
    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    a, b, c = (5.0, 5.0, 5.0, 5.0), (1.0, 1.0, 1.0, 1.0), (9.0, 0.0, 0.0, 0.0)
    front = stage11.pareto_front([a, b, c])
    criterion("N1", 1 in front,
              "a dominated point is excluded from the front"
              if 1 not in front else "a dominated point appears on the front")
    criterion("N2", not (0 in front and 2 in front),
              "both non-dominated points are on the front"
              if 0 in front and 2 in front else f"front was {front}, expected 0 and 2")

    total = sum(attrition.values()) + sum(1 for r in rows if r.survived)
    criterion("N3", total != len(decisions),
              f"attrition accounts for {total} of {len(decisions)}")

    # N4 re-specified. The original -- no weighted sum is emitted -- was written
    # when ranking and gating were one step, and what it protected was that a
    # strong efficacy score must not buy off a weak safety one. The two-level
    # structure now carries that guarantee: a candidate over the applied ceiling
    # never reaches Level B. What N4 forbids is therefore restated as the two
    # things that would actually break it.
    rescued = [r.gene for r in rows if not r.survived and r.overall is not None]
    imputed = [f"{r.gene}.{k}" for r in rows if r.scorecard
               for k, c in r.scorecard.components.items()
               if c.state != scoring.MEASURED and c.value is not None]
    criterion("N4", bool(rescued or imputed),
              f"no gate failure carries a score ({sum(1 for r in rows if not r.survived)} "
              f"failed, {sum(1 for r in rows if r.overall is not None)} scored) and "
              f"no component is imputed"
              if not (rescued or imputed) else
              f"rescued {rescued[:3]}; imputed {imputed[:3]}")

    criterion("N5",
              status == stage11.RANKED and not any(r.survived for r in rows),
              f"status {status} matches the survivor count "
              f"{sum(1 for r in rows if r.survived)}")

    criterion("N6", len(rows) != manifest["pool_size"],
              f"{len(rows)} rows against the {manifest['pool_size']} the Stage 4 "
              "manifest records")

    # W8 and W9 land with the decision column rather than at step 6, because
    # they are that column's criteria and it ships now.
    classes = [validation.CONSERVATIVE, validation.INNOVATIVE]
    undeclared = sorted({r.decision for r in rows} - set(stage11.DECISIONS))
    overlap = [f"{d} <-> {k}" for d in stage11.DECISIONS for k in classes
               if d in k or k in d]
    criterion("W8", bool(undeclared or overlap),
              f"{len(set(r.decision for r in rows))} decision value(s) all "
              f"declared, and none overlaps either design class by substring"
              if not (undeclared or overlap) else
              f"undeclared {undeclared}; overlapping {overlap}")

    # Recomputed from gate status and the front alone. Design class is not read,
    # which is the independence W9 exists to assert.
    #
    # The mapping is written out here rather than read from stage11.GATE_DECISION
    # or recomputed by stage11.decision_for. Calling either would compare the
    # subject against itself and clear whatever it produced -- the shape this
    # repository has now recorded thirteen times. These literals are the pin: a
    # change to either the gate tokens or the gate-to-decision mapping trips
    # this, which is the point.
    RECOVERABLE = {"NO_BINDER_RETRIEVED", "NO_CONSTRUCT_ASSEMBLED"}
    TERMINAL = {"BLOCKED_ON_NORMAL_TISSUE_RISK", "NO_DESIGN_RECOMMENDED",
                "OVER_PAYLOAD_BUDGET"}

    def expected(r) -> str:
        """The decision this row must carry, derived without stage11's mapping."""
        if not r.survived:
            if r.gate_status in TERMINAL:
                return "EXCLUDED"
            if r.gate_status in RECOVERABLE:
                return "REQUIRES_EVIDENCE"
            return f"UNMAPPED_GATE:{r.gate_status}"
        return "ADVANCE" if r.on_front else "BACKUP"

    mismatched = [f"{r.gene} carries {r.decision}, expected {expected(r)}"
                  for r in rows if r.decision != expected(r)]
    survivors_ = [r for r in rows if r.survived]
    splits = len({r.on_front for r in survivors_}) > 1
    undistinguished = splits and len({r.decision for r in survivors_}) == 1
    criterion("W9", bool(mismatched or undistinguished),
              f"every decision recomputes from gate status and the front alone; "
              f"{sum(1 for r in survivors_ if r.on_front)} of {len(survivors_)} "
              f"survivors are on the front and the decisions distinguish them"
              if not (mismatched or undistinguished) else
              f"{len(mismatched)} decision(s) do not recompute: "
              + "; ".join(mismatched[:3])
              if mismatched else
              "the survivors split on the front but all carry one decision")

    # ---------------- Level B, the scoring frame ----------------
    cards = [r.scorecard for r in rows if r.scorecard]

    w_sum = sum(scoring.WEIGHTS.values())
    declared = set(scoring.WEIGHTS) | {"evidence_confidence", "prediction_uncertainty"}
    w1_bad = []
    if abs(w_sum - 1.0) > 1e-12:
        w1_bad.append(f"weights sum to {w_sum!r}, not 1.0")
    if len(declared) != 11:
        w1_bad.append(f"{len(declared)} components declared, expected 11")
    if not scoring.WEIGHT_VERSION:
        w1_bad.append("the weight set carries no version")
    criterion("W1", bool(w1_bad),
              f"nine weights sum to {w_sum:.12g}, all eleven components carry a "
              f"declared coefficient, version {scoring.WEIGHT_VERSION}"
              if not w1_bad else "; ".join(w1_bad))

    w2_bad = []
    for card in cards:
        for key, c in card.components.items():
            if c.state not in scoring.STATES:
                w2_bad.append(f"{card.gene}.{key} state {c.state!r}")
            if c.state != scoring.MEASURED:
                if c.value is not None:
                    w2_bad.append(f"{card.gene}.{key} carries a value while {c.state}")
                if not c.source:
                    w2_bad.append(f"{card.gene}.{key} is {c.state} with no reason")
        # The denominator is the measured weight; nothing else may enter it.
        recomputed = sum(c.weight for c in card.components.values()
                         if c.state == scoring.MEASURED)
        if abs(recomputed - card.measured_weight) > 1e-12:
            w2_bad.append(f"{card.gene} denominator {card.measured_weight} "
                          f"against {recomputed} over MEASURED alone")
    criterion("W2", bool(w2_bad),
              f"{len(cards)} scorecard(s): every component is exactly one of the "
              f"three states, every non-measured one names its reason and carries "
              f"no value, and the denominator is the measured weight alone"
              if not w2_bad else "; ".join(w2_bad[:3]))

    w3_bad = []
    for card in cards:
        if card.overall is None:
            continue
        num = sum(c.weight * c.value for c in card.components.values()
                  if c.state == scoring.MEASURED)
        den = sum(c.weight for c in card.components.values()
                  if c.state == scoring.MEASURED)
        again = (num / den) * card.confidence_adjustment
        if abs(again - card.overall) > 1e-12:
            w3_bad.append(f"{card.gene} records {card.overall} against {again}")
    criterion("W3", bool(w3_bad),
              f"{sum(1 for c in cards if c.overall is not None)} score(s) "
              f"recompute from the components, weights and adjustment recorded "
              f"on the candidate, to within 1e-12"
              if not w3_bad else "; ".join(w3_bad[:3]))

    w4_bad = [r.gene for r in rows
              if not r.survived and (r.overall is not None or r.scorecard is not None)]
    criterion("W4", bool(w4_bad),
              f"{sum(1 for r in rows if not r.survived)} gate failure(s) carry "
              f"neither a score nor a scorecard; scoring is reached by survivors "
              f"only, which is what makes the weighted sum safe"
              if not w4_bad else f"{len(w4_bad)} scored despite failing: {w4_bad[:3]}")

    w5_bad = []
    for card in cards:
        below = card.fraction < scoring.MINIMUM_SCORED_FRACTION
        if below and card.overall is not None:
            w5_bad.append(f"{card.gene} scored at fraction {card.fraction:.4f}")
        if not below and card.overall is None and card.confidence_adjustment is not None:
            w5_bad.append(f"{card.gene} unscored at fraction {card.fraction:.4f}")
    criterion("W5", bool(w5_bad),
              f"every candidate above the {scoring.MINIMUM_SCORED_FRACTION} floor "
              f"carries a number and every candidate below it carries null "
              f"({sum(1 for c in cards if c.overall is not None)} scored of {len(cards)})"
              if not w5_bad else "; ".join(w5_bad[:3]))

    # W6: normal-tissue risk and evidence confidence are never combined. Both
    # directions are perturbed rather than argued from the call signature.
    probe = next(r for r in rows if r.survived)
    base = probe.scorecard
    shifted_conf = dataclasses.replace(
        stage3_by_gene[probe.gene],
        confidence=0.5 if base.evidence_confidence != 0.5 else 0.25)
    card_conf = scoring.score(probe, shifted_conf, gated.get(probe.gene),
                              constructs.get(probe.gene), binders.get(probe.gene),
                              stage6.BUDGET_BP)
    safety_before = base.components["normal_tissue_safety"].value
    safety_after = card_conf.components["normal_tissue_safety"].value

    shifted_risk = dataclasses.replace(
        gated[probe.gene], risk=min(gated[probe.gene].risk * 0.5,
                                    gated[probe.gene].ceiling))
    card_risk = scoring.score(probe, stage3_by_gene[probe.gene], shifted_risk,
                              constructs.get(probe.gene), binders.get(probe.gene),
                              stage6.BUDGET_BP)
    w6_bad = []
    if safety_before != safety_after:
        w6_bad.append(f"halving confidence moved the safety component "
                      f"{safety_before} -> {safety_after}")
    if card_risk.confidence_adjustment != base.confidence_adjustment:
        w6_bad.append(f"halving risk moved the adjustment "
                      f"{base.confidence_adjustment} -> {card_risk.confidence_adjustment}")
    if card_risk.components["normal_tissue_safety"].value == safety_before:
        w6_bad.append("halving risk did not move the safety component, so the "
                      "probe proves nothing")
    criterion("W6", bool(w6_bad),
              f"confidence moved {base.evidence_confidence} -> "
              f"{card_conf.evidence_confidence} and the safety component held at "
              f"{safety_before:.4f}; risk moved the safety component to "
              f"{card_risk.components['normal_tissue_safety'].value:.4f} and the "
              f"adjustment held at {base.confidence_adjustment}"
              if not w6_bad else "; ".join(w6_bad))

    # W7: the front is computed from component values, never from the weights,
    # so no choice of weights can move its membership.
    ALT_WEIGHTS = {
        "tumour_coverage": 0.05, "malignant_specificity": 0.05,
        "normal_tissue_safety": 0.35, "binder_quality": 0.05,
        "manufacturability": 0.35, "developability": 0.05,
        "structural_feasibility": 0.05, "functional_prediction": 0.03,
        "pairing_robustness": 0.02,
    }
    original_weights = dict(scoring.WEIGHTS)
    front_before = {r.gene for r in rows if r.on_front}
    order_before = [r.gene for r in sorted(
        (r for r in rows if r.overall is not None), key=lambda r: -r.overall)]
    try:
        scoring.WEIGHTS.clear()
        scoring.WEIGHTS.update(ALT_WEIGHTS)
        alt_rows, _a, _s2 = stage11.rank(
            decisions, binders, constructs, gated, liabilities, composites,
            ceiling, indication_key=PDAC.key, stage3_rows=stage3_by_gene,
            budget_bp=stage6.BUDGET_BP)
        front_after = {r.gene for r in alt_rows if r.on_front}
        order_after = [r.gene for r in sorted(
            (r for r in alt_rows if r.overall is not None), key=lambda r: -r.overall)]
        scores_after = {r.gene: r.overall for r in alt_rows if r.overall is not None}
    finally:
        scoring.WEIGHTS.clear()
        scoring.WEIGHTS.update(original_weights)

    scores_before = {r.gene: r.overall for r in rows if r.overall is not None}
    # Without this the criterion clears whenever the swap silently fails to
    # take effect, which would make it a check that cannot fail.
    moved = [g for g in scores_before
             if abs(scores_before[g] - scores_after.get(g, scores_before[g])) > 1e-9]
    w7_bad = []
    if front_before != front_after:
        w7_bad.append(f"the front moved {sorted(front_before)} -> "
                      f"{sorted(front_after)} when only the weights changed")
    if not moved:
        w7_bad.append("no score changed under the alternative weights, so the "
                      "swap did not take effect and this criterion proves nothing")
    criterion("W7", bool(w7_bad),
              f"the front is {sorted(front_before)} under both weight sets while "
              f"{len(moved)} of {len(scores_before)} score(s) moved; the score "
              + ("order held at " + str(order_before) if order_before == order_after
                 else f"order moved {order_before} -> {order_after}")
              + " -- order is free to change, membership is not"
              if not w7_bad else "; ".join(w7_bad))

    hash_before = stage11.configuration_hash(stage9_hash, genes)
    try:
        scoring.WEIGHTS["tumour_coverage"] = 0.19
        hash_after = stage11.configuration_hash(stage9_hash, genes)
    finally:
        scoring.WEIGHTS.clear()
        scoring.WEIGHTS.update(original_weights)
    criterion("W11", hash_before == hash_after,
              f"changing one weight moves the Stage 11 hash "
              f"{hash_before} -> {hash_after}"
              if hash_before != hash_after else
              f"the hash stayed {hash_before} when a weight changed, so a run "
              "under different weights compares equal to this one")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHERE THE TWO HUNDRED WENT")
    print("=" * 72)
    running = len(rows)
    for gate in stage11.GATES:
        n = attrition[gate]
        running -= n
        print(f"    {gate:34s} -{n:4d}    {running:4d} remain")
    print(f"    {'reached the end':34s}        {sum(1 for r in rows if r.survived):4d}")

    print()
    if status == stage11.NO_DESIGN_REACHES_THE_END:
        print("  " + "=" * 68)
        print("  NO_DESIGN_REACHES_THE_END")
        print("  " + "=" * 68)
        print("    Reported as a status, not as an empty ranking. An empty table")
        print("    reads as 'nothing ranked highly'; the true statement is that")
        print("    nothing arrived to be ranked, and that difference is the whole")
        print("    result of this pipeline for this indication.")
        print()
        print("    Each drop above is a measurement, not a failure of the stage")
        print("    that made it. The safety ceiling is Stage 1's and it is doing")
        print("    what it exists for; the binder gap is a statement about the")
        print("    literature; the budget overrun is arithmetic against a vector")
        print("    capacity fixed before any of this ran.")
    else:
        front = [r for r in rows if r.on_front]
        print(f"  Pareto front: {len(front)} design(s), no weighted total")
        for r in front:
            print(f"    {r.gene:10s} attractiveness {r.attractiveness}  "
                  f"margin {r.safety_margin}  binders {r.binder_count}  "
                  f"cleanliness {r.cleanliness}")

    print()
    print(f"  configuration hash "
          f"{stage11.configuration_hash(stage9_hash, [r.gene for r in rows])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
