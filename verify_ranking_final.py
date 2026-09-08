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
    binder_check, scoring, stage3, stage4, stage5, stage6, stage9, stage10, stage11,
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
    surface_by_gene = {r.gene: r for r in surface if r.gene}
    rows, attrition, status = stage11.rank(
        decisions, binders, constructs, gated, liabilities, composites,
        ceiling, indication_key=PDAC.key, stage3_rows=stage3_by_gene,
        budget_bp=stage6.BUDGET_BP, surface_records=surface_by_gene)

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
            budget_bp=stage6.BUDGET_BP, surface_records=surface_by_gene)
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

    # ------------------------------------------------------------------
    # W12-W18: counting target-matched binders rather than database hits.
    # The criteria were fixed in specs/binder-count-correction.md before any
    # of this was written, and four of them assert that something did NOT move.
    # ------------------------------------------------------------------

    survivors = [r for r in rows if r.survived]

    # W12: the three counts are arithmetically closed, so no binder is lost or
    # double-counted between what the search returned and what the ranking
    # used. The pin matters: with no wrong-antigen binder anywhere in the pool
    # the closure holds trivially and proves nothing.
    w12_bad = [f"{r.gene}: {r.retrieved_binder_count} retrieved != "
               f"{r.binder_count} counted + {r.wrong_antigen_binder_count} wrong"
               for r in rows
               if r.retrieved_binder_count
               != r.binder_count + r.wrong_antigen_binder_count]
    flagged = [r for r in rows if r.wrong_antigen_binder_count >= 1]
    if not flagged:
        w12_bad.append("no candidate carries a wrong-antigen binder, so the "
                       "closure holds trivially and this criterion is vacuous")
    criterion("W12", bool(w12_bad),
              f"retrieved = counted + wrong-antigen on all {len(rows)} "
              f"candidate(s), with {len(flagged)} carrying a wrong-antigen "
              f"binder: "
              + ", ".join(f"{r.gene} {r.retrieved_binder_count}="
                          f"{r.binder_count}+{r.wrong_antigen_binder_count}"
                          for r in flagged[:3])
              if not w12_bad else "; ".join(w12_bad))

    # W13: the basis string reaches the configuration hash. Without it the
    # payload carries nothing describing how an objective is computed, and a
    # cached run under the superseded basis compares equal to this one.
    basis_before = stage11.BINDER_COUNT_BASIS
    hash_basis_now = stage11.configuration_hash(stage9_hash, genes)
    try:
        stage11.BINDER_COUNT_BASIS = "every row the search returned"
        hash_basis_alt = stage11.configuration_hash(stage9_hash, genes)
    finally:
        stage11.BINDER_COUNT_BASIS = basis_before
    criterion("W13", hash_basis_now == hash_basis_alt,
              f"changing the binder-count basis moves the Stage 11 hash "
              f"{hash_basis_now} -> {hash_basis_alt}"
              if hash_basis_now != hash_basis_alt else
              f"the hash stayed {hash_basis_now} when the basis changed, so a "
              "run counting hits compares equal to a run counting matches")

    # W14: both decisions are carried and the correction actually reached the
    # ranking. If nothing moved, the change had no effect and saying so is the
    # result, not a pass.
    moved_decision = [r for r in rows if r.decision_changed]
    criterion("W14", not moved_decision,
              "; ".join(f"{r.gene} {r.decision_under_retrieved_count} -> "
                        f"{r.decision}" for r in moved_decision)
              + " -- both decisions carried on the record"
              if moved_decision else
              "no candidate decision differs between the retrieved and the "
              "counted objective, so the correction did not reach the ranking")

    # W15/W16: rank the same inputs again with the counting rule swapped back
    # to counting every hit, and require the gates and the Level B scores to be
    # untouched. The swap is verified to have taken effect before either
    # criterion is allowed to clear on an identity.
    original_rule = binder_check.counts_towards_objective
    try:
        binder_check.counts_towards_objective = lambda row: True
        hits_rows, hits_attrition, _hs = stage11.rank(
            decisions, binders, constructs, gated, liabilities, composites,
            ceiling, indication_key=PDAC.key, stage3_rows=stage3_by_gene,
            budget_bp=stage6.BUDGET_BP, surface_records=surface_by_gene)
    finally:
        binder_check.counts_towards_objective = original_rule

    hits_by_gene = {r.gene: r for r in hits_rows}
    front_counted = sorted(r.gene for r in rows if r.on_front)
    front_hits = sorted(r.gene for r in hits_rows if r.on_front)
    swap_took_effect = front_counted != front_hits

    w15_bad = []
    if not swap_took_effect:
        w15_bad.append(
            f"the front is {front_counted} under both counting rules, so the "
            "swap did not take effect and these criteria prove nothing")
    if hits_attrition != attrition:
        w15_bad.append(f"attrition moved {attrition} -> {hits_attrition}")
    gate_moved = [f"{r.gene} {r.gate_status} -> {hits_by_gene[r.gene].gate_status}"
                  for r in rows
                  if r.gene in hits_by_gene
                  and r.gate_status != hits_by_gene[r.gene].gate_status]
    if gate_moved:
        w15_bad.append("gate status moved: " + "; ".join(gate_moved))
    criterion("W15", bool(w15_bad),
              f"the gate did not move: attrition identical and all "
              f"{len(rows)} gate statuses identical, while the front moved "
              f"{front_hits} -> {front_counted}"
              if not w15_bad else "; ".join(w15_bad))

    w16_bad = []
    if not swap_took_effect:
        w16_bad.append("the counting swap did not take effect")
    for r in survivors:
        other = hits_by_gene.get(r.gene)
        if other is None or r.scorecard is None or other.scorecard is None:
            continue
        for label, a, b in (
                ("overall", r.overall, other.overall),
                ("fraction", r.scorecard.fraction, other.scorecard.fraction),
                ("applicable", r.scorecard.applicable, other.scorecard.applicable),
                ("measured", r.scorecard.measured_weight,
                 other.scorecard.measured_weight)):
            if (a is None) != (b is None):
                w16_bad.append(f"{r.gene} {label} {a} vs {b}")
            elif a is not None and abs(a - b) > 1e-12:
                w16_bad.append(f"{r.gene} {label} moved {a} -> {b}")
    criterion("W16", bool(w16_bad),
              f"Level B did not move: overall, scored fraction, applicable and "
              f"measured weight identical across all {len(survivors)} "
              f"survivor(s) under both counting rules"
              if not w16_bad else "; ".join(w16_bad))

    # W17: no estimated affinity is introduced anywhere. The correction is
    # about which binders count, not about inventing a measurement for them.
    numeric_affinity = []
    for gene, record in binders.items():
        for c in list(record.structure) + list(record.sequence):
            value = getattr(c, "affinity", None)
            if isinstance(value, (int, float)):
                numeric_affinity.append(f"{gene} {c.identifier} affinity={value}")
    criterion("W17", bool(numeric_affinity),
              "no binder carries a numeric affinity; every one reads the "
              "not-connected token, and fold error and correlation are not "
              "computed from it"
              if not numeric_affinity else
              "estimated affinity present: " + "; ".join(numeric_affinity[:3]))

    # W18: a design that counts a binder says which path counted it, so a count
    # resting on absent annotation is not mistaken for one resting on a match.
    w18_bad = []
    for r in rows:
        total = sum(r.binder_count_paths.values())
        if total != r.binder_count:
            w18_bad.append(f"{r.gene}: paths sum to {total}, "
                           f"binder_count is {r.binder_count}")
        if r.binder_count >= 1 and not r.binder_count_paths:
            w18_bad.append(f"{r.gene} counts {r.binder_count} binder(s) with "
                           "no path breakdown emitted")
    counting = [r for r in rows if r.binder_count >= 1]
    criterion("W18", bool(w18_bad),
              f"the path breakdown sums to the counted total on all "
              f"{len(rows)} candidate(s); {len(counting)} count at least one "
              f"binder and each names the path"
              if not w18_bad else "; ".join(w18_bad))

    # ------------------------------------------------------------------
    # M1-M8: tolerant antigen matching. Criteria fixed in
    # specs/tolerant-antigen-matching.md before the implementation existed.
    # Most of them assert what the matcher must REFUSE, because a tolerant
    # matcher that over-matches is worse than the error it fixes.
    # ------------------------------------------------------------------

    def exact_only(element, lowered, _usable):
        """The superseded rule: exact equality on the whole element."""
        return (binder_check.MATCHED
                if element.strip().lower() in lowered else None)

    def verdict(gene, antigen, tolerant=True):
        """One verdict for one target, under either matching rule."""
        record = surface_by_gene.get(gene)
        if record is None:
            return None
        if tolerant:
            return binder_check.evaluate(antigen, binder_check.name_set(record))
        original = binder_check.match_element
        try:
            binder_check.match_element = exact_only
            return binder_check.evaluate(antigen, binder_check.name_set(record))
        finally:
            binder_check.match_element = original

    ACCEPT = [
        ("CLDN18", "Isoform A2 of Claudin-18"),
        ("NT5E", "5'-nucleotidase, ecto (CD73), isoform CRA_a"),
        ("CDH1", "Ubiquitin-like protein SMT3,Cadherin-1"),
        ("MUC16", "Maltose/maltodextrin-binding periplasmic protein,Mucin-16"),
        ("MUC1", "MUC1 Peptide Fragment"),
        ("TMPRSS2", "Inactive TMPRSS2 construct"),
    ]
    REFUSE = [
        ("MET", "Hepatocyte growth factor beta chain",
         "the ligand, not the receptor"),
        ("ERBB3", "Receptor tyrosine-protein kinase erbB-2", "a sibling receptor"),
        ("ITGB6", "Integrin alpha-V heavy chain", "the partner chain"),
        ("ITGA6", "Integrin beta-1", "the partner chain"),
        ("CLDN18", "Claudin-1", "a different claudin in the same pool"),
        ("MUC1", "Mucin-16", "a different mucin in the same pool"),
        ("GPR35", "Guanine nucleotide-binding protein subunit alpha-13",
         "the stabilising G-protein"),
        # Found by the sibling sweep below, and the reason word-order
        # tolerance was given up: each has the same word set as one of the
        # target own names and belongs to a different protein.
        ("SLC38A1", "N-system amino acid transporter 1",
         "a paralogue whose name reorders the same words"),
        ("LRRN2", "Leucine-rich repeat neuronal 2 protein",
         "a paralogue whose name reorders the same words"),
        ("SIGLEC12", "Sialic acid-binding Ig-like lectin 1",
         "a paralogue whose name drops a repeated word"),
        ("CELSR2", "Multiple EGF-like domains protein 2",
         "a paralogue with offset family numbering"),
        ("AOC2", "Amine oxidase [copper-containing] 3",
         "a family member reached through a bracket-truncated stem"),
    ]

    # M1: every accept-case moves, and was genuinely failing before. Asserting
    # only that they pass now would clear against a matcher that always passes.
    m1_bad = []
    for gene, antigen in ACCEPT:
        before = verdict(gene, antigen, tolerant=False)
        after = verdict(gene, antigen)
        if before is None or after is None:
            m1_bad.append(f"{gene}: no reference record")
        elif before["target_match"] != binder_check.FAIL:
            m1_bad.append(f"{gene}: was {before['target_match']}, not FAIL, so "
                          "this case proves nothing")
        elif after["target_match"] != binder_check.PASS:
            m1_bad.append(f"{gene}: still {after['target_match']}")
    # Word-order variants were an accepted shape in the specification and are
    # no longer one. Set-based matching bought them at the cost of conflating
    # paralogues across whole families, so the shape was given up rather than
    # the conflation tolerated. The case is asserted as still failing, so the
    # limitation stays visible and cannot quietly reverse.
    word_order = verdict("ERBB2", "Receptor protein-tyrosine kinase erbB-2")
    if word_order and word_order["target_match"] != binder_check.FAIL:
        m1_bad.append("the word-order shape matches again; it was given up "
                      "deliberately and its return needs review, not silence")

    criterion("M1", bool(m1_bad),
              f"all {len(ACCEPT)} label shapes move FAIL -> PASS: isoform, "
              f"fusion partner, fragment and construct labels; the word-order "
              f"shape stays FAIL by decision"
              if not m1_bad else "; ".join(m1_bad))

    # M2: the refusals. This is the criterion the change exists to satisfy.
    m2_bad = [f"{gene} now matches {antigen!r} ({why})"
              for gene, antigen, why in REFUSE
              if (verdict(gene, antigen) or {}).get("target_match")
              != binder_check.FAIL]
    criterion("M2", bool(m2_bad),
              f"all {len(REFUSE)} near-miss cases stay FAIL, the ligand of a "
              f"receptor and two sibling receptors among them"
              if not m2_bad else "; ".join(m2_bad))

    # M3: degenerate reference names are excluded, and the guard is shown to
    # fire. Without the pin this clears on a pool that happens to have none.
    excluded_total, excluded_example = 0, None
    for gene, record in surface_by_gene.items():
        names = binder_check.name_set(record)
        usable = {n for n, _t in binder_check.usable_names(names)}
        dropped = {n.strip() for n in names} - usable
        if dropped:
            excluded_total += len(dropped)
            if excluded_example is None:
                excluded_example = (gene, sorted(dropped)[0])
    leaked = [n for _g, r in surface_by_gene.items()
              for n, t in binder_check.usable_names(binder_check.name_set(r))
              if len(n) < binder_check.MINIMUM_NAME_LENGTH
              or not any(c.isalpha() for c in n)]
    m3_bad = []
    if leaked:
        m3_bad.append(f"{len(leaked)} degenerate name(s) reached the matcher, "
                      f"e.g. {leaked[0]!r}")
    if not excluded_total:
        m3_bad.append("no reference name was excluded anywhere in the pool, so "
                      "the guard is untested")
    criterion("M3", bool(m3_bad),
              f"{excluded_total} degenerate reference name(s) excluded from "
              f"tolerant matching across the pool, e.g. {excluded_example[1]!r} "
              f"on {excluded_example[0]}; none reached the matcher"
              if not m3_bad else "; ".join(m3_bad))

    # M4: sibling names differing only by a trailing number. Found in the pool
    # rather than written down, so the criterion cannot go stale against it.
    def digitless(tokens):
        return frozenset(t for t in tokens if not t.isdigit())

    by_shape = {}
    for gene, record in surface_by_gene.items():
        for name in binder_check.name_set(record):
            tokens = binder_check.name_tokens(name)
            digits = frozenset(t for t in tokens if t.isdigit())
            if not digits or len(tokens) < 2:
                continue
            by_shape.setdefault(digitless(tokens), []).append(
                (gene, name.strip(), digits))

    # This probe was wrong twice before it was right, and both faults are
    # recorded rather than quietly removed.
    #
    # It first paired enzyme classification codes -- "EC 7.6.2.1" against
    # "EC 6.2.1.7" -- as sibling names. They are not names differing by a
    # trailing number, they are structured codes, and the matcher was fixed for
    # that: a reference name carrying more than one digit-only word is no
    # longer usable for tolerant matching.
    #
    # It then still paired them, because two enzymes may legitimately share a
    # classification code. "EC 7.6.2.1" is ABCB1's OWN name, so the probe was
    # testing a target against itself and calling the correct answer a failure.
    #
    # Narrowing a criterion after it trips is the move that deserves the most
    # suspicion, so the sweep is pinned twice below: a shared stem must be a
    # real word, a candidate name must not belong to the target, and the sweep
    # must still find the claudins and the mucins it exists for.
    def real_stem(shape):
        return any(len(t) >= 4 and t.isalpha() for t in shape)

    siblings = []
    for shape, members in by_shape.items():
        if len({g for g, _n, _d in members}) < 2 or not real_stem(shape):
            continue
        for gene_a, _name_a, dig_a in members:
            own = {n.strip().lower() for n in
                   binder_check.name_set(surface_by_gene[gene_a])}
            for gene_b, name_b, dig_b in members:
                if gene_a == gene_b or dig_a == dig_b:
                    continue
                if name_b.lower() in own:
                    continue
                siblings.append((gene_a, gene_b, name_b))

    m4_bad = []
    for gene_a, _gene_b, name_b in siblings:
        got = verdict(gene_a, name_b)
        if got and got["target_match"] == binder_check.PASS:
            m4_bad.append(f"{gene_a} matched {name_b!r}, a sibling name")

    swept = {g for g, _b, _n in siblings}
    if len(siblings) < 20:
        m4_bad.append(f"only {len(siblings)} sibling pair(s) found, too few "
                      "for the sweep to mean anything")
    for wanted in ("CLDN18", "MUC1"):
        if wanted not in swept:
            m4_bad.append(f"{wanted} is not in the sweep, so the narrowing has "
                          "removed a case this criterion exists for")
    criterion("M4", bool(m4_bad),
              f"{len(siblings)} sibling-name pair(s) across {len(swept)} "
              f"target(s), none cross-matching; the claudins and the mucins "
              f"are among them"
              if not m4_bad else "; ".join(sorted(set(m4_bad))[:3]))

    # M5: a relational word in the surplus refuses the match. Constructed on
    # every pool target from its own reference name, not just on MET.
    # The first version of this probe appended a relational word to one of the
    # target's own names and required a refusal. It tripped on targets where
    # the constructed string was a legitimate name of that same target: the
    # gonadotropin-releasing hormone receptor is recorded as both "GnRH-R" and
    # "GnRH receptor", so "GnRH-R receptor" is the target and matching it is
    # correct. It also tripped where the chosen name already contained a comma,
    # which made the probe a fusion label rather than a modified name. Both
    # were faults in the probe, not in the matcher, and the probe is narrowed
    # rather than the result explained. The live case it exists for -- a
    # receptor and its own ligand -- is asserted by name so the narrowing
    # cannot quietly remove what the criterion is for.
    m5_bad = []
    m5_tested = []
    for gene, record in surface_by_gene.items():
        names = binder_check.name_set(record)
        usable = binder_check.usable_names(names)
        if not usable:
            continue
        for word in ("receptor", "ligand", "binding"):
            if any(word in n.lower() for n in names):
                continue
            name, tokens = next(
                ((n, t) for n, t in usable if "," not in n), (None, None))
            if name is None:
                continue
            m5_tested.append(gene)
            got = verdict(gene, f"{name} {word}")
            if got and got["target_match"] == binder_check.PASS:
                m5_bad.append(f"{gene} matched {name + ' ' + word!r}")
            break

    # The named instance, which is why this criterion exists at all.
    met = verdict("MET", "Hepatocyte growth factor beta chain")
    if met is None:
        m5_bad.append("MET is absent, so the receptor-and-its-ligand case "
                      "could not be tested")
    elif met["target_match"] != binder_check.FAIL:
        m5_bad.append("MET matched its own ligand, hepatocyte growth factor "
                      "beta chain")
    if len(m5_tested) < 50:
        m5_bad.append(f"only {len(m5_tested)} target(s) reached the sweep, too "
                      "few for it to mean anything")
    criterion("M5", bool(m5_bad),
              f"a relational word added to a target's own reference name "
              f"refuses the match across all {len(m5_tested)} target(s) swept, "
              f"and MET does not match hepatocyte growth factor beta chain, "
              f"which is its ligand rather than the receptor"
              if not m5_bad else "; ".join(m5_bad[:3]))

    # M6/M7: over the whole pool, under both rules.
    def survey(tolerant):
        """Verdicts for every structural binder row in the pool."""
        original = binder_check.match_element
        if not tolerant:
            binder_check.match_element = exact_only
        try:
            out = {}
            for row in decisions:
                gene = row["gene"]
                record = binders.get(gene)
                if record is None:
                    continue
                out[gene] = binder_check.check(
                    [{"route": "structure", "identifier": c.identifier,
                      "antigen_name": c.antigen_name} for c in record.structure],
                    surface_by_gene.get(gene))
            return out
        finally:
            binder_check.match_element = original

    before_survey = survey(False)
    after_survey = survey(True)
    flat_before = [r for rows_ in before_survey.values() for r in rows_]
    flat_after = [r for rows_ in after_survey.values() for r in rows_]
    fail_before = sum(1 for r in flat_before if r["target_match"] == binder_check.FAIL)
    fail_after = sum(1 for r in flat_after if r["target_match"] == binder_check.FAIL)
    genes_before = sum(1 for rows_ in before_survey.values()
                       if any(r["target_match"] == binder_check.FAIL for r in rows_))
    genes_after = sum(1 for rows_ in after_survey.values()
                      if any(r["target_match"] == binder_check.FAIL for r in rows_))

    # M6: monotonicity. Tolerance may only add matches.
    lost = []
    for gene, rows_ in before_survey.items():
        for i, row in enumerate(rows_):
            if row["target_match"] != binder_check.PASS:
                continue
            now = after_survey[gene][i]
            if now["target_match"] != binder_check.PASS:
                lost.append(f"{gene} {row['identifier']} "
                            f"PASS -> {now['target_match']}")
    criterion("M6", bool(lost),
              f"every one of the {sum(1 for r in flat_before if r['target_match'] == binder_check.PASS)} "
              f"row(s) matching under exact equality still matches; tolerance "
              f"only adds"
              if not lost else "; ".join(lost[:3]))

    # M7: the flagged count moved, reported against the prediction.
    PREDICTED_ROWS, PREDICTED_GENES = 57, 21
    m7_bad = []
    if fail_after == fail_before:
        m7_bad.append(f"the flagged count stayed at {fail_before}, so the rule "
                      "did not take effect")
    criterion("M7", bool(m7_bad),
              f"flagged rows {fail_before} -> {fail_after} "
              f"(predicted {PREDICTED_ROWS}), affected targets "
              f"{genes_before} -> {genes_after} (predicted {PREDICTED_GENES})"
              if not m7_bad else "; ".join(m7_bad))

    # M8: nothing moves a decision without carrying both.
    m8_bad = [f"{r.gene} carries no superseded decision"
              for r in rows if r.decision_changed
              and not r.decision_under_retrieved_count]
    moved_now = [r.gene for r in rows if r.decision_changed]
    criterion("M8", bool(m8_bad),
              f"every decision that moved carries both: {moved_now}"
              if not m8_bad else "; ".join(m8_bad))

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
