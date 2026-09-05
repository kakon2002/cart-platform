"""Runs the final ranking and tests it against its criteria."""

from __future__ import annotations

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
    stage3, stage4, stage5, stage6, stage9, stage10, stage11, validation,
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

    rows, attrition, status = stage11.rank(
        decisions, binders, constructs, gated, liabilities, composites,
        ceiling, indication_key=PDAC.key)

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

    summed = any(hasattr(r, "total_score") or hasattr(r, "weighted") for r in rows)
    criterion("N4", summed,
              "no weighted or summed score across objectives is emitted")

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
