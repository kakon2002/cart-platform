"""Runs the safety gate and tests it against the criteria fixed in the spec.

Criteria before biology, and the two positive pins first. Stage 5's lesson was
that a retrieval route can be dead while every negative check passes.
"""

from __future__ import annotations

import sys

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
from car_pipeline.data.trials import TrialSource
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import stage3, stage4, stage5, stage6, stage9
from car_pipeline.stages.stage1 import build_spec

#: Pinned. Both antigens have many registered trials; zero would mean the route is
#: dead rather than the antigen untried.
PINNED_TRIALS = ["MSLN", "CLDN18"]
#: Pinned. Both are `-xi-` names and must classify chimeric.
PINNED_ORIGIN = {"Amatuximab": "chimeric", "Zolbetuximab": "chimeric"}


def main() -> int:
    print("loading upstream", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    records = stage5.retrieve(decisions, AntibodySource(), progress=False)
    binders = {r.gene: r for r in records}

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
    ranked, _model, _ = stage3.rank(
        coverage_rows, {r.accession: r for r in surface}, by_acc, by_sym,
        cells, cell_index, gtex_profiles, gtex_tissues, cohort, cohort_join,
        dependency, dep_index, overrides, ceiling, calibration)
    risks = {r.gene: (r.risk, r.risk_organ) for r in ranked if r.gene}

    stage5_hash = stage5.configuration_hash(
        manifest["stage4_hash"], [r.gene for r in records])
    stage6_hash = stage6.configuration_hash(
        stage5_hash, [d["gene"] for d in decisions])
    genes = [d["gene"] for d in decisions]
    trials = TrialSource(antigens=genes).load()
    gated = stage9.gate(decisions, binders, risks, trials, ceiling)
    by_gene = {g.gene: g for g in gated}

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    empty = [g for g in PINNED_TRIALS
             if g not in trials or trials[g].total == 0]
    criterion("S1", bool(empty),
              "registry returns trials for both pins ("
              + ", ".join(f"{g}={trials[g].total}" for g in PINNED_TRIALS
                          if g in trials) + ")"
              if not empty else f"no trials returned for {empty}")

    origin_bad = [f"{n}->{stage9.binder_origin(n)}"
                  for n, want in PINNED_ORIGIN.items()
                  if stage9.binder_origin(n) != want]
    criterion("S2", bool(origin_bad),
              "pinned names classify as expected ("
              + ", ".join(f"{n}={stage9.binder_origin(n)}" for n in PINNED_ORIGIN)
              + ")" if not origin_bad else "; ".join(origin_bad))

    missing_risk = [g.gene for g in gated
                    if binders.get(g.gene)
                    and (binders[g.gene].sequence or binders[g.gene].structure)
                    and g.risk is None]
    criterion("S3", bool(missing_risk),
              "every target with a binder carries a Stage 3 risk"
              if not missing_risk else f"no risk for {missing_risk[:5]}")

    contradicts = [g.gene for g in gated
                   if g.risk is not None and g.risk > ceiling
                   and g.verdict != stage9.BLOCKED]
    criterion("S4", bool(contradicts),
              "no target over the ceiling escapes BLOCKED"
              if not contradicts else f"{contradicts[:5]} over the ceiling, not blocked")

    leaked = [g.gene for g in gated
              if g.epitope_immunogenicity != stage9.NOT_CONNECTED]
    criterion("S5", bool(leaked),
              "epitope immunogenicity is NOT_CONNECTED on every row"
              if not leaked else f"{len(leaked)} rows carry a value")

    # Against the manifest's recorded pool size, not against this run's own input.
    # Comparing output to input cannot fail, and the specification pins 200.
    expected_rows = manifest["pool_size"]
    out_genes = {g.gene for g in gated}
    criterion("S6",
              len(gated) != expected_rows or len(out_genes) != expected_rows,
              f"{len(gated)} rows and {len(out_genes)} distinct genes against the "
              f"{expected_rows} the Stage 4 manifest records")

    unflagged = [g.gene for g in gated
                 if g.trials_stopped and g.verdict == stage9.PASSES]
    criterion("S7", bool(unflagged),
              "every target with stopped trials is flagged or blocked"
              if not unflagged else f"{unflagged[:5]} passed with stopped trials")

    print("=" * 72)
    print(f"  {7 - len(tripped)}/7 criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    # ---------------- the biology ---------------------------------------
    print()
    print("=" * 72)
    print("WHAT THE GATE SAYS")
    print("=" * 72)
    counts: dict[str, int] = {}
    for g in gated:
        counts[g.verdict] = counts.get(g.verdict, 0) + 1
    for verdict in (stage9.PASSES, stage9.FLAGGED, stage9.BLOCKED, stage9.NO_GATE):
        n = counts.get(verdict, 0)
        print(f"    {verdict:22s} {n:4d}  ({n / len(gated):.0%})")
    # Which of this gate's own questions actually ran. Reported because a
    # criterion that passes on a code path nothing reached has not been tested,
    # and Stage 5 has already shown what that looks like.
    reached = [g for g in gated if g.verdict in (stage9.PASSES, stage9.FLAGGED)]
    print()
    print("  WHICH QUESTIONS THIS RUN ACTUALLY EXERCISED")
    print("  " + "-" * 68)
    print(f"    reached the immunogenicity and trials questions   {len(reached)} "
          f"of {len(gated)}")
    print(f"    stopped at carried Stage 3 risk                   "
          f"{sum(1 for g in gated if g.verdict == stage9.BLOCKED)}")
    print(f"    stopped for want of a binder                      "
          f"{sum(1 for g in gated if g.verdict == stage9.NO_GATE)}")
    if not reached:
        print()
        print("    So this gate's own two questions decided nothing here. Every")
        print("    target was settled by the risk Stage 3 already measured, which")
        print("    is the correct precedence and leaves the new logic untested on")
        print("    live data. S7 passes vacuously and should be read that way.")
        print("    S2 does not: it pins the origin rule against known names rather")
        print("    than against this run's output, which is why it still has force.")
        print()
        print("    The origin table below is computed for every target carrying a")
        print("    binder, including blocked ones, so it is measured rather than")
        print("    vacuous — it simply did not change any verdict.")

    print()
    print("    PASSES_STATED_CHECKS is not a safety claim. It means three named")
    print("    questions failed to show a problem, and one of them — epitope-level")
    print("    immunogenicity — was not asked at all.")

    print()
    print("  The five known targets")
    print("  " + "-" * 68)
    for gene in ("CEACAM5", "CEACAM6", "CLDN18", "MSLN", "MUC1"):
        g = by_gene.get(gene)
        if g is None:
            continue
        risk = f"{g.risk:.4f}" if g.risk is not None else "n/a"
        print(f"    {gene:9s} {g.verdict:22s} risk {risk:>7s} ({g.risk_organ}) "
              f"binder {g.binder_name or '-'} [{g.binder_origin}]")
        print(f"      trials {g.trials_total:4d}, stopped {g.trials_stopped}")
        for reason in g.reasons:
            print(f"      - {reason}")

    origins: dict[str, int] = {}
    for g in gated:
        if g.binder_name:
            origins[g.binder_origin] = origins.get(g.binder_origin, 0) + 1
    print()
    print("  Binder origin across targets that have one, read from the name stem")
    for origin, n in sorted(origins.items(), key=lambda kv: -kv[1]):
        print(f"    {origin:22s} {n:4d}")
    print("    A naming convention, not a sequence measurement. Stage 5 records")
    print("    humanisation_state NOT_CONNECTED from the deposited taxonomy; this")
    print("    is a second, equally indirect reading reported beside it.")

    print()
    print(f"  configuration hash "
          f"{stage9.configuration_hash(stage6_hash, [g.gene for g in gated], ceiling)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
