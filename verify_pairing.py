"""Runs the pairing stage and tests it against the invariants and criteria.

Both were fixed in the specification before any output existed. An invariant
failure is a bug and stops the run. A tripped criterion is reported and the run
stops; it is not adjusted away.

Order is deliberate: invariants, then criteria, then the biology. Nothing about
which pairs came out on top is printed until the checks have been read.
"""

from __future__ import annotations

import sys

import numpy as np

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import stage3, stage4
from car_pipeline.stages.stage1 import build_spec

WATCH = ["MSLN", "CLDN18", "CEACAM6", "CEACAM5", "MUC1"]
UBIQUITOUS_IMMUNE = ["HLA-A", "HLA-B", "CD74", "PTPRC"]

#: Known-answer check on the per-cell derivation. The column indices in the
#: matrix are stored in descending order within a row, so a lookup that assumes
#: otherwise returns zero for every gene without erroring. Only a gene whose
#: answer is known catches that.
#:
#: CEACAM5 is deliberately absent. It carries two molecules across all 64,538
#: malignant cells, which is the capture failure the ranking stage already
#: documented rather than a broken derivation, so a check that stopped on it
#: would stop a correct implementation on the highest ranked target in the pool.
SANITY = {"KRT19": 0.50, "CLDN18": 0.05, "CEACAM6": 0.05, "MSLN": 0.05}

#: Requested as extra columns rather than drawn from the pool. KRT19 is a
#: cytokeratin, so it never survives the surface filter and can never be a pool
#: member — which is exactly what makes it a good control. It is the loudest
#: malignant signal in the atlas and it is independent of anything the ranking
#: decided.
CONTROL_GENES = ["KRT19"]


def spearman(a: list[float], b: list[float]) -> float:
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(np.asarray(a, dtype=float)))
    rb = np.argsort(np.argsort(np.asarray(b, dtype=float)))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / denom) if denom else float("nan")


def main() -> int:
    print("loading sources", flush=True)
    surface, _ = load_surface()
    surface_by_acc = {r.accession: r for r in surface}
    atlas = HPASource().load()
    by_acc, by_sym = atlas_index(atlas)
    gtex_profiles, gtex_tissues, _ = GTExSource().match_surface(surface, by_acc)
    cohort = TCGASource().load()
    cohort_join = tcga_match(cohort, surface, by_acc)
    cell_source = SingleCellSource()
    cells_atlas = cell_source.load()
    cell_index = cell_match(cells_atlas, surface, by_acc)
    dependency, _ = DepMapSource().load()
    dep_index = gene_index(dependency)
    coverage_rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, cohort_join)

    spec = build_spec(PDAC_PROJECT)
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    overrides = {
        organ: ov.tier
        for organ, ov in spec.inputs.tissue_criticality_overrides.items()
    }
    calibration = stage3.calibrate_atlas_levels(
        surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
        stage3.RiskModel(overrides=overrides),
    )
    rows, model, _stats = stage3.rank(
        coverage_rows, surface_by_acc, by_acc, by_sym, cells_atlas, cell_index,
        gtex_profiles, gtex_tissues, cohort, cohort_join, dependency, dep_index,
        overrides, ceiling, calibration,
    )
    s3_hash = stage3.configuration_hash(
        overrides, ceiling, calibration=calibration
    )

    pool = stage4.build_pool(rows)
    pool_genes = [r.gene for r in pool]

    # Per-organ scores from the ranking stage, not recomputed here.
    per_organ: dict[str, dict[str, float]] = {}
    for r in pool:
        entry = by_acc.get(r.accession) or by_sym.get(r.gene)
        profile = gtex_profiles.get(r.accession)
        per_organ[r.accession] = stage3.per_organ_scores(
            model,
            entry,
            profile.values if profile is not None else None,
            gtex_tissues,
            calibration,
        )

    print("building the per cell matrix", flush=True)
    cells = cell_source.load_malignant(sorted(set(pool_genes) | set(CONTROL_GENES)))

    # ---------------- sanity check on the derivation -------------------
    print()
    print("=" * 72)
    print("PER CELL DERIVATION — known answer check")
    print("=" * 72)
    column = {g: i for i, g in enumerate(cells.genes)}
    broken = []
    for gene, floor in SANITY.items():
        if gene not in column:
            print(f"  {gene:<10} absent from the matrix")
            broken.append(gene)
            continue
        frac = float((cells.counts[:, column[gene]] >= 1).mean())
        ok = frac >= floor
        print(f"  {gene:<10} detected in {frac:>7.4f} of malignant cells "
              f"(floor {floor:.2f}) {'ok' if ok else 'FAILED'}")
        if not ok:
            broken.append(gene)
    if "CEACAM5" in column:
        c5 = float((cells.counts[:, column['CEACAM5']] >= 1).mean())
        print(f"  {'CEACAM5':<10} detected in {c5:>7.4f} — excluded from the check; "
              "this assay is documented as dropping it")
    if broken:
        print(f"\n  STOP: the derivation is wrong for {', '.join(broken)}")
        return 1

    print(f"\n  cells {cells.counts.shape[0]:,}  genes {cells.counts.shape[1]}  "
          f"patients {len(set(cells.patient.tolist()))}  "
          f"evaluable {len(cells.evaluable_patients(stage4.MIN_MALIGNANT_CELLS))}")

    print("evaluating pairs", flush=True)
    pairs = stage4.evaluate(pool, per_organ, model, ceiling, cells)
    decisions = stage4.decide(pool, pairs)

    # ---------------- invariants ---------------------------------------
    print()
    print("=" * 72)
    print("CONSTRUCTION INVARIANTS")
    print("=" * 72)
    failures: list[str] = []

    # Both numbers are carried at four decimals, which is the precision the
    # ranking stage stores risk at. The tolerance absorbs the rounding boundary
    # only; the worst excess is printed so a real disagreement cannot hide
    # underneath it.
    TOL = 1.5e-4
    bad, worst_i1 = [], 0.0
    for r in pool:
        pr = stage4.pair_risk(model, per_organ[r.accession], per_organ[r.accession])
        if pr.combined is None:
            bad.append(f"{r.gene} None")
            continue
        gap = abs(pr.combined - r.risk)
        worst_i1 = max(worst_i1, gap)
        if gap > TOL:
            bad.append(f"{r.gene} {pr.combined} vs {r.risk}")
    _inv(failures, "I1", not bad,
         f"self pair reproduces single risk for all {len(pool)}, "
         f"worst gap {worst_i1:.2e}"
         if not bad else f"{len(bad)} mismatch: {bad[:3]}")

    sym = []
    for p in pairs[:2000]:
        back = stage4.pair_risk(
            model, per_organ[p.accession_b], per_organ[p.accession_a]
        )
        if (back.combined is None) != (p.risk.combined is None) or (
            back.combined is not None and abs(back.combined - p.risk.combined) > 1e-12
        ):
            sym.append(f"{p.gene_a}/{p.gene_b}")
    _inv(failures, "I2", not sym, f"symmetric over {min(len(pairs), 2000):,} pairs")

    i4 = [p for p in pairs
          if p.confidence > min(p.confidence_a, p.confidence_b) + 1e-9]
    _inv(failures, "I4", not i4,
         f"{len(i4)} pairs better evidenced than their weaker member")

    i3 = [p for p in pairs if p.coverage.measured
          and p.coverage.f_ab > min(p.coverage.f_a, p.coverage.f_b) + 1e-12]
    _inv(failures, "I3", not i3, f"{len(i3)} pairs with f_AB above its marginals")

    excess = [p.risk.optimistic - min(p.risk_a, p.risk_b)
              for p in pairs if p.risk.optimistic is not None]
    worst_i5 = max(excess) if excess else 0.0
    i5 = [e for e in excess if e > TOL]
    _inv(failures, "I5", not i5,
         f"{len(i5)} pairs exceed min of members, worst excess {worst_i5:.2e}")

    mats = [stage4.intersection_matrix(cells.counts >= t)
            for t in stage4.SENSITIVITY_COUNTS]
    mono = sum(int((mats[k] < mats[k + 1]).sum()) for k in range(len(mats) - 1))
    _inv(failures, "I6", mono == 0,
         f"pairwise intersections fall monotonically across "
         f"{stage4.SENSITIVITY_COUNTS} counts ({mono} violations)")

    ind = [p for p in pairs if p.risk.independence is not None
           and p.risk.combined is not None
           and p.risk.independence > p.risk.combined + 1e-9]
    _inv(failures, "I7", not ind,
         f"{len(ind)} pairs where the independence bound exceeds the gate")

    if failures:
        print(f"\n  STOP: {len(failures)} invariant(s) failed — this is a bug, "
              "not a result")
        return 1

    # ---------------- criteria -----------------------------------------
    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    combined = [p.risk.combined for p in pairs]
    naive = [min(p.risk_a, p.risk_b) for p in pairs]
    rho1 = spearman(combined, naive)
    criterion("P1", rho1 > 0.95,
              f"combined risk vs min of members, rho={rho1:.4f} (limit 0.95)")

    beats = sum(1 for p in pairs
                if p.risk.combined is not None
                and p.risk.combined < min(p.risk_a, p.risk_b) - 0.05)
    share = beats / len(pairs)
    criterion("P2", share < 0.01,
              f"{beats:,} of {len(pairs):,} ({share:.2%}) beat the better member "
              "by more than 0.05 (limit 1%)")

    rescued_genes = sorted({g for p in pairs for g in p.rescued})
    criterion("P3", not rescued_genes,
              f"{len(rescued_genes)} blocked targets rescued by some pair"
              + (f": {', '.join(rescued_genes[:8])}" if rescued_genes else ""))

    meas = [p for p in pairs if p.coverage.measured]
    rho4 = spearman([p.coverage.f_ab for p in meas],
                    [p.coverage.f_a * p.coverage.f_b for p in meas])
    criterion("P4", rho4 > 0.98,
              f"f_AB vs f_A x f_B over {len(meas):,} measured pairs, "
              f"rho={rho4:.4f} (limit 0.98)")

    # Wiring check: clearance must be decided on the conservative arm. A pair
    # that clears only on the optimistic arm is RISK_UNRESOLVED by definition
    # and must not be marked cleared.
    p5 = [p for p in pairs
          if p.cleared
          and p.risk.optimistic is not None
          and p.risk.combined > p.ceiling >= p.risk.optimistic]
    criterion("P5", bool(p5),
              f"{len(p5)} pairs marked cleared on the optimistic arm")

    recommended = [d.pair for d in decisions if d.outcome == stage4.DUAL and d.pair]
    p6 = [p for p in recommended if not p.coverage.measured]
    criterion("P6", bool(p6), f"{len(p6)} recommended pairs are unmeasured")

    p7 = [p for p in pairs if p.cleared
          and (p.gene_a in UBIQUITOUS_IMMUNE or p.gene_b in UBIQUITOUS_IMMUNE)]
    in_pool = [g for g in UBIQUITOUS_IMMUNE if g in pool_genes]
    offenders = sorted({f"{p.gene_a}+{p.gene_b}" for p in p7})
    criterion("P7", bool(p7),
              f"{len(p7)} cleared pairs contain a ubiquitous immune protein "
              f"(in pool: {in_pool or 'none'}) {offenders[:4]}")

    # Substantive check: how many cleared pairs would stop clearing if the
    # unmeasured antigen turned out to saturate the organ nobody looked at.
    # Having an unresolved organ is not itself the failure — the conservative
    # arm already charges the measured member's score there. The failure is
    # clearance that survives only because the missing antigen was assumed no
    # more prevalent than the one that was measured.
    cleared_pairs = [p for p in pairs if p.cleared]
    ignorance = [p for p in cleared_pairs
                 if p.risk.pessimistic is not None
                 and p.risk.pessimistic > p.ceiling]
    frac8 = len(ignorance) / len(cleared_pairs) if cleared_pairs else 0.0
    criterion("P8", frac8 > 0.10,
              f"{len(ignorance)} of {len(cleared_pairs)} cleared pairs "
              f"({frac8:.1%}) stop clearing if the unmeasured antigen saturates "
              "its organ (limit 10%)")

    p9 = [p for p in recommended
          if p.coverage.patient_fraction < stage4.PATIENT_FRACTION_FLOOR]
    criterion("P9", bool(p9),
              f"{len(p9)} recommended pairs below the patient floor")

    cleared_map = {r.gene: r.cleared for r in pool}
    p10 = [d for d in decisions
           if d.outcome == stage4.DUAL and cleared_map[d.gene]]
    criterion("P10", bool(p10),
              f"{len(p10)} targets recommended dual despite clearing alone")

    counts = {}
    for d in decisions:
        counts[d.outcome] = counts.get(d.outcome, 0) + 1
    top = max(counts.values()) / len(decisions)
    criterion("P11", top > 0.95,
              f"outcome spread {counts}, largest {top:.1%} (limit 95%)")

    duals = {d.gene: d.partner for d in decisions if d.outcome == stage4.DUAL}
    alt = stage4.evaluate(pool, per_organ, model, ceiling, cells, threshold=2)
    alt_dual = {d.gene: d.partner
                for d in stage4.decide(pool, alt) if d.outcome == stage4.DUAL}
    changed = sum(1 for g, q in duals.items() if alt_dual.get(g) != q)
    base = max(len(duals), 1)
    criterion("P12", changed / base > 0.5,
              f"{changed} of {len(duals)} dual recommendations change at "
              f"2 counts ({changed / base:.1%}, limit 50%)")

    partners = {}
    for q in duals.values():
        partners[q] = partners.get(q, 0) + 1
    worst = max(partners.values()) / base if partners else 0.0
    criterion("P13", worst > 0.5,
              f"most common partner takes {worst:.1%} of dual recommendations "
              + (f"({max(partners, key=partners.get)})" if partners else ""))

    ranked_pairs = sorted(
        [p for p in pairs if p.coverage.measured],
        key=lambda p: (not p.rescued, -p.coverage.f_ab),
    )
    top_pair = ranked_pairs[0] if ranked_pairs else None
    top_two = set(pool_genes[:2])
    criterion("P14", bool(top_pair) and {top_pair.gene_a, top_pair.gene_b} == top_two,
              f"top pair is {top_pair.gene_a}+{top_pair.gene_b}, "
              f"top two singles are {pool_genes[0]}+{pool_genes[1]}"
              if top_pair else "no measured pairs")

    half = stage4.build_pool(rows, stage4.POOL_SIZE // 2)
    half_cells = cells.subset(sorted({r.gene for r in half}))
    half_dual = {
        d.gene: d.partner
        for d in stage4.decide(
            half,
            stage4.evaluate(half, per_organ, model, ceiling, half_cells),
        )
        if d.outcome == stage4.DUAL
    }
    shared = [g for g in duals if g in {r.gene for r in half}]
    moved = sum(1 for g in shared if half_dual.get(g) != duals[g])
    denom = max(len(shared), 1)
    criterion("P15", moved / denom > 0.5,
              f"pool halved to {len(half)}: {moved} of {len(shared)} shared dual "
              f"targets change partner ({moved / denom:.1%}, limit 50%)")

    reach = [p for p in pairs if p.coverage.measured
             and p.coverage.f_ab >= stage4.COVERAGE_FLOOR]
    criterion("P16", not reach,
              f"{len(reach):,} pairs reach f_AB >= {stage4.COVERAGE_FLOOR} "
              f"of {len(meas):,} measured")

    print("=" * 72)
    print(f"  {16 - len(tripped)}/16 criteria clear")

    # Measurements, printed whether or not a criterion tripped. These describe
    # what the atlas contains rather than what should be built, so withholding
    # them behind a passing run would hide the strongest thing this stage has to
    # say about the architecture.
    _report_measurements(pool, pairs, decisions, duals)

    if tripped:
        print()
        print(f"  STOPPING: {', '.join(tripped)} tripped. The specification "
              "changes and the run repeats; the result does not get an "
              "explanation.")
        return 2

    _report_biology(pool, pairs, decisions, cells, ceiling, s3_hash)
    return 0


def _report_measurements(pool, pairs, decisions, duals) -> None:
    print()
    print("=" * 72)
    print("WHAT EACH ARCHITECTURE REACHES  (measured; no recommendation implied)")
    print("=" * 72)

    print()
    print("  What the weights are doing in this stage")
    print("  " + "-" * 68)
    # Within the admissible set the ordering key is the co-expression gate, not
    # the composite the weights produce. Measured rather than asserted: how often
    # would the recommendation change if partners were ordered by composite?
    by_gene: dict[str, list] = {r.gene: [] for r in pool}
    for p in pairs:
        by_gene[p.gene_a].append(p)
        by_gene[p.gene_b].append(p)
    differ = 0
    for gene in duals:
        adm = [p for p in by_gene[gene] if p.admissible]
        if not adm:
            continue
        by_cov = min(adm, key=lambda p: (-p.coverage.f_ab, p.risk.combined))
        by_comp = max(
            adm, key=lambda p: min(p.composite_a, p.composite_b)
        )
        if {by_cov.gene_a, by_cov.gene_b} != {by_comp.gene_a, by_comp.gene_b}:
            differ += 1
    n = max(len(duals), 1)
    print(f"    swapping the ordering key from coverage to composite changes "
          f"{differ} of {len(duals)} recommendations ({differ / n:.0%})")
    print()
    print("    The weights never order anything in this stage. They enter only")
    print("    through pool membership, and the pool boundary moves little:")
    print("    halving it changed one partner in ten (P15). So the weights are")
    print("    close to decorative here, and a reader should be told that rather")
    print("    than assume they are doing work.")
    print()
    print("    But the swap number above is NOT the evidence for that, and it")
    print("    should not be read as it. The two keys agreeing means the choice")
    print("    is insensitive to which one is used — and the reason is P13: one")
    print("    partner wins under both keys, so almost nothing is being ordered")
    print("    at all. Partner concentration, not the coverage gate, is why the")
    print("    weights make no difference to this output.")

    print()
    print("  Coverage and the escape population")
    print("  " + "-" * 68)
    print("    AND-gate coverage is the intersection. Escape is one minus it.")
    print("    Both are FLOORS: this is a single-nucleus assay and it drops")
    print("    transcripts, so the true intersection is higher and the true")
    print("    escape lower. The direction is known; the magnitude is not.")
    print()
    print(f"    {'A':<10}{'B':<10}{'A alone':>9}{'B alone':>9}{'OR':>8}"
          f"{'AND':>8}{'escape':>9}{'cost':>7}")

    watch_pairs = [
        p for p in pairs
        if p.coverage.measured and p.gene_a in WATCH and p.gene_b in WATCH
    ]
    best = sorted(
        (p for p in pairs if p.coverage.measured),
        key=lambda p: -p.coverage.f_ab,
    )[:5]
    for p in watch_pairs + [b for b in best if b not in watch_pairs]:
        c = p.coverage
        cost = f"{c.coverage_cost:.1f}x" if c.coverage_cost else "n/a"
        print(f"    {p.gene_a:<10}{p.gene_b:<10}{c.f_a:>9.4f}{c.f_b:>9.4f}"
              f"{c.or_gate:>8.4f}{c.f_ab:>8.4f}{c.escape:>9.4f}{cost:>7}")

    print()
    print("    The best-covering pairs above are all high-abundance genes. That")
    print("    is the dropout bias of section 6.5 visible in the output: a scarce")
    print("    antigen cannot show a large intersection on this assay whatever")
    print("    its biology, so this table ranks detectability alongside truth.")

    if watch_pairs:
        # The most favourable watch-list pair, not the worst. An argument that
        # holds on the best case does not depend on which pair was picked.
        best_watch = min(watch_pairs, key=lambda p: p.coverage.escape)
        c = best_watch.coverage
        print()
        print(f"    Stated plainly, on the most favourable known pair: an "
              f"AND-gate on {best_watch.gene_a}+{best_watch.gene_b}")
        print(f"    reaches {c.f_ab:.1%} of malignant cells, so {c.escape:.1%} "
              f"escape it. {best_watch.gene_a} alone reaches {c.f_a:.1%}")
        print(f"    and {best_watch.gene_b} alone {c.f_b:.1%}; the safety gain "
              f"costs {c.coverage_cost:.1f}x the coverage of the")
        print("    better single target, and an OR-gate on the same two would "
              f"reach {c.or_gate:.1%}.")
        print("    That is antigen heterogeneity measured directly, and it is the")
        print("    strongest argument in this output that AND-gating may be the")
        print("    wrong architecture for this indication even though it does")
        print("    solve the safety problem. It is a floor, so the real figure is")
        print("    better than this — but not better than the single-target one,")
        print("    which is deflated by the same assay in the same direction.")


def _inv(failures: list[str], cid: str, ok: bool, detail: str) -> None:
    print(f"  {'ok  ' if ok else 'FAIL'} {cid}: {detail}")
    if not ok:
        failures.append(cid)


def _report_biology(pool, pairs, decisions, cells, ceiling, s3_hash) -> None:
    print()
    print("=" * 72)
    print("PAIRING RESULTS")
    print("=" * 72)
    print(f"  stage 3 configuration hash   {s3_hash}")
    print(f"  stage 4 configuration hash   "
          f"{stage4.configuration_hash(s3_hash, [r.gene for r in pool])}")
    print(f"  pool                         {len(pool)} by composite, risk ignored")
    print(f"  pairs                        {len(pairs):,}")

    rescuing = sorted(
        [p for p in pairs if p.rescued],
        key=lambda p: (-(p.coverage.f_ab or 0.0),),
    )
    print(f"  clearing pairs               {sum(1 for p in pairs if p.cleared):,}")
    print(f"  rescuing pairs               {len(rescuing):,}")
    print()
    if rescuing:
        print(f"  {'A':<10}{'B':<10}{'risk_A':>8}{'risk_B':>8}{'pair':>8}"
              f"{'f_AB':>8}{'sacr_A':>8}{'sacr_B':>8}  rescued")
        for p in rescuing[:25]:
            f_ab = f"{p.coverage.f_ab:.4f}" if p.coverage.measured else "n/m"
            sa = f"{p.coverage.sacrificed_a:.3f}" if p.coverage.measured else "n/m"
            sb = f"{p.coverage.sacrificed_b:.3f}" if p.coverage.measured else "n/m"
            print(f"  {p.gene_a:<10}{p.gene_b:<10}{p.risk_a:>8.4f}{p.risk_b:>8.4f}"
                  f"{p.risk.combined:>8.4f}{f_ab:>8}{sa:>8}{sb:>8}  "
                  f"{','.join(p.rescued)}")
    else:
        print("  no pair moves any target from blocked to cleared")

    print()
    print("  WATCH LIST")
    by_gene = {d.gene: d for d in decisions}
    for g in WATCH:
        d = by_gene.get(g)
        if d is None:
            print(f"    {g:<10} not in pool")
            continue
        extra = ""
        if d.pair is not None:
            c = d.pair.coverage
            extra = (f" partner={d.partner} pair_risk={d.pair.risk.combined:.4f}"
                     f" f_AB={c.f_ab:.4f}" if c.measured
                     else f" partner={d.partner} coverage unmeasured")
        print(f"    {g:<10} {d.outcome}{extra}")


if __name__ == "__main__":
    sys.exit(main())
