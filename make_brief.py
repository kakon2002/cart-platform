"""Emit the comparison tables behind reports/msln-cldn18-brief.md."""

from __future__ import annotations

import json
import sys

import numpy as np

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.genespan import GeneSpanSource
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import (
    PRIMARY_TUMOUR,
    TCGASource,
    match_surface as tcga_match,
)
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import stage3, stage4
from car_pipeline.stages.stage1 import build_spec

NAMED = ("MSLN", "CLDN18")
PIPELINE_PAIR = ("NPSR1", "PTPRN2")
TOP_N = 20


def load():
    """Everything the pairing stage needs, loaded as its verifier loads it."""
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
    coverage_rows = build_coverage(
        surface, by_acc, by_sym, gtex_profiles, cohort_join)

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

    pool = stage4.build_pool(rows)
    per_organ = {}
    for r in pool:
        entry = by_acc.get(r.accession) or by_sym.get(r.gene)
        profile = gtex_profiles.get(r.accession)
        per_organ[r.accession] = stage3.per_organ_scores(
            model, entry,
            profile.values if profile is not None else None,
            gtex_tissues, calibration,
        )

    print("building the per cell matrix", flush=True)
    cells = cell_source.load_malignant(sorted({r.gene for r in pool}))
    pairs = stage4.evaluate(pool, per_organ, model, ceiling, cells)
    try:
        spans = GeneSpanSource().load()
    except Exception:                                          # noqa: BLE001
        spans = {}
    annotated = stage4.annotate_span_context(pairs, spans)
    print(f"  {len(pairs):,} pairs, span context on {annotated:,}")

    primary = cohort.sample_types == PRIMARY_TUMOUR
    tumour_tpm: dict[str, float] = {}
    for r in pool:
        cj = cohort_join.get(r.accession)
        if cj is None:
            continue
        value = float(np.median(cohort.values[primary, cj[0]]))
        if np.isfinite(value):
            tumour_tpm[r.gene] = value
    return pool, pairs, ceiling, tumour_tpm


def _partner(p, gene):
    """The far side of a pair."""
    return p.gene_b if p.gene_a == gene else p.gene_a


def cell(value, fmt="%.4f"):
    """A number, or the words that say it was not computed."""
    return "not computed" if value is None else fmt % value


def row_for(p, ceiling):
    """One comparison row on the columns the brief already uses."""
    cov = p.coverage
    return {
        "pair": f"{p.gene_a} + {p.gene_b}",
        "coexpression": None if cov.f_ab is None else 100.0 * cov.f_ab,
        "span_kb": cov.span_geomean_kb,
        "span_pct": None if cov.span_percentile is None
        else 100.0 * cov.span_percentile,
        "risk_a": p.risk_a,
        "risk_b": p.risk_b,
        "combined": p.risk.combined,
        "cleared": p.cleared,
        "unresolved": len(p.risk.unresolved_organs or []),
        "confidence": p.confidence,
        "organs_resolved": p.organs_resolved,
        "organs_total": p.organs_total,
        "rescued": ",".join(p.rescued) or "-",
    }


def table(title, rows):
    """Print one table."""
    print()
    print("=" * 118)
    print(title)
    print("=" * 118)
    print("%-26s %11s %8s %9s %19s %11s %8s %6s"
          % ("pair", "co-expr", "span kb", "span %ile",
             "member risks", "combined", "s4 rec", "unres"))
    for r in rows:
        risks = "%s / %s" % (cell(r["risk_a"]), cell(r["risk_b"]))
        print("%-26s %10s %8s %9s %19s %11s %8s %6s"
              % (r["pair"],
                 cell(r["coexpression"], "%.2f%%"),
                 cell(r["span_kb"], "%.0f"),
                 cell(r["span_pct"], "%.0f"),
                 risks,
                 cell(r["combined"]),
                 "yes" if r.get("stage4_would_recommend") else "no",
                 r["unresolved"]))


def main() -> int:
    """Build every table the brief needs."""
    pool, pairs, ceiling, tumour_tpm = load()
    by_gene = {r.gene: r for r in pool}
    blocked = {g for g, r in by_gene.items() if not r.cleared}

    def recommendable(pair, gene):
        """Stage 4's own rule: measured coverage and an eligible partner."""
        return (
            pair.coverage.measured
            and tumour_tpm.get(_partner(pair, gene), 0.0)
            >= stage4.PARTNER_MIN_TUMOUR_TPM
        )

    best_for: dict[str, object] = {}
    best_eligible: dict[str, object] = {}
    for p in pairs:
        for g in p.rescued:
            key = (p.risk.combined, _partner(p, g))
            cur = best_for.get(g)
            if cur is None or key < (cur.risk.combined, _partner(cur, g)):
                best_for[g] = p
            if recommendable(p, g):
                cur_e = best_eligible.get(g)
                if cur_e is None or key < (
                        cur_e.risk.combined, _partner(cur_e, g)):
                    best_eligible[g] = p

    rescued = sorted(best_for)
    measured = sum(1 for g in rescued if best_for[g].coverage.measured)
    print()
    print(f"pool {len(pool)}, blocked alone {len(blocked)}, "
          f"rescued by some cleared pair {len(rescued)}, "
          f"of which {measured} have measured co-expression")

    named = next(
        (p for p in pairs
         if {p.gene_a, p.gene_b} == set(NAMED)), None)
    pipeline_pair = next(
        (p for p in pairs
         if {p.gene_a, p.gene_b} == set(PIPELINE_PAIR)), None)

    ranked = sorted(
        rescued,
        key=lambda g: (best_for[g].risk.combined, _partner(best_for[g], g)),
    )

    rows = []
    for g in ranked:
        best = best_for[g]
        partner = _partner(best, g)
        entry = row_for(best, ceiling)
        entry["pair"] = g + " + " + partner
        entry["rescued_target"] = g
        entry["partner_tumour_tpm"] = tumour_tpm.get(partner)
        entry["stage4_would_recommend"] = bool(
            best.coverage.measured
            and tumour_tpm.get(partner, 0.0) >= stage4.PARTNER_MIN_TUMOUR_TPM
        )
        rows.append(entry)
    table(f"PAIRING-VIABLE SET — {len(rescued)} targets that cannot clear alone "
          f"and can clear paired; top {min(TOP_N, len(rows))} by combined risk",
          rows[:TOP_N])

    elig_rows = []
    for g in sorted(
            best_eligible,
            key=lambda x: (best_eligible[x].risk.combined,
                           _partner(best_eligible[x], x))):
        best = best_eligible[g]
        partner = _partner(best, g)
        entry = row_for(best, ceiling)
        entry["pair"] = g + " + " + partner
        entry["rescued_target"] = g
        entry["partner_tumour_tpm"] = tumour_tpm.get(partner)
        entry["stage4_would_recommend"] = True
        elig_rows.append(entry)
    table(f"THE SAME SET, RESTRICTED TO PAIRS STAGE 4 WOULD RECOMMEND - "
          f"{len(elig_rows)} of {len(rescued)} have a partner above the "
          f"{stage4.PARTNER_MIN_TUMOUR_TPM:.0f} TPM eligibility floor",
          elig_rows[:TOP_N])

    extra = []
    for p_obj in (pipeline_pair, named):
        if p_obj is None:
            continue
        e = row_for(p_obj, ceiling)
        e["stage4_would_recommend"] = bool(
            p_obj.cleared and p_obj.coverage.measured)
        extra.append(e)
    table("THE TWO PAIRS UNDER DISCUSSION", extra)

    print()
    print("=" * 118)
    print("EVIDENCE FIELDS NOW PERSISTED")
    print("=" * 118)
    for g in NAMED + PIPELINE_PAIR:
        r = by_gene.get(g)
        if r is None:
            print(f"  {g:<10} not in the pool")
            continue
        measured = sum(1 for c in r.components.values() if c.measured)
        print(f"  {g:<10} evidence_class={r.evidence_class} "
              f"measured_weight={r.measured_weight} "
              f"components={measured}/{len(r.components)} "
              f"target_confidence={r.confidence} "
              f"risk={r.risk} organ={r.risk_organ} cleared={r.cleared}")
    for label, p in (("MSLN+CLDN18", named), ("NPSR1+PTPRN2", pipeline_pair)):
        if p is None:
            print(f"  {label:<14} pair not present")
            continue
        print(f"  {label:<14} pair_confidence={p.confidence} "
              f"organs_resolved={p.organs_resolved}/{p.organs_total} "
              f"coverage_measured={p.coverage.measured} "
              f"combined={p.risk.combined} cleared={p.cleared} "
              f"unresolved_organs={len(p.risk.unresolved_organs or [])}")

    out = {
        "pool": len(pool),
        "blocked_alone": len(blocked),
        "rescued": len(rescued),
        "rescued_genes": rescued,
        "ranked": rows,
        "ranked_eligible": elig_rows,
        "named_pair": row_for(named, ceiling) if named else None,
        "pipeline_pair": (
            row_for(pipeline_pair, ceiling) if pipeline_pair else None),
    }
    path = "reports/brief-tables.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print()
    print(f"written to {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
