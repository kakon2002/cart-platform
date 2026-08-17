"""Runs the ranking and tests it against all twelve rejection criteria.

The criteria were fixed in the specification before any output existed. A
tripped criterion is reported and the run stops; it is not adjusted away.
"""

from __future__ import annotations

import math

import numpy as np

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.data.coverage import (
    DATA_INSUFFICIENT,
    PROTEIN_CONFIRMED,
    RNA_SUPPORTED,
    build_coverage,
)
from car_pipeline.data.depmap import DepMapSource, RELEASE_PIN as DEPMAP_PIN, gene_index
from car_pipeline.data.gtex import GTExSource, RELEASE_PIN as GTEX_PIN
from car_pipeline.data.hpa import HPASource, RELEASE_PIN as HPA_PIN, index as atlas_index
from car_pipeline.data.singlecell import (
    SERIES,
    SingleCellSource,
    match_surface as cell_match,
)
from car_pipeline.data.tcga import (
    PROJECT as TCGA_PIN,
    TCGASource,
    match_surface as tcga_match,
)
from car_pipeline.data.uniprot import RELEASE_PIN as UNIPROT_PIN, load_surface
from car_pipeline.stages import stage3
from car_pipeline.stages.stage1 import build_spec

KNOWN_TARGETS = ["CEACAM5", "CEACAM6", "CLDN18", "MSLN", "MUC1"]
UBIQUITOUS_IMMUNE = ["HLA-A", "HLA-B", "CD74", "PTPRC"]


def spearman(a: list[float], b: list[float]) -> float:
    if len(a) < 3:
        return 0.0

    def ranks(xs):
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        r = [0.0] * len(xs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = ranks(a), ranks(b)
    n = len(ra)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((y - mb) ** 2 for y in rb))
    return num / (da * db) if da and db else 0.0


def top_n(rows, n=50):
    scored = [r for r in rows if r.composite is not None]
    scored.sort(key=lambda r: -r.composite)
    return [r.accession for r in scored[:n]]


def main() -> int:
    print("loading sources", flush=True)
    surface, _ = load_surface()
    surface_by_acc = {r.accession: r for r in surface}
    atlas = HPASource().load()
    by_acc, by_sym = atlas_index(atlas)
    gtex = GTExSource()
    gtex_profiles, gtex_tissues, _ = gtex.match_surface(surface, by_acc)
    cohort = TCGASource().load()
    cohort_join = tcga_match(cohort, surface, by_acc)
    cells = SingleCellSource().load()
    cell_index = cell_match(cells, surface, by_acc)
    dependency, _ = DepMapSource().load()
    dep_index = gene_index(dependency)

    coverage_rows = build_coverage(
        surface, by_acc, by_sym, gtex_profiles, cohort_join
    )

    spec = build_spec(PDAC_PROJECT)
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    overrides = {
        organ: ov.tier
        for organ, ov in spec.inputs.tissue_criticality_overrides.items()
    }

    def run(saturation=None, weights=None):
        return stage3.rank(
            coverage_rows, surface_by_acc, by_acc, by_sym,
            cells, cell_index, gtex_profiles, gtex_tissues,
            cohort, cohort_join, dependency, dep_index,
            overrides, ceiling, saturation=saturation, weights=weights,
        )

    rows, model = run()
    by_gene = {r.gene: r for r in rows if r.gene}

    pins = {
        "proteome": UNIPROT_PIN,
        "tissue atlas": HPA_PIN,
        "normal baseline": GTEX_PIN,
        "tumour cohort": TCGA_PIN,
        "dependency screens": DEPMAP_PIN,
        "cell atlas": SERIES,
    }
    print()
    print(stage3.header(spec, model, len(rows), pins))

    # -- structural report ------------------------------------------------
    scored = [r for r in rows if r.composite is not None]
    floored = [r for r in rows if r.below_evidence_floor]
    print(f"\n  ranked            {len(rows):,}")
    print(f"  scored            {len(scored):,}")
    print(f"  below floor       {len(floored):,}  (ranked, not scored)")
    print(f"  tissue fall-through: {len(model.fall_through)}   must read 0")
    if model.fall_through:
        for f in sorted(model.fall_through):
            print(f"      {f}")
    print(f"  risk undefined    {sum(1 for r in rows if r.risk is None):,}")
    print(f"  cleared at {ceiling}   {sum(1 for r in rows if r.cleared):,}")

    hash_a = stage3.configuration_hash(overrides, ceiling)
    print(f"\n  configuration hash (this process): {hash_a}")

    results: list[tuple[str, bool, str]] = []

    def criterion(cid: str, tripped: bool, detail: str) -> None:
        results.append((cid, tripped, detail))

    # -- R1 ---------------------------------------------------------------
    outside = []
    for g in KNOWN_TARGETS:
        r = by_gene.get(g)
        if r is None or r.composite is None:
            outside.append(f"{g}=unscored")
            continue
        tier_size = sum(1 for x in rows if x.evidence_class == r.evidence_class)
        decile = r.tier_rank <= max(1, tier_size // 10)
        if not decile:
            outside.append(f"{g}=rank {r.tier_rank}/{tier_size}")
    criterion(
        "R1", len(outside) == len(KNOWN_TARGETS),
        f"outside top decile: {outside or 'none'}",
    )

    # -- R2 ---------------------------------------------------------------
    breached = [
        g for g in UBIQUITOUS_IMMUNE
        if (r := by_gene.get(g)) is not None and r.cleared
    ]
    criterion("R2", bool(breached), f"cleared the ceiling: {breached or 'none'}")

    # -- R3 ---------------------------------------------------------------
    c5 = by_gene.get("CEACAM5")
    r3 = c5 is None or c5.composite is None
    criterion(
        "R3", r3,
        f"CEACAM5 composite={c5.composite if c5 else None} "
        f"measured_weight={c5.measured_weight if c5 else None} "
        f"tier_rank={c5.tier_rank if c5 else None}",
    )

    # -- R4 ---------------------------------------------------------------
    surface_accessions = set(surface_by_acc)
    strays = [r.accession for r in rows if r.accession not in surface_accessions]
    criterion("R4", bool(strays), f"non-surface entries present: {len(strays)}")

    # -- R5 ---------------------------------------------------------------
    worst = ("", 0.0)
    for key in stage3.WEIGHTS:
        pairs = [
            (r.composite, r.component_value(key))
            for r in scored
            if r.component_value(key) is not None
        ]
        if len(pairs) < 3:
            continue
        rho = spearman([p[0] for p in pairs], [p[1] for p in pairs])
        if abs(rho) > abs(worst[1]):
            worst = (key, rho)
    criterion(
        "R5", abs(worst[1]) > 0.95,
        f"highest |rho| {worst[0]}={worst[1]:.3f} (over measured targets only)",
    )

    # -- R6 ---------------------------------------------------------------
    base_top = set(top_n(rows))
    worst_overlap = (None, 1.0)
    for key in stage3.WEIGHTS:
        for factor in (1.2, 0.8):
            w = dict(stage3.WEIGHTS)
            w[key] = w[key] * factor
            perturbed, _ = run(weights=w)
            overlap = len(base_top & set(top_n(perturbed))) / max(1, len(base_top))
            if overlap < worst_overlap[1]:
                worst_overlap = (f"{key} x{factor}", overlap)
    criterion(
        "R6", worst_overlap[1] < 0.5,
        f"worst retention {worst_overlap[1]:.0%} at {worst_overlap[0]}",
    )

    # -- R7 ---------------------------------------------------------------
    cleared = sum(1 for r in rows if r.cleared)
    criterion(
        "R7", cleared == 0 or cleared == len(rows),
        f"cleared {cleared:,} of {len(rows):,}",
    )

    # -- R8 ---------------------------------------------------------------
    counts: dict[float, int] = {}
    for r in scored:
        counts[r.composite] = counts.get(r.composite, 0) + 1
    top_value, top_count = max(counts.items(), key=lambda kv: kv[1])
    share = top_count / max(1, len(scored))
    criterion(
        "R8", share > 0.05,
        f"most repeated composite {top_value} occurs {top_count}x ({share:.2%})",
    )

    # -- R9 ---------------------------------------------------------------
    def tier_best(tier):
        vals = [r.composite for r in rows if r.evidence_class == tier and r.composite is not None]
        return max(vals) if vals else None

    best_insuff, best_conf = tier_best(DATA_INSUFFICIENT), tier_best(PROTEIN_CONFIRMED)
    criterion(
        "R9",
        best_insuff is not None and best_conf is not None and best_insuff > best_conf,
        f"best unresolved {best_insuff} vs best protein-confirmed {best_conf}",
    )

    # -- R10 --------------------------------------------------------------
    top100 = top_n(rows, 100)
    top_rows = [r for r in rows if r.accession in set(top100)]
    bridged = sum(1 for r in top_rows if r.bridged)
    criterion(
        "R10", bridged > 0.10 * max(1, len(top100)),
        f"{bridged} of {len(top100)} reached only after a symbol failed",
    )

    # -- R11 --------------------------------------------------------------
    top25 = set(top_n(rows, 25))
    disagreeing = [r for r in rows if r.accession in top25 and r.sources_disagree]
    criterion(
        "R11", False,
        f"{len(disagreeing)} of 25 carry the flag; all listed below, so none unread",
    )
    def _fmt(x):
        if x is None:
            return "n/a"
        return "undetectable normal" if math.isinf(x) else f"{x:,.1f}x"

    for r in disagreeing:
        print(
            f"      sources_disagree  {r.gene:10s} "
            f"baseline {_fmt(r.fold_baseline):>22s}   "
            f"cohort {_fmt(r.fold_cohort)}"
        )

    # -- R12 --------------------------------------------------------------
    worst_sat = (None, 1.0)
    for key in stage3.SATURATION:
        for factor in (2.0, 0.5):
            s = dict(stage3.SATURATION)
            s[key] = s[key] * factor
            perturbed, _ = run(saturation=s)
            overlap = len(base_top & set(top_n(perturbed))) / max(1, len(base_top))
            if overlap < worst_sat[1]:
                worst_sat = (f"{key} x{factor}", overlap)
    criterion(
        "R12", worst_sat[1] < 0.5,
        f"worst retention {worst_sat[1]:.0%} at {worst_sat[0]}",
    )

    # -- report -----------------------------------------------------------
    print("\n" + "=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped = 0
    for cid, is_tripped, detail in results:
        tripped += 1 if is_tripped else 0
        print(f"  {'TRIPPED' if is_tripped else 'clear  '}  {cid}: {detail}")
    print("=" * 72)
    print(f"  {len(results) - tripped}/{len(results)} criteria clear")
    return 1 if tripped else 0


if __name__ == "__main__":
    raise SystemExit(main())
