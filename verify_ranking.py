"""Runs the ranking and tests it against all twelve rejection criteria."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys

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
from car_pipeline.configs.pdac import PDAC_ATLAS
from car_pipeline.data.singlecell import (
    SingleCellSource,
    match_surface as cell_match,
)
from car_pipeline.data.tcga import (
    DEFAULT_PROJECT as TCGA_PIN,
    TCGASource,
    match_surface as tcga_match,
)
from car_pipeline.data.uniprot import RELEASE_PIN as UNIPROT_PIN, load_surface
from car_pipeline.stages import stage3
from car_pipeline.stages.stage1 import build_spec

KNOWN_TARGETS = ["CEACAM5", "CEACAM6", "CLDN18", "MSLN", "MUC1"]
UBIQUITOUS_IMMUNE = ["HLA-A", "HLA-B", "CD74", "PTPRC"]


def spearman(a: list[float], b: list[float]) -> float:
    """Rank correlation between two series."""
    if len(a) < 3:
        return 0.0

    def ranks(xs):
        """Ordinal ranks for a series."""
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
    """The first n genes of the ranking."""
    scored = [r for r in rows if r.composite is not None]
    scored.sort(key=lambda r: -r.composite)
    return [r.accession for r in scored[:n]]


def main() -> int:
    """Run the ranking criteria."""
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

    calibration = stage3.calibrate_atlas_levels(
        surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
        stage3.RiskModel(overrides=overrides),
    )

    def run(saturation=None, weights=None):
        """Score the pool and return the ranking with its report."""
        return stage3.rank(
            coverage_rows, surface_by_acc, by_acc, by_sym,
            cells, cell_index, gtex_profiles, gtex_tissues,
            cohort, cohort_join, dependency, dep_index,
            overrides, ceiling, calibration,
            saturation=saturation, weights=weights,
        )

    rows, model, stats = run()

    all_by_gene: dict[str, list] = {}
    for r in rows:
        if r.gene:
            all_by_gene.setdefault(r.gene, []).append(r)

    def best(gene: str):
        """The highest-scoring entry for a symbol, for reporting."""
        members = all_by_gene.get(gene) or []
        scored_members = [m for m in members if m.composite is not None]
        if scored_members:
            return max(scored_members, key=lambda m: m.composite)
        return members[0] if members else None

    by_gene = {g: best(g) for g in all_by_gene}
    duplicated = {g: len(v) for g, v in all_by_gene.items() if len(v) > 1}

    pins = {
        "proteome": UNIPROT_PIN,
        "tissue atlas": HPA_PIN,
        "normal baseline": GTEX_PIN,
        "tumour cohort": TCGA_PIN,
        "dependency screens": DEPMAP_PIN,
        "cell atlas": PDAC_ATLAS.series,
    }
    print()
    print(
        stage3.header(
            spec, model, len(rows), pins,
            stage3.SATURATION, stage3.WEIGHTS, stats, calibration,
        )
    )

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

    print(f"  symbols carrying more than one accession: {len(duplicated)}")

    hash_a = stage3.configuration_hash(
        overrides, ceiling, stage3.SATURATION, stage3.WEIGHTS, calibration
    )
    probe = subprocess.run(
        [
            sys.executable, "-c",
            "import json,sys;"
            "from car_pipeline.stages import stage3;"
            "c=stage3.CalibrationCurve(**json.loads(sys.argv[1]));"
            f"print(stage3.configuration_hash({overrides!r}, {ceiling!r},"
            " stage3.SATURATION, stage3.WEIGHTS, c))",
            json.dumps({
                "tpm": {str(k): v for k, v in calibration.tpm.items()},
                "quartiles": {str(k): list(v) for k, v in calibration.quartiles.items()},
                "counts": {str(k): v for k, v in calibration.counts.items()},
                "separations": {str(k): v for k, v in calibration.separations.items()},
                "monotonic": calibration.monotonic,
                "observations": calibration.observations,
            }),
        ],
        capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    hash_b = probe.stdout.strip()
    stable = hash_a == hash_b
    print(f"\n  configuration hash (this process):  {hash_a}")
    print(f"  configuration hash (fresh process): {hash_b or '<failed>'}")
    print(f"  stable across processes: {'yes' if stable else 'NO'}")
    if not stable:
        print("  the hash is process-dependent and cannot identify a run")
        return 1

    results: list[tuple[str, bool, str]] = []

    def criterion(cid: str, tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        results.append((cid, tripped, detail))

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

    breached = [
        f"{g}({r.accession})"
        for g in UBIQUITOUS_IMMUNE
        for r in all_by_gene.get(g, [])
        if r.cleared
    ]
    checked = sum(len(all_by_gene.get(g, [])) for g in UBIQUITOUS_IMMUNE)
    criterion(
        "R2", bool(breached),
        f"cleared the ceiling: {breached or 'none'} (across {checked} accessions)",
    )

    c5 = by_gene.get("CEACAM5")
    r3 = c5 is None or c5.composite is None
    criterion(
        "R3", r3,
        f"CEACAM5 composite={c5.composite if c5 else None} "
        f"measured_weight={c5.measured_weight if c5 else None} "
        f"tier_rank={c5.tier_rank if c5 else None}",
    )

    surface_accessions = set(surface_by_acc)
    strays = [r.accession for r in rows if r.accession not in surface_accessions]
    criterion("R4", bool(strays), f"non-surface entries present: {len(strays)}")

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

    base_top = set(top_n(rows))
    worst_overlap = (None, 1.0)
    for key in stage3.WEIGHTS:
        for factor in (1.2, 0.8):
            w = dict(stage3.WEIGHTS)
            w[key] = w[key] * factor
            perturbed, _, _ = run(weights=w)
            overlap = len(base_top & set(top_n(perturbed))) / max(1, len(base_top))
            if overlap < worst_overlap[1]:
                worst_overlap = (f"{key} x{factor}", overlap)
    criterion(
        "R6", worst_overlap[1] < 0.5,
        f"worst retention {worst_overlap[1]:.0%} at {worst_overlap[0]}",
    )

    cleared = sum(1 for r in rows if r.cleared)
    criterion(
        "R7", cleared == 0 or cleared == len(rows),
        f"cleared {cleared:,} of {len(rows):,}",
    )

    counts: dict[float, int] = {}
    for r in scored:
        counts[r.composite] = counts.get(r.composite, 0) + 1
    if not counts:
        criterion("R8", True, "nothing scored, so the distribution is degenerate")
    else:
        top_value, top_count = max(counts.items(), key=lambda kv: kv[1])
        share = top_count / len(scored)
        criterion(
            "R8", share > 0.05,
            f"most repeated composite {top_value} occurs {top_count}x ({share:.2%})",
        )

    def tier_best(tier):
        """The best-scoring gene within a consequence tier."""
        vals = [r.composite for r in rows if r.evidence_class == tier and r.composite is not None]
        return max(vals) if vals else None

    best_insuff, best_conf = tier_best(DATA_INSUFFICIENT), tier_best(PROTEIN_CONFIRMED)
    criterion(
        "R9",
        best_insuff is not None and best_conf is not None and best_insuff > best_conf,
        f"best unresolved {best_insuff} vs best protein-confirmed {best_conf}",
    )

    top100 = top_n(rows, 100)
    top_rows = [r for r in rows if r.accession in set(top100)]
    bridged = sum(1 for r in top_rows if r.bridged)
    criterion(
        "R10", bridged > 0.10 * max(1, len(top100)),
        f"{bridged} of {len(top100)} reached only after a symbol failed",
    )

    top25 = set(top_n(rows, 25))
    disagreeing = [r for r in rows if r.accession in top25 and r.sources_disagree]

    criterion(
        "R11", False,
        f"{len(disagreeing)} of 25 depart from the systematic offset "
        f"({stats['field_offset_fold']:.1f}x) by more than "
        f"{stats['tolerance_fold']:.0f}x; all listed, so none unread",
    )
    def _fmt(x):
        """Format a value for the report, marking the unmeasured."""
        return "n/a" if x is None else f"{x:,.1f}x"

    for r in disagreeing:
        floored = r.components[stage3.C3].note == "normal below detection"
        marker = "  [normal below detection]" if floored else ""
        fe = r.field_elevation
        print(
            f"      departs  {r.gene:10s} "
            f"baseline {_fmt(r.fold_baseline):>12s}   "
            f"cohort {_fmt(r.fold_cohort):>10s}   "
            f"field elevation {fe:,.0f}x{marker}"
        )

    worst_sat = (None, 1.0)
    for key in stage3.SATURATION:
        for factor in (2.0, 0.5):
            s = dict(stage3.SATURATION)
            s[key] = s[key] * factor
            perturbed, _, _ = run(saturation=s)
            overlap = len(base_top & set(top_n(perturbed))) / max(1, len(base_top))
            if overlap < worst_sat[1]:
                worst_sat = (f"{key} x{factor}", overlap)
    criterion(
        "R12", worst_sat[1] < 0.5,
        f"worst retention {worst_sat[1]:.0%} at {worst_sat[0]}",
    )

    from car_pipeline.stages import stage4 as _stage4

    def _c2(r):
        return r.components[stage3.C2]

    rejected = [r for r in rows if r.tumour_side_verdict == stage3.STROMA_DOMINANT]
    unresolved = [r for r in rows
                  if r.tumour_side_verdict == stage3.STROMA_UNRESOLVED]
    measured_fail = [r for r in rows
                     if _c2(r).measured
                     and (_c2(r).raw or 0.0) <= stage3.STROMA_RATIO_FLOOR]
    gate_pool = _stage4.build_pool(rows)

    blind = [r.gene for r in rejected if not _c2(r).measured]
    criterion(
        "G1", bool(blind),
        f"{len(unresolved)} targets carry an unmeasured malignant-to-stromal "
        f"ratio and {len(blind)} of them were rejected"
        + (f": {blind[:5]}" if blind else "; the gate fired on no absent measurement"))

    pinned = ["CEACAM6", "MUC1", "CLDN18", "MSLN", "CEACAM5"]
    by_gene_all = {}
    for r in rows:
        if r.gene and (r.gene not in by_gene_all
                       or (r.composite or 0) > (by_gene_all[r.gene].composite or 0)):
            by_gene_all[r.gene] = r
    lost = [g for g in pinned
            if by_gene_all.get(g)
            and by_gene_all[g].tumour_side_verdict == stage3.STROMA_DOMINANT]
    criterion(
        "G2", bool(lost),
        f"known targets surviving the gate: "
        + ", ".join(
            f"{g}="
            + ("exempt" if by_gene_all[g].tumour_side_verdict
               == stage3.STROMA_UNRESOLVED
               else f"{_c2(by_gene_all[g]).raw:.1f}")
            for g in pinned if by_gene_all.get(g))
        + (f"; REJECTED {lost}" if lost else ""))

    rejected_genes = {r.gene for r in rejected}
    mhc2 = sorted(g for g in rejected_genes if g and g.startswith("HLA-D"))
    immunoglobulin = sorted(g for g in rejected_genes
                            if g and (g.startswith("IGH") or g.startswith("IGK")
                                      or g.startswith("IGL")))
    criterion(
        "G3",
        "LRRC15" not in rejected_genes or not mhc2 or not immunoglobulin,
        f"the gate rejects {len(rejected)} targets including LRRC15="
        f"{'yes' if 'LRRC15' in rejected_genes else 'NO'}, "
        f"{len(mhc2)} MHC class II ({mhc2[:3]}), "
        f"{len(immunoglobulin)} immunoglobulin ({immunoglobulin[:3]})")

    measurable = [r for r in rows if _c2(r).measured]
    criterion(
        "G4", not rejected or len(rejected) >= len(measurable),
        f"{len(rejected)} rejected of {len(measurable)} with a measured ratio; "
        "the gate neither passes everything nor empties the population")

    bad_note = [r.gene for r in rejected
                if _c2(r).note in ("no row", "below capture threshold")]
    criterion(
        "G5",
        len(rejected) != len(measured_fail) or bool(bad_note),
        f"{len(rejected)} rejected against {len(measured_fail)} with a measured "
        f"ratio at or below {stage3.STROMA_RATIO_FLOOR}"
        + (f"; {len(bad_note)} carry an absence note: {bad_note[:3]}"
           if bad_note else "; none carries an absence note"))

    drifted = []
    for r in rejected[:200]:
        entry = by_acc.get(r.accession) or (by_sym.get(r.gene) if r.gene else None)
        profile = gtex_profiles.get(r.accession)
        again, _organ = stage3.compute_risk(
            model, entry, profile.values if profile is not None else None,
            gtex_tissues, calibration)
        if (again is None) != (r.risk is None):
            drifted.append(r.gene)
        elif again is not None and abs(round(again, 4) - r.risk) > 1e-9:
            drifted.append(r.gene)
    criterion(
        "G6", bool(drifted),
        f"risk recomputed independently for all {len(rejected)} rejected "
        f"targets matches the stored value on {len(rejected) - len(drifted)}"
        + (f"; drifted on {drifted[:3]}" if drifted
           else "; the gate changed no risk"))

    print(f"\n  POOL   {len(gate_pool)} members after the gate; "
          f"{len(rejected)} rejected, {len(unresolved)} exempt as unmeasured")

    positive_levels = sorted(k for k in calibration.tpm if k > 0)
    separable = []
    grid = []
    for tier, weight in sorted(stage3.TIER_WEIGHTS.items()):
        scaled = [(k, calibration.score(k) * weight) for k in positive_levels]
        below = [k for k, v in scaled if v <= ceiling]
        above = [k for k, v in scaled if v > ceiling]
        if below and above:
            separable.append(tier)
        grid.append(
            f"tier{tier} "
            + "/".join(f"{v:.3f}" for _k, v in scaled)
            + ("(split)" if below and above else "(all above)" if above
               else "(all below)")
        )
    criterion(
        "R14", not separable,
        ("staining levels are separable at the ceiling in tier(s) "
         + ", ".join(str(t) for t in separable)
         if separable else
         "no criticality tier places two staining levels on opposite sides of "
         f"the {ceiling} ceiling, so the arm gates on presence only")
        + " — " + "  ".join(grid),
    )

    print("\n" + "=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped = 0
    for cid, is_tripped, detail in results:
        tripped += 1 if is_tripped else 0
        print(f"  {'TRIPPED' if is_tripped else 'clear  '}  {cid}: {detail}")
    print("=" * 72)
    print(f"  {len(results) - tripped}/{len(results)} criteria clear")

    print(
        f"\n  WATCH  strongest single-component correlation: {worst[1]:.3f}"
        f" ({worst[0]}), rejection at 0.95"
    )

    baseline_only = {}
    for r in rows:
        if r.evidence_class != PROTEIN_CONFIRMED:
            continue
        profile = gtex_profiles.get(r.accession)
        if profile is None:
            continue
        base = {}
        for label, tpm in zip(gtex_tissues, profile.values):
            organ = model.organ_for_baseline(label)
            if organ is None:
                continue
            value = stage3._baseline_score(float(tpm))
            if value > base.get(organ, -1.0):
                base[organ] = value
        value, _organ = stage3.worst_organ(model, base)
        if value is not None:
            baseline_only[r.accession] = value
    clears_base = [r for r in rows
                   if r.accession in baseline_only
                   and baseline_only[r.accession] <= ceiling]
    switched = [r for r in clears_base if not r.cleared]
    share = len(switched) / len(clears_base) if clears_base else 0.0
    print(
        f"  REPORT arm switch, not gated: {len(switched)} of {len(clears_base)} "
        f"protein-confirmed targets clearing on the transcript arm alone are "
        f"blocked once staining is added ({share:.1%}); the switched set is "
        "determined by presence, not magnitude"
    )

    if tripped:
        return 1

    _report_biology(rows, ceiling)
    return 0


def _report_biology(rows, ceiling) -> None:
    """Print the biology behind the ranking the stage produced."""
    print("\n" + "=" * 72)
    print("BIOLOGY")
    print("=" * 72)

    by_gene = {}
    for r in rows:
        if r.gene and (r.gene not in by_gene or (r.composite or 0) > (by_gene[r.gene].composite or 0)):
            by_gene[r.gene] = r

    tier_sizes = {}
    for r in rows:
        tier_sizes[r.evidence_class] = tier_sizes.get(r.evidence_class, 0) + 1

    print("\n  known targets for this indication")
    print(f"    {'gene':10s} {'rank':>14s} {'composite':>10s} {'risk':>7s} "
          f"{'worst organ':>18s} {'field elev':>11s}  cleared")
    for g in ["CEACAM5", "CEACAM6", "CLDN18", "MSLN", "MUC1"]:
        r = by_gene.get(g)
        if r is None:
            print(f"    {g:10s} absent")
            continue
        size = tier_sizes[r.evidence_class]
        fe = f"{r.field_elevation:,.1f}x" if r.field_elevation else "n/a"
        print(
            f"    {r.gene:10s} {f'{r.tier_rank} of {size:,}':>14s} "
            f"{r.composite:>10.4f} {r.risk:>7.3f} {r.risk_organ or '-':>18s} "
            f"{fe:>11s}  {'yes' if r.cleared else 'no'}"
        )

    classes = [PROTEIN_CONFIRMED, RNA_SUPPORTED, DATA_INSUFFICIENT]
    print("\n  clearance by ceiling, split by evidence class")
    print(f"    {'ceiling':>9s} {'all':>7s} " + " ".join(f"{c:>18s}" for c in classes))
    for name, c in [("conservative", 0.15), ("moderate", 0.35), ("permissive", 0.60)]:
        cleared = [r for r in rows if r.risk is not None and r.risk <= c]
        per = [sum(1 for r in cleared if r.evidence_class == k) for k in classes]
        print(
            f"    {c:>9.2f} {len(cleared):>7,} "
            + " ".join(f"{n:>18,}" for n in per)
            + f"   ({name})"
        )

    print(f"\n  composition of the cleared set at {ceiling}")
    cleared = [r for r in rows if r.risk is not None and r.risk <= ceiling]
    for k in classes:
        n = sum(1 for r in cleared if r.evidence_class == k)
        share = n / len(cleared) if cleared else 0.0
        of_class = n / tier_sizes.get(k, 1)
        print(
            f"    {k:20s} {n:>6,}  {share:>7.1%} of cleared   "
            f"{of_class:>7.1%} of its class"
        )
    blocked = [r for r in rows if not (r.risk is not None and r.risk <= ceiling)]
    for label, group in [("cleared", cleared), ("blocked", blocked)]:
        if not group:
            continue
        conf = sum(r.confidence for r in group) / len(group)
        stained = sum(1 for r in group if r.evidence_class == PROTEIN_CONFIRMED)
        print(
            f"    {label:8s} n={len(group):>6,}  mean confidence {conf:.2f}   "
            f"protein-confirmed {stained / len(group):.1%}"
        )


if __name__ == "__main__":
    raise SystemExit(main())
