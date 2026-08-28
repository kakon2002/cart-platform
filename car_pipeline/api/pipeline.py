"""The pipeline behind the API, run once per job."""

from __future__ import annotations

from car_pipeline.configs.registry import resolve as resolve_indication
from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.genespan import GeneSpanSource
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import PRIMARY_TUMOUR, TCGASource, match_surface as tcga_match
from car_pipeline.data.trials import TrialSource
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import (
    routing, stage3, stage4, stage5, stage6, stage9, stage10, stage11,
)
from car_pipeline.stages.stage1 import build_spec

import numpy as np

STAGES = ("sources", "screen", "pairing", "binders", "constructs",
          "safety", "developability", "ranking")


USABLE = "USABLE"
NOT_USABLE = "NOT_USABLE"


def project_for(cancer_type: str):
    """The configured project for an indication, or a refusal naming what exists."""
    _indication, project = resolve_indication(cancer_type)
    return project


def run(cancer_type: str, progress=lambda stage, note="": None) -> dict:
    """Run the whole pipeline and return everything the API serves."""
    progress("sources", "loading cached sources")
    indication, project = resolve_indication(cancer_type)
    spec = build_spec(project)
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    overrides = {o: ov.tier
                 for o, ov in spec.inputs.tissue_criticality_overrides.items()}

    surface, _ = load_surface()
    surface_by_acc = {r.accession: r for r in surface}
    atlas = HPASource().load()
    by_acc, by_sym = atlas_index(atlas)
    gtex_profiles, gtex_tissues, _ = GTExSource().match_surface(surface, by_acc)

    unavailable: list[str] = []

    missing_required = [
        name for name, value in (
            ("tumour cohort", indication.tcga_project),
            ("normal-tissue denominator", indication.gtex_bulk_label),
        ) if not value
    ]
    if missing_required:
        raise ValueError(
            f"{indication.cancer_type} declares no "
            + " and no ".join(missing_required)
            + ". These have no default: falling back would screen this "
              "indication against another one's data and report it under this "
              "name."
        )
    cohort = TCGASource(indication.tcga_project).load()
    cohort_join = tcga_match(cohort, surface, by_acc)

    if indication.atlas is None:
        cells_atlas, cell_index = None, {}
        unavailable.append(
            "single-cell atlas: none is connected for this indication, so "
            "malignant_expression and malignant_vs_stroma cannot be measured"
        )
    else:
        cells_atlas = SingleCellSource(indication.atlas).load()
        cell_index = cell_match(cells_atlas, surface, by_acc)

    try:
        dependency, _ = DepMapSource(indication.depmap_lineage).load()
        dep_index = gene_index(dependency)
    except Exception as exc:
        dependency, dep_index = None, {}
        unavailable.append(
            f"dependency lineage {indication.depmap_lineage!r}: "
            f"{type(exc).__name__}: {exc}"
        )
    coverage_rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, cohort_join)

    progress("screen", "ranking the surface proteome")
    calibration = stage3.calibrate_atlas_levels(
        surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
        stage3.RiskModel(overrides=overrides))
    ranked, model, _ = stage3.rank(
        coverage_rows, surface_by_acc, by_acc, by_sym, cells_atlas, cell_index,
        gtex_profiles, gtex_tissues, cohort, cohort_join, dependency, dep_index,
        overrides, ceiling, calibration,
        margin_label=indication.gtex_bulk_label)
    s3_hash = stage3.configuration_hash(
        overrides, ceiling, calibration=calibration,
        margin_label=indication.gtex_bulk_label)

    if indication.atlas is None:
        return {
            "indication": indication,
            "usability": NOT_USABLE,
            "unavailable": unavailable,
            "spec": spec,
            "ceiling": ceiling,
            "ranked": ranked,
            "pool": [],
            "pairs": [],
            "decisions": [],
            "binders": {},
            "constructs": [],
            "gated": [],
            "developability": [],
            "developability_status": "NOT_SCORED",
            "final": [],
            "attrition": {},
            "status": NOT_USABLE,
            "stage3_hash": s3_hash,
            "stage4_hash": None,
            "reasons": [
                "No single-cell atlas is connected for this indication, so "
                "malignant_expression and malignant_vs_stroma are unmeasured.",
                "malignant_vs_stroma is the only component that rejects "
                "stromal and immune genes. Without it nothing distinguishes a "
                "tumour antigen from a gene expressed by the lymphocytes beside "
                "the tumour, and the ranking fills with immunoglobulin, T-cell "
                "receptor and MHC-II genes while still scoring ~3,400 targets.",
                "That is a structural loss, not a missing 0.45 of weight, so no "
                "ranking is served for this indication.",
            ],
        }

    progress("pairing", "evaluating pairs")
    pool = stage4.build_pool(ranked)
    per_organ = {
        r.accession: stage3.per_organ_scores(
            model, by_acc.get(r.accession) or by_sym.get(r.gene),
            gtex_profiles[r.accession].values
            if r.accession in gtex_profiles else None,
            gtex_tissues, calibration)
        for r in pool
    }

    cells = (SingleCellSource(indication.atlas)
             .load_malignant(sorted({r.gene for r in pool}))
             if indication.atlas is not None else None)
    pairs = stage4.evaluate(pool, per_organ, model, ceiling, cells)
    try:
        stage4.annotate_span_context(pairs, GeneSpanSource().load())
    except Exception:
        pass
    primary = cohort.sample_types == PRIMARY_TUMOUR
    tumour_tpm = {}
    for r in pool:
        cj = cohort_join.get(r.accession)
        if cj is not None:
            value = float(np.median(cohort.values[primary, cj[0]]))
            if np.isfinite(value):
                tumour_tpm[r.gene] = value

    tolerances = routing.Tolerances(
        persistent=ceiling,
        terminable=spec.design_constraints.terminable_risk_ceiling,
    )
    decisions_obj = stage4.decide(pool, pairs, tumour_tpm, tolerances)
    decisions = stage4.decision_rows(decisions_obj)
    s4_hash = stage4.configuration_hash(
        s3_hash, [r.gene for r in pool], tolerances)

    progress("binders", "retrieving binders")

    records = stage5.load_or_retrieve(decisions, AntibodySource(), s4_hash)
    binders = {r.gene: r for r in records}

    progress("constructs", "assembling")
    constructs = stage6.build(decisions, binders)
    by_construct = {c.gene: c for c in constructs}

    progress("safety", "gating")
    risks = {r.gene: (r.risk, r.risk_organ) for r in ranked if r.gene}
    trials = TrialSource(antigens=[d["gene"] for d in decisions]).load()
    gated = stage9.gate(decisions, binders, risks, trials, ceiling)
    by_gate = {g.gene: g for g in gated}

    progress("developability", "scoring sequences")
    dev_rows, dev_status = stage10.assess(binders)
    liabilities: dict[str, list] = {}
    for row in dev_rows:
        liabilities.setdefault(row.gene, []).append(row)

    progress("ranking", "attribution and Pareto front")
    composites = {r.gene: r.composite for r in ranked if r.gene}
    final, attrition, status = stage11.rank(
        decisions, binders, by_construct, by_gate, liabilities, composites, ceiling)

    usability = NOT_USABLE if indication.atlas is None else USABLE

    return {
        "indication": indication,
        "usability": usability,
        "unavailable": unavailable,
        "spec": spec,
        "ceiling": ceiling,
        "ranked": ranked,
        "pool": pool,
        "pairs": pairs,
        "decisions": decisions,
        "binders": binders,
        "constructs": constructs,
        "gated": gated,
        "developability": dev_rows,
        "developability_status": dev_status,
        "final": final,
        "attrition": attrition,
        "status": status,
        "stage3_hash": s3_hash,
        "stage4_hash": s4_hash,
    }


UNKNOWN_TARGET = "UNKNOWN_TARGET"
SUITABLE = "SUITABLE"
CONDITIONAL = "CONDITIONAL"
UNSUITABLE = "UNSUITABLE"


NOT_ASSESSED = "NOT_ASSESSED"


def validate(cancer_type: str, target: str, progress=lambda s, n="": None) -> dict:
    """Evaluate one supplied target, rather than screening for targets."""
    symbol = (target or "").strip().upper()
    if not symbol:
        raise ValueError("a target symbol is required in validation mode")

    screen = run(cancer_type, progress)
    ranked = [x for x in screen["ranked"] if x.gene]
    scored = sorted(
        (x for x in ranked if x.composite is not None),
        key=lambda x: (-x.composite, x.gene, x.accession),
    )
    rank_of = {x.gene: i + 1 for i, x in enumerate(scored)}

    entry = next((x for x in ranked if x.gene.upper() == symbol), None)
    base = {
        "mode": "A",
        "cancer_type": screen["indication"].cancer_type,
        "target": symbol,
        "screen": screen,
    }

    if entry is None:
        from car_pipeline.data.uniprot import UniProtSource
        every = UniProtSource().load()
        match = next((r for r in every if (r.gene or "").upper() == symbol), None)
        if match is None:
            return {
                **base,
                "verdict": UNKNOWN_TARGET,
                "reasons": [
                    f"{symbol} is not a gene symbol in the reviewed human "
                    f"proteome ({len(every):,} entries, release pinned).",
                ],
            }
        return {
            **base,
            "verdict": UNSUITABLE,
            "accession": match.accession,
            "reasons": [
                f"{symbol} ({match.accession}) is in the reviewed proteome but "
                "is not surface-accessible, so a CAR binder has nothing to "
                "engage. This is a rejection on topology, not an absence of "
                "evidence.",
                f"Membrane class: {match.membrane_class}; "
                f"outward-facing: {match.outward}; attached: {match.attached}.",
            ],
        }

    decision = next(
        (d for d in screen["decisions"] if d["gene"] == entry.gene), None)
    binder = screen["binders"].get(entry.gene)
    construct = next(
        (c for c in screen["constructs"] if c.gene == entry.gene), None)

    in_pool = decision is not None
    architecture = (decision or {}).get("architecture") or "NOT_EVALUATED"
    if not in_pool:
        verdict = NOT_ASSESSED
    elif architecture == "CONVENTIONAL":
        verdict = SUITABLE
    elif architecture in ("AND_GATE", "ADAPTOR"):
        verdict = CONDITIONAL
    else:
        verdict = UNSUITABLE

    reasons = []
    if entry.composite is None:
        reasons.append(
            "Below the evidence floor: too little of the score is measured for "
            "this target to be ranked at all.")
    else:
        reasons.append(
            f"Ranks {rank_of.get(entry.gene, '?')} of {len(scored)} on tumour "
            f"attractiveness (composite {entry.composite:.4f}, "
            f"{entry.measured_weight:.2f} of the evidence measured).")
    if entry.risk is not None:
        reasons.append(
            f"Normal-tissue risk {entry.risk:.4f}, peak organ "
            f"{entry.risk_organ}.")
    if not in_pool:
        reasons.append(
            f"Not carried into pairing: only the top {len(screen['pool'])} of "
            f"{len(scored)} ranked targets are, and this one is not among them. "
            "No architecture was routed for it, which is not the same as no "
            "architecture fitting it.")
    elif decision.get("route_reason"):
        reasons.append(decision["route_reason"])
    if binder is not None:
        n = len(binder.sequence) + len(binder.structure)
        reasons.append(
            f"{n} binder candidate(s) retrieved"
            + (f"; {len(binder.sequence)} carry a usable sequence."
               if binder.sequence else ", none carrying a usable sequence."))
    if construct is not None and construct.total_bp:
        reasons.append(
            f"Assembles at {construct.total_bp} bp against the budget "
            f"({construct.verdict}).")

    return {
        **base,
        "verdict": verdict,
        "accession": entry.accession,
        "rank": rank_of.get(entry.gene),
        "of": len(scored),
        "composite": entry.composite,
        "measured_weight": entry.measured_weight,
        "risk": entry.risk,
        "risk_organ": entry.risk_organ,
        "evidence_class": entry.evidence_class,
        "architecture": architecture,
        "outcome": (decision or {}).get("outcome"),
        "reasons": reasons,
    }
