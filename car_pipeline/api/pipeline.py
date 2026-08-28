"""The pipeline behind the API, run once per job.

Kept apart from the HTTP layer so the run can be exercised without a server.

**A run that produces no buildable design is a completed run**, not a failure.
Every stage returns what it measured and the job records that; nothing here maps
an empty result onto an error, because the emptiness is the finding.
"""

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

#: An indication with no single-cell atlas can still be scored, and the result
#: must not be presented as a ranking. Measured: dropping the atlas costs C1 and
#: C2, leaves 3,399 of 3,466 targets scored, and fills the top of the pool with
#: immunoglobulin, TCR and MHC-II genes -- because C2 is the ONLY component that
#: rejects stromal and immune expression. Losing 0.45 of weight is survivable
#: arithmetic; losing the only discriminator against stroma is not, and the
#: renormalised mean hides that by rescaling what remains as the whole score.
USABLE = "USABLE"
NOT_USABLE = "NOT_USABLE"


def project_for(cancer_type: str):
    """The configured project for an indication, or a refusal naming what exists."""
    _indication, project = resolve_indication(cancer_type)
    return project


def run(cancer_type: str, progress=lambda stage, note="": None) -> dict:
    """Run the whole pipeline and return everything the API serves.

    `progress` is called with a stage name so a polling client can see movement;
    the screen takes minutes and a request/response shape would time out.
    """
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
    # Every tumour-side source is resolved from the indication. Nothing below
    # reaches for a module constant naming one cohort or one atlas.
    unavailable: list[str] = []
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
    except Exception as exc:                       # noqa: BLE001
        # A lineage with no screened cell lines used to reach np.vstack and
        # raise "need at least one array to concatenate". Six of the 36 lineages
        # in the cached model table qualify. Escape resistance is 0.05 of the
        # weight and its absence changes almost nothing, so this degrades and
        # names the source rather than ending the run.
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
    cells = SingleCellSource().load_malignant(sorted({r.gene for r in pool}))
    pairs = stage4.evaluate(pool, per_organ, model, ceiling, cells)
    try:
        stage4.annotate_span_context(pairs, GeneSpanSource().load())
    except Exception:
        pass                      # annotation only; never loses the run
    primary = cohort.sample_types == PRIMARY_TUMOUR
    tumour_tpm = {}
    for r in pool:
        cj = cohort_join.get(r.accession)
        if cj is not None:
            value = float(np.median(cohort.values[primary, cj[0]]))
            if np.isfinite(value):
                tumour_tpm[r.gene] = value
    # Stage 4a. Both ceilings come from the project spec; neither is invented
    # here. A project with no declared terminable ceiling gets no adaptor row.
    tolerances = routing.Tolerances(
        persistent=ceiling,
        terminable=spec.design_constraints.terminable_risk_ceiling,
    )
    decisions_obj = stage4.decide(pool, pairs, tumour_tpm, tolerances)
    decisions = stage4.decision_rows(decisions_obj)
    s4_hash = stage4.configuration_hash(
        s3_hash, [r.gene for r in pool], tolerances)

    progress("binders", "retrieving binders")
    # The Stage 4 hash, which is the slot this artifact records and what every
    # downstream hash chains from. Passing the Stage 3 hash here would corrupt
    # the provenance of a cache three other drivers read.
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

    # An atlas-less run is refused as a ranking rather than served with a
    # caveat. A number that looks like an answer is worse than a refusal,
    # because only one of the two gets checked.
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
