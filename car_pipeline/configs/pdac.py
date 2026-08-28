"""Pancreatic ductal adenocarcinoma project definition.

``target_antigen`` is deliberately absent. That null is what selects discovery
mode, and filling it in would turn the screen into a validation of something
already assumed.
"""

from car_pipeline.schemas.project import (
    CARFormat,
    MalignancyType,
    ManufacturingConstraints,
    ProductType,
    ProjectInput,
    SafetyTolerance,
    TissueCriticalityOverride,
)

PANCREAS_RATIONALE = (
    "The tumour arises in this organ and the intended population is surgically "
    "resected, so expression in normal pancreas is a weaker objection here than "
    "the same expression in lung or heart would be. This relaxation is specific "
    "to this indication; the platform default holds pancreas at tier 1 for every "
    "other indication."
)

PDAC_PROJECT = ProjectInput(
    cancer_type="Pancreatic Ductal Adenocarcinoma",
    malignancy_type=MalignancyType.SOLID,
    product_type=ProductType.AUTOLOGOUS,
    car_format=CARFormat.AUTO,
    safety_tolerance=SafetyTolerance.CONSERVATIVE,
    # A policy input, fixed before the run and pinned in
    # specs/stage4a-architecture-routing.md §3. It is the risk this project
    # accepts from an exposure that can be stopped — an adaptor design, where
    # activation needs a separately dosed protein — as distinct from the 0.15 it
    # accepts from a T cell that cannot be withdrawn.
    #
    # This pipeline cannot measure it. Criterion A9 therefore reports the
    # admitted count across the whole sweep so the choice is visible, and A10
    # trips if this value ever stops matching the one recorded in the spec.
    terminable_risk_ceiling=0.35,
    manufacturing=ManufacturingConstraints(
        vector_payload_limit_kb=4.7,
        max_genetic_edits=2,
    ),
    tissue_criticality_overrides={
        "pancreas": TissueCriticalityOverride(
            tier=2,
            rationale=PANCREAS_RATIONALE,
        )
    },
)
