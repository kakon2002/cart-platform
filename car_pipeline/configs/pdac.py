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
