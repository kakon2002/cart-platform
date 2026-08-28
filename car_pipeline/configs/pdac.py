"""Pancreatic ductal adenocarcinoma project definition.

``target_antigen`` is deliberately absent. That null is what selects discovery
mode, and filling it in would turn the screen into a validation of something
already assumed.
"""

from car_pipeline.configs.indication import AtlasSchema, Indication, geo_url
from car_pipeline.schemas.project import (
    CARFormat,
    MalignancyType,
    ManufacturingConstraints,
    ProductType,
    ProjectInput,
    SafetyTolerance,
    TissueCriticalityOverride,
)

#: The reference atlas, described rather than assumed. Every value here was a
#: module-level constant inside the loader until a second indication needed a
#: different one; the loader now reads them from this object.
PDAC_ATLAS = AtlasSchema(
    series="GSE202051",
    archive="GSE202051_totaldata-final-toshare.h5ad.gz",
    url=geo_url("GSE202051", "GSE202051_totaldata-final-toshare.h5ad.gz"),
    level1_column="Level 1 Annotation",
    level3_column="Level 3 Annotation",
    malignant_label="Epithelial (malignant)",
    compartment_map={
        "Epithelial (malignant)": "malignant",
        "Cancer-associated fibroblast": "fibroblast",
        "Epithelial (non-malignant)": "epithelial non-malignant",
        "Lymphoid": "immune",
        "Myeloid": "immune",
        "Endothelial": "endothelial",
    },
    treatment_column="treatment_status",
    untreated_label="Untreated",
    patient_column="pid",
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


#: The reference indication, as an Indication object. Everything that was a
#: module constant inside a loader now lives here.
PDAC = Indication(
    key="pdac",
    cancer_type="Pancreatic Ductal Adenocarcinoma",
    tcga_project="TCGA-PAAD",
    depmap_lineage="Pancreas",
    # Four GTEx pancreas columns exist; three are cell-sorted fractions. The
    # bulk one is the denominator, which was a judgement call recorded in
    # stage 3 rather than an obvious choice.
    gtex_bulk_label="Pancreas",
    atlas=PDAC_ATLAS,
    tissue_overrides={"pancreas": (2, PANCREAS_RATIONALE)},
    terminable_ceiling=0.35,
)
