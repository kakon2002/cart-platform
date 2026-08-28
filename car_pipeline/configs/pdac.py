"""Pancreatic ductal adenocarcinoma project definition."""

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


PDAC = Indication(
    key="pdac",
    cancer_type="Pancreatic Ductal Adenocarcinoma",
    tcga_project="TCGA-PAAD",
    depmap_lineage="Pancreas",
    gtex_bulk_label="Pancreas",
    atlas=PDAC_ATLAS,
    tissue_overrides={"pancreas": (2, PANCREAS_RATIONALE)},
    terminable_ceiling=0.35,
)
