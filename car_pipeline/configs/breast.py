"""Invasive breast carcinoma — the second indication."""

from car_pipeline.configs.indication import AtlasSchema, Indication
from car_pipeline.schemas.project import (
    CARFormat,
    MalignancyType,
    ManufacturingConstraints,
    ProductType,
    ProjectInput,
    SafetyTolerance,
)


BREAST_ATLAS = AtlasSchema(
    series="GSE176078",
    archive="GSE176078_breast_cellxgene.h5ad",
    url=("https://datasets.cellxgene.cziscience.com/"
         "22a27631-aecf-463b-86c6-a8334a2f2cf2.h5ad"),
    level1_column="celltype_major",
    level3_column="celltype_minor",
    malignant_label="Cancer Epithelial",
    compartment_map={
        "Cancer Epithelial": "malignant",
        "CAFs": "fibroblast",
        "PVL": "fibroblast",
        "Normal Epithelial": "epithelial non-malignant",
        "T-cells": "immune",
        "B-cells": "immune",
        "Myeloid": "immune",
        "Plasmablasts": "immune",
        "Endothelial": "endothelial",
    },
    counts_path="raw/X",
    symbol_field="feature_name",
    ensembl_field="_index",
    treatment_column="treatment_status",
    untreated_label="Naïve",
    patient_column="donor_id",
)

BREAST = Indication(
    key="brca",
    cancer_type="Invasive Breast Carcinoma",
    tcga_project="TCGA-BRCA",
    depmap_lineage="Breast",
    gtex_bulk_label="Breast_Mammary_Tissue",
    atlas=BREAST_ATLAS,
    tissue_overrides={},
    terminable_ceiling=0.35,
)

BREAST_PROJECT = ProjectInput(
    cancer_type=BREAST.cancer_type,
    malignancy_type=MalignancyType.SOLID,
    product_type=ProductType.AUTOLOGOUS,
    car_format=CARFormat.AUTO,
    safety_tolerance=SafetyTolerance.CONSERVATIVE,
    manufacturing=ManufacturingConstraints(
        vector_payload_limit_kb=4.7,
        max_genetic_edits=2,
    ),
    tissue_criticality_overrides={},
    terminable_risk_ceiling=BREAST.terminable_ceiling,
)
