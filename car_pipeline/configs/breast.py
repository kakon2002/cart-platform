"""Invasive breast carcinoma — the second indication.

Chosen over colorectal and lung because its atlas is a structural drop-in for the
loader rather than a rewrite: the matrix is CSR in the same layout, the
normalisation matches the reference submission (per-cell sum of expm1 ~9,974
against ~10,000), the variable index is 100% Ensembl, and malignancy sits at the
*coarse* annotation level, which is the shape the malignant-cell reader requires.

Two things it brings that the reference indication does not:

* **113 solid-tissue normal samples against pancreas's 4.** The tumour-versus-
  normal margin is computed from that arm, so this indication supports it 28x
  better than the one the platform was designed on.
* **No criticality override.** Breast defaults to tier 3, so nothing has to be
  declared. Pancreas needed a tier-2 override because it defaults to tier 1, and
  that override was the single hardest judgement in its config.

``target_antigen`` is absent here as it is for every indication. The null selects
discovery mode; Mode A supplies a target explicitly and does not touch this file.
"""

from car_pipeline.configs.indication import AtlasSchema, Indication
from car_pipeline.schemas.project import (
    CARFormat,
    MalignancyType,
    ManufacturingConstraints,
    ProductType,
    ProjectInput,
    SafetyTolerance,
)

#: Wu et al. 2021, GEO GSE176078, taken as the CELLxGENE standardised export
#: rather than the GEO supplementary files: the export is a single h5ad with the
#: annotations already in `obs`, where the GEO deposit is a matrix and metadata
#: that would have to be reassembled. 843,892,052 bytes, verified reachable with
#: HDF5 magic at byte 0.
#:
#: The URL is pinned to a specific asset. CELLxGENE re-issues asset UUIDs when a
#: dataset version changes, so this is a version pin, and a fetch that 404s is a
#: signal the dataset moved rather than a transient failure.
BREAST_ATLAS = AtlasSchema(
    series="GSE176078",
    archive="GSE176078_breast_cellxgene.h5ad",
    url=("https://datasets.cellxgene.cziscience.com/"
         "22a27631-aecf-463b-86c6-a8334a2f2cf2.h5ad"),
    # Malignancy is a category of the coarse level, which is what the malignant
    # reader requires. Corroborated independently: cross-tabulating
    # celltype_major against normal_cell_call gives (Cancer Epithelial, cancer)
    # = 24,489 exactly and (Normal Epithelial, normal) = 4,355.
    level1_column="celltype_major",
    level3_column="celltype_minor",
    malignant_label="Cancer Epithelial",
    # The submitter's nine coarse branches. Two immune branches collapse into
    # one compartment, exactly as the reference atlas does with Lymphoid and
    # Myeloid. PVL (perivascular-like) maps to fibroblast: it is a mural
    # stromal population, and the alternative is "other", which would exclude
    # 5,423 stromal cells from the comparison a tumour antigen has to win.
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
    treatment_column="treatment_status",
    untreated_label="Naive",
    patient_column="donor_id",
)

BREAST = Indication(
    key="brca",
    cancer_type="Invasive Breast Carcinoma",
    tcga_project="TCGA-BRCA",
    depmap_lineage="Breast",
    # GTEx v10 declares exactly one breast column. Pancreas has four, three of
    # them cell-sorted fractions, which is why the reference indication's
    # denominator was itself a judgement call.
    gtex_bulk_label="Breast_Mammary_Tissue",
    atlas=BREAST_ATLAS,
    # Deliberately empty. Breast is tier 3 by platform default; nothing about
    # this indication argues for relaxing it, and inventing an override to make
    # more targets clear is the move this project exists to refuse.
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
