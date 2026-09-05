"""Every indication the platform is configured for."""

from __future__ import annotations

from car_pipeline.configs.breast import BREAST, BREAST_PROJECT
from car_pipeline.configs.indication import Indication
from car_pipeline.configs.pdac import PDAC, PDAC_PROJECT


INDICATIONS: dict[str, Indication] = {
    PDAC.cancer_type.lower(): PDAC,
    BREAST.cancer_type.lower(): BREAST,
}

PROJECTS = {
    PDAC.cancer_type.lower(): PDAC_PROJECT,
    BREAST.cancer_type.lower(): BREAST_PROJECT,
}


ALIASES = {
    "pdac": PDAC.cancer_type.lower(),
    "pancreatic": PDAC.cancer_type.lower(),
    "pancreatic cancer": PDAC.cancer_type.lower(),
    "brca": BREAST.cancer_type.lower(),
    "breast": BREAST.cancer_type.lower(),
    "breast cancer": BREAST.cancer_type.lower(),
    "breast carcinoma": BREAST.cancer_type.lower(),
}
# The bare organ names are aliases because they are what a reader types. They
# are unambiguous only while one indication per organ is configured; a second
# breast or pancreatic indication makes "breast" a choice between two, and the
# alias must then be removed rather than silently pointing at whichever was
# registered first.


def resolve(cancer_type: str) -> tuple[Indication, object]:
    """The indication and project spec for a cancer type, or a named refusal."""
    key = (cancer_type or "").strip().lower()
    key = ALIASES.get(key, key)
    if key not in INDICATIONS:
        raise ValueError(
            f"no configuration for {cancer_type!r}. Configured indications are: "
            + "; ".join(sorted(i.cancer_type for i in INDICATIONS.values()))
            + ". An indication needs a tumour cohort, a single-cell atlas, a "
            "dependency lineage and a normal-tissue denominator declared before "
            "it can be screened; none of those is derivable from the name."
        )
    return INDICATIONS[key], PROJECTS[key]


def registered() -> list[dict]:
    """What each configured indication has, for the availability endpoint."""
    return [
        {
            "cancer_type": ind.cancer_type,
            "key": ind.key,
            "cohort": ind.tcga_project,
            "atlas": ind.atlas.series if ind.atlas else None,
            "dependency_lineage": ind.depmap_lineage,
            "normal_denominator": ind.gtex_bulk_label,
        }
        for ind in sorted(INDICATIONS.values(), key=lambda i: i.cancer_type)
    ]
