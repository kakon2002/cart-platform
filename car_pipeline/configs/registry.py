"""Every indication the platform is configured for.

The platform's *design* is cancer-agnostic; its *configuration* is not, and that
distinction is the honest one. An indication needs a cohort, an atlas, a lineage
and a normal-tissue denominator to be named, and none of those can be derived
from a cancer type string. So a request for something unregistered is refused by
name rather than answered with another indication's results.
"""

from __future__ import annotations

from car_pipeline.configs.breast import BREAST, BREAST_PROJECT
from car_pipeline.configs.indication import Indication
from car_pipeline.configs.pdac import PDAC, PDAC_PROJECT

#: Keyed by the lowercased cancer type as a caller would send it.
INDICATIONS: dict[str, Indication] = {
    PDAC.cancer_type.lower(): PDAC,
    BREAST.cancer_type.lower(): BREAST,
}

PROJECTS = {
    PDAC.cancer_type.lower(): PDAC_PROJECT,
    BREAST.cancer_type.lower(): BREAST_PROJECT,
}

#: Short aliases, so a caller need not reproduce the full oncological name.
ALIASES = {
    "pdac": PDAC.cancer_type.lower(),
    "pancreatic cancer": PDAC.cancer_type.lower(),
    "brca": BREAST.cancer_type.lower(),
    "breast cancer": BREAST.cancer_type.lower(),
    "breast carcinoma": BREAST.cancer_type.lower(),
}


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
