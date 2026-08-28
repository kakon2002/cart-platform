"""What makes one cancer type different from another, in one object.

Before this, four module-level constants pinned the whole platform to pancreas:
a TCGA project, a GEO accession, a DepMap lineage and a GTEx column. Each lived
inside the loader that used it, so a second indication meant editing library
code rather than supplying a configuration.

**Nothing here describes the human body.** UniProt, GTEx, HPA, GENCODE and the
antibody sets are the same measurements whichever cancer is being screened, and
they stay shared. What varies is the tumour: which cohort, which atlas, which
cell lines, and which normal organ the tumour is compared against.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AtlasSchema:
    """How to read one single-cell submission.

    Every field here exists because a GEO submitter chose a name and a later
    submitter chose a different one. None of it is derivable from the accession,
    which is why it is declared rather than inferred: guessing a column name and
    getting silence is how an atlas gets read with the wrong annotations.
    """

    series: str
    #: The submitter's own filename. Not derivable from the accession.
    archive: str
    #: Full URL. The GEO bucket is derived rather than hardcoded, because the
    #: previous literal "GSE202nnn" produced a wrong URL for any other series.
    url: str
    #: The coarse annotation, which must carry the malignant call as one of its
    #: categories, and the fine annotation nested strictly inside it.
    level1_column: str
    level3_column: str
    #: The category naming malignant cells inside `level1_column`.
    malignant_label: str
    #: Coarse label -> compartment. Anything unmapped becomes "other", which is
    #: reported rather than silently dropped.
    compartment_map: dict[str, str]
    #: Raw counts layer, and the var column carrying Ensembl identifiers.
    counts_layer: str = "counts"
    ensembl_column: str = "ensg"
    #: Optional. An atlas without a treatment split still works; the untreated
    #: subset is simply absent rather than fabricated.
    treatment_column: str | None = None
    untreated_label: str | None = None
    #: Optional per-cell donor identifier, used for patient prevalence.
    patient_column: str | None = None

    @property
    def slug(self) -> str:
        """What namespaces this atlas's cache entries."""
        return self.series


@dataclass(frozen=True)
class Indication:
    """One cancer type, and everything that differs because of it."""

    key: str
    cancer_type: str
    #: GDC project. None means no bulk cohort is connected for this indication.
    tcga_project: str | None
    #: DepMap OncotreeLineage. None, or a lineage with no screened lines, means
    #: the dependency component is unavailable rather than zero.
    depmap_lineage: str | None
    #: The GTEx column the tumour-versus-normal margin is measured against.
    gtex_bulk_label: str | None
    #: None means no single-cell atlas is connected. That is not a missing 0.45
    #: of weight, it is the loss of the only component that rejects stromal and
    #: immune genes, and it makes the screen NOT_USABLE. See specs/multi-indication.md.
    atlas: AtlasSchema | None
    #: Declared, never derived. Empty where the platform default already holds.
    tissue_overrides: dict[str, tuple[int, str]] = field(default_factory=dict)
    #: Policy input for the adaptor row; None disables it for this indication.
    terminable_ceiling: float | None = None

    @property
    def slug(self) -> str:
        return self.key


def geo_url(series: str, archive: str) -> str:
    """The GEO supplementary URL for a series.

    The bucket is derived. The previous code carried "GSE202nnn" as a literal,
    which silently produced a wrong URL for every accession outside that block.
    """
    digits = series[3:]
    bucket = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else f"GSE{digits}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{series}/suppl/{archive}"
