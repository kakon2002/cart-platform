"""What makes one cancer type different from another, in one object."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AtlasSchema:
    """How to read one single-cell submission."""

    series: str

    archive: str

    url: str

    level1_column: str
    level3_column: str

    malignant_label: str

    compartment_map: dict[str, str]

    counts_path: str = "layers/counts"

    symbol_field: str = "_index"
    ensembl_field: str = "ensg"

    treatment_column: str | None = None
    untreated_label: str | None = None

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

    tcga_project: str | None

    depmap_lineage: str | None

    gtex_bulk_label: str | None

    atlas: AtlasSchema | None

    tissue_overrides: dict[str, tuple[int, str]] = field(default_factory=dict)

    terminable_ceiling: float | None = None

    @property
    def slug(self) -> str:
        """The indication's cache-namespacing tag."""
        return self.key


def geo_url(series: str, archive: str) -> str:
    """The GEO supplementary URL for a series."""
    digits = series[3:]
    bucket = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else f"GSE{digits}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{series}/suppl/{archive}"
