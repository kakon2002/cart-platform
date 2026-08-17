"""Project definition input schema.

Every model forbids unknown fields. A mistyped field name has to fail loudly:
silently accepting ``target_antigens`` would leave a run in discovery mode while
appearing to carry a supplied target, and nothing downstream could tell the
difference.
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, field_validator


class MalignancyType(str, Enum):
    SOLID = "solid"
    HEMATOLOGICAL = "hematological"


class ProductType(str, Enum):
    AUTOLOGOUS = "autologous"
    ALLOGENEIC = "allogeneic"


class CARFormat(str, Enum):
    AUTO = "auto"
    CONVENTIONAL = "conventional"
    DUAL_TARGET = "dual_target"
    LOGIC_GATED = "logic_gated"
    SWITCHABLE = "switchable"
    ARMORED = "armored"


class SafetyTolerance(str, Enum):
    CONSERVATIVE = "conservative"
    MODERATE = "moderate"
    PERMISSIVE = "permissive"


class BinderFormat(str, Enum):
    SCFV = "scFv"
    VH_VL = "VH_VL"
    VHH = "VHH"
    LIGAND = "ligand"


class DiscoveryMode(str, Enum):
    """A validates a supplied antigen; B screens for one."""

    VALIDATE = "A"
    DISCOVER = "B"


class ExistingBinder(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str
    format: BinderFormat


class ManufacturingConstraints(BaseModel):
    model_config = ConfigDict(extra="forbid")

    vector_payload_limit_kb: float = 4.7
    max_genetic_edits: int = 2


class TissueCriticalityOverride(BaseModel):
    """Relaxes or tightens the platform default criticality tier for one tissue.

    The rationale is required and travels into the output header, so a reader can
    see which safety defaults were moved and why.
    """

    model_config = ConfigDict(extra="forbid")

    tier: int = Field(ge=1, le=3)
    rationale: str = Field(min_length=20)


class ProjectInput(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cancer_type: str
    malignancy_type: MalignancyType

    cancer_subtype: str | None = None
    target_antigen: str | None = None
    patient_subgroup: str | None = None

    product_type: ProductType = ProductType.AUTOLOGOUS
    car_format: CARFormat = CARFormat.AUTO
    safety_tolerance: SafetyTolerance = SafetyTolerance.CONSERVATIVE

    manufacturing: ManufacturingConstraints = Field(
        default_factory=ManufacturingConstraints
    )
    existing_binder: ExistingBinder | None = None
    tissue_criticality_overrides: dict[str, TissueCriticalityOverride] = Field(
        default_factory=dict
    )

    @field_validator("cancer_type")
    @classmethod
    def _require_cancer_type(cls, value: str) -> str:
        trimmed = value.strip()
        if not trimmed:
            raise ValueError("cancer_type must not be blank")
        return trimmed

    @field_validator("target_antigen")
    @classmethod
    def _blank_antigen_is_absent(cls, value: str | None) -> str | None:
        """An empty string from a form is an absent target, not a supplied one.

        Left as-is it would count as a target and skip the entire screen.
        """
        if value is None:
            return None
        trimmed = value.strip()
        return trimmed or None

    @property
    def discovery_mode(self) -> DiscoveryMode:
        """Computed, never stored. A field would let the two disagree."""
        if self.target_antigen is None:
            return DiscoveryMode.DISCOVER
        return DiscoveryMode.VALIDATE
