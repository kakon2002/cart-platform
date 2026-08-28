"""Resolved project specification produced by stage 1."""

from __future__ import annotations

from datetime import datetime
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

from car_pipeline.schemas.project import CARFormat, DiscoveryMode, ProjectInput


class DatasetStatus(str, Enum):
    """Three distinct states, deliberately not two.

    ``not_configured`` means no connector exists for this source at all.
    ``unreachable`` means the connector exists but the data is not readable.
    Collapsing them hides whether the gap is unbuilt or merely unfetched.
    """

    AVAILABLE = "available"
    UNREACHABLE = "unreachable"
    NOT_CONFIGURED = "not_configured"


class RequiredDataset(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    purpose: str
    stages: list[int]
    required: bool = True
    status: DatasetStatus = DatasetStatus.NOT_CONFIGURED


class DesignConstraints(BaseModel):
    model_config = ConfigDict(extra="forbid")

    max_construct_kb: float
    max_genetic_edits: int
    require_safety_switch: bool
    allowed_car_formats: list[CARFormat]
    normal_tissue_risk_ceiling: float = Field(ge=0.0, le=1.0)
    #: The ceiling for an exposure that can be stopped. An adaptor design does
    #: not make an antigen safer - the adaptor still binds it - it makes the
    #: exposure terminable, because activation needs a separately dosed protein.
    #: Magnitude and reversibility are different axes, so this is a second
    #: number rather than an adjustment to the first.
    #:
    #: Optional and defaulted to None on purpose. Both ceilings are policy
    #: inputs, not measurements; a default here would be this code quietly
    #: setting clinical policy, so its absence disables the adaptor row instead.
    terminable_risk_ceiling: float | None = Field(default=None, ge=0.0, le=1.0)


class ProjectSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    created_at: datetime
    inputs: ProjectInput
    discovery_mode: DiscoveryMode
    required_datasets: list[RequiredDataset]
    design_constraints: DesignConstraints

    @property
    def data_availability_score(self) -> float:
        """Fraction of blocking datasets that are actually readable.

        Counts only ``required=True`` entries. Including optional ones caps an
        otherwise complete run below 1.0 and makes a real blocking gap look
        milder than it is.
        """
        blocking = [d for d in self.required_datasets if d.required]
        if not blocking:
            return 0.0
        available = sum(1 for d in blocking if d.status is DatasetStatus.AVAILABLE)
        return available / len(blocking)
