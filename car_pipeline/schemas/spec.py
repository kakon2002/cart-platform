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
