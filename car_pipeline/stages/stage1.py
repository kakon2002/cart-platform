"""Stage 1 — turn a project definition into a resolved specification."""

from __future__ import annotations

import copy
import importlib.util
import re
import uuid
from datetime import datetime, timezone

from car_pipeline.schemas.project import (
    CARFormat,
    DiscoveryMode,
    ProjectInput,
    SafetyTolerance,
)
from car_pipeline.schemas.spec import (
    DesignConstraints,
    ProjectSpec,
    RequiredDataset,
)


BACKBONE_OVERHEAD_KB = 1.2


TOLERANCE_RULES: dict[SafetyTolerance, tuple[bool, float]] = {
    SafetyTolerance.CONSERVATIVE: (True, 0.15),
    SafetyTolerance.MODERATE: (False, 0.35),
    SafetyTolerance.PERMISSIVE: (False, 0.60),
}


_SCREENING_DATASETS: list[tuple[str, str, list[int]]] = [
    ("UniProt", "membrane topology", [3]),
    ("GDC TCGA", "bulk tumour expression", [3]),
    ("DepMap", "cancer dependency", [3]),
]


_SHARED_DATASETS: list[tuple[str, str, list[int]]] = [
    ("Human Protein Atlas", "surface localisation, normal tissue", [3, 9]),
    ("Single-cell tumour atlas", "per cell type expression", [3, 9]),
    ("GTEx", "normal tissue baseline", [3, 9]),
]


_DOWNSTREAM_DATASETS: list[tuple[str, str, list[int], bool]] = [
    ("PDB", "binder structures", [5], True),
    ("SAbDab", "antibody chain and numbering annotation", [5], True),
    ("IEDB", "immunogenicity", [9], False),
    ("ClinicalTrials.gov", "prior outcomes", [9], False),
]

_RESOLVER_MODULE = "car_pipeline.data.availability"


KNOWN_DATASET_NAMES = frozenset(
    [n for n, _, _ in _SCREENING_DATASETS]
    + [n for n, _, _ in _SHARED_DATASETS]
    + [n for n, _, _, _ in _DOWNSTREAM_DATASETS]
)


def _slug(text: str) -> str:
    """Lowercase the text into a hyphenated identifier."""
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def _project_id(cancer_type: str, created_at: datetime) -> str:
    """Slug, timestamp and a random suffix."""
    stamp = created_at.strftime("%Y%m%dT%H%M%S")
    return f"{_slug(cancer_type)}-{stamp}-{uuid.uuid4().hex[:8]}"


def _datasets(mode: DiscoveryMode) -> list[RequiredDataset]:
    """The datasets this mode needs, screening sources included only in discovery."""
    rows: list[RequiredDataset] = []
    if mode is DiscoveryMode.DISCOVER:
        for name, purpose, stages in _SCREENING_DATASETS:
            rows.append(RequiredDataset(name=name, purpose=purpose, stages=stages))
    for name, purpose, stages in _SHARED_DATASETS:
        rows.append(RequiredDataset(name=name, purpose=purpose, stages=stages))
    for name, purpose, stages, required in _DOWNSTREAM_DATASETS:
        rows.append(
            RequiredDataset(
                name=name, purpose=purpose, stages=stages, required=required
            )
        )
    return rows


def _design_constraints(inputs: ProjectInput) -> DesignConstraints:
    """Turn the declared tolerance and manufacturing limits into design constraints."""
    require_switch, ceiling = TOLERANCE_RULES[inputs.safety_tolerance]
    budget = inputs.manufacturing.vector_payload_limit_kb - BACKBONE_OVERHEAD_KB
    return DesignConstraints(
        max_construct_kb=round(budget, 6),
        max_genetic_edits=inputs.manufacturing.max_genetic_edits,
        require_safety_switch=require_switch,
        allowed_car_formats=[f for f in CARFormat if f is not CARFormat.AUTO],
        normal_tissue_risk_ceiling=ceiling,
        terminable_risk_ceiling=inputs.terminable_risk_ceiling,
    )


def build_spec(inputs: ProjectInput, resolve_sources: bool = True) -> ProjectSpec:
    """Resolve a project definition into a specification."""
    local = copy.deepcopy(inputs)
    created_at = datetime.now(timezone.utc)
    mode = local.discovery_mode
    datasets = _datasets(mode)

    if resolve_sources:
        datasets = _resolve(datasets)

    return ProjectSpec(
        project_id=_project_id(local.cancer_type, created_at),
        created_at=created_at,
        inputs=local,
        discovery_mode=mode,
        required_datasets=datasets,
        design_constraints=_design_constraints(local),
    )


def _resolve(datasets: list[RequiredDataset]) -> list[RequiredDataset]:
    """Ask the cache what is present, if a resolver has been built yet."""
    if importlib.util.find_spec(_RESOLVER_MODULE) is None:
        return datasets
    from car_pipeline.data.availability import resolve_dataset_statuses

    return resolve_dataset_statuses(datasets)
