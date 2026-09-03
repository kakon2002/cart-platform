"""Resolves dataset status against the cache rather than the network."""

from __future__ import annotations

from car_pipeline.data.antibodies import (
    AntibodySource, RELEASE_PIN as ANTIBODY_PIN)
from car_pipeline.data.depmap import DepMapSource, RELEASE_PIN as DEPMAP_PIN
from car_pipeline.data.genespan import GeneSpanSource, RELEASE_PIN as SPAN_PIN
from car_pipeline.data.gtex import GTExSource, RELEASE_PIN as GTEX_PIN
from car_pipeline.data.hpa import HPASource, RELEASE_PIN as HPA_PIN
from car_pipeline.data.singlecell import SingleCellSource
from car_pipeline.data.source import DataSource
from car_pipeline.data.tcga import TCGASource
from car_pipeline.data.trials import TrialSource, RELEASE_PIN as TRIALS_PIN
from car_pipeline.data.uniprot import UniProtSource, RELEASE_PIN as UNIPROT_PIN
from car_pipeline.schemas.spec import DatasetStatus, RequiredDataset
from car_pipeline.stages.stage1 import KNOWN_DATASET_NAMES


CONNECTORS: dict[str, type[DataSource]] = {
    UniProtSource.name: UniProtSource,
    TCGASource.name: TCGASource,
    DepMapSource.name: DepMapSource,
    HPASource.name: HPASource,
    SingleCellSource.name: SingleCellSource,
    GTExSource.name: GTExSource,
}


_ORPHANED = set(CONNECTORS) - KNOWN_DATASET_NAMES
if _ORPHANED:
    raise RuntimeError(
        "connectors registered under names no stage emits: "
        + ", ".join(sorted(_ORPHANED))
    )


def resolve_status(name: str) -> DatasetStatus:
    """Whether a named dataset is connected, readable, or neither."""
    connector = CONNECTORS.get(name)
    if connector is None:
        return DatasetStatus.NOT_CONFIGURED
    try:
        source = connector()
    except Exception:
        return DatasetStatus.UNREACHABLE

    return (
        DatasetStatus.AVAILABLE if source.is_cached() else DatasetStatus.UNREACHABLE
    )


def resolve_dataset_statuses(
    datasets: list[RequiredDataset],
) -> list[RequiredDataset]:
    """The same datasets with each one's status filled in."""
    return [d.model_copy(update={"status": resolve_status(d.name)}) for d in datasets]


if __name__ == "__main__":
    from car_pipeline.configs.pdac import PDAC_PROJECT
    from car_pipeline.stages.stage1 import build_spec

    spec = build_spec(PDAC_PROJECT)

    print("dataset status")
    for d in spec.required_datasets:
        flag = "blocking" if d.required else "optional"
        print(f"  {d.name:24s} {d.status.value:16s} {flag}")

    blocking = [d for d in spec.required_datasets if d.required]
    available = [d for d in blocking if d.status is DatasetStatus.AVAILABLE]
    score = spec.data_availability_score
    print(
        f"\n  availability score: {score:.3f}  expected 0.750"
        f"   ({len(available)} of {len(blocking)} blocking)"
    )
    missing = [d.name for d in blocking if d.status is not DatasetStatus.AVAILABLE]
    print(f"  outstanding: {', '.join(missing) or 'none'}")


def release_pins(indication) -> list[dict]:
    """Every connected source with the release this run was pinned to."""
    out = [
        {"source": UniProtSource.name, "role": "reviewed human proteome",
         "release": UNIPROT_PIN, "per_indication": False},
        {"source": HPASource.name, "role": "normal-tissue and pathology atlas",
         "release": HPA_PIN, "per_indication": False},
        {"source": GTExSource.name, "role": "bulk normal-tissue baseline",
         "release": GTEX_PIN, "per_indication": False},
        {"source": GeneSpanSource.name, "role": "gene span annotation",
         "release": SPAN_PIN, "per_indication": False},
        {"source": AntibodySource.name,
         "role": "antibody structures and named therapeutics",
         "release": ANTIBODY_PIN, "per_indication": False},
        {"source": TrialSource.name, "role": "trial registry",
         "release": TRIALS_PIN, "per_indication": False},
        {"source": TCGASource.name, "role": "tumour cohort",
         "release": indication.tcga_project, "per_indication": True},
        {"source": DepMapSource.name, "role": "dependency lineage",
         "release": f"{DEPMAP_PIN} / {indication.depmap_lineage}",
         "per_indication": True},
        {"source": SingleCellSource.name, "role": "single-cell tumour atlas",
         "release": indication.atlas.series if indication.atlas else None,
         "per_indication": True},
    ]
    return out
