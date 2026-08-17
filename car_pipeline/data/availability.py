"""Resolves dataset status against the cache rather than the network.

The honest question for a run about to start is whether the data is on this
machine, not whether a server is answering. A source that responds but has not
been fetched cannot feed a stage, and one that has been fetched works whether or
not anything is reachable right now.

Only manifests are read. The deep integrity check — row counts, gene axes,
checksums — belongs at load, where a failure can be attributed to a specific
file rather than to the whole source.
"""

from __future__ import annotations

from car_pipeline.data.depmap import DepMapSource
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource
from car_pipeline.data.singlecell import SingleCellSource
from car_pipeline.data.source import DataSource
from car_pipeline.data.tcga import TCGASource
from car_pipeline.data.uniprot import UniProtSource
from car_pipeline.schemas.spec import DatasetStatus, RequiredDataset

# Datasets with a connector. Anything absent from this map has none, which is a
# different finding from a connector whose data will not read.
CONNECTORS: dict[str, type[DataSource]] = {
    UniProtSource.name: UniProtSource,
    TCGASource.name: TCGASource,
    DepMapSource.name: DepMapSource,
    HPASource.name: HPASource,
    SingleCellSource.name: SingleCellSource,
    GTExSource.name: GTExSource,
}


def resolve_status(name: str) -> DatasetStatus:
    connector = CONNECTORS.get(name)
    if connector is None:
        return DatasetStatus.NOT_CONFIGURED
    try:
        source = connector()
    except Exception:
        return DatasetStatus.UNREACHABLE
    # A partial cache is not a usable source: every entry the source declares
    # has to be present and its manifest has to match.
    return (
        DatasetStatus.AVAILABLE if source.is_cached() else DatasetStatus.UNREACHABLE
    )


def resolve_dataset_statuses(
    datasets: list[RequiredDataset],
) -> list[RequiredDataset]:
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
        f"\n  availability score: {score:.3f}  expected 0.857"
        f"   ({len(available)} of {len(blocking)} blocking)"
    )
    missing = [d.name for d in blocking if d.status is not DatasetStatus.AVAILABLE]
    print(f"  outstanding: {', '.join(missing) or 'none'}")
