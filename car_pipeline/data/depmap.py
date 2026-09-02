"""Cancer dependency screens for the indication lineage."""

from __future__ import annotations

import csv
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from car_pipeline.data.source import CacheEntry, DataSource, stream_to_file

RELEASE_PIN = "24Q4"
ARTICLE_ID = 27993248
ARTICLE_ENDPOINT = f"https://api.figshare.com/v2/articles/{ARTICLE_ID}/files?page_size=500"
DOWNLOAD_HOST = "https://ndownloader.figstatic.com"

WANTED_FILES = {
    "gene_effect": "CRISPRGeneEffect.csv",
    "model": "Model.csv",
}


DEFAULT_LINEAGE = "Pancreas"
DEPENDENCY_THRESHOLD = -0.5


@dataclass
class Dependency:
    genes: np.ndarray
    models: np.ndarray
    diseases: np.ndarray
    values: np.ndarray

    def median_effect(self, gene_index: int) -> float:
        """Median gene effect across screened lines, NaN if none was screened."""
        column = self.values[:, gene_index]
        if np.all(np.isnan(column)):
            return float("nan")
        return float(np.nanmedian(column))

    def dependent_lines(self, gene_index: int) -> int:
        """How many lines sit at or below the dependency threshold."""
        column = self.values[:, gene_index]
        return int(np.sum(column <= DEPENDENCY_THRESHOLD))

    def screened_lines(self, gene_index: int) -> int:
        """How many lines carry a measurement for this gene."""
        return int(np.sum(~np.isnan(self.values[:, gene_index])))


class DepMapSource(DataSource):
    name = "DepMap"
    namespace = "depmap"

    def __init__(self, lineage: str | None = None, root=None) -> None:
        """Bind this source to one lineage."""
        super().__init__(root=root)
        self.lineage = lineage or DEFAULT_LINEAGE

    def cache_entries(self) -> Iterable[CacheEntry]:
        """The release files and the per-lineage matrix derived from them."""
        return [
            CacheEntry(
                key=key,
                filename=name,
                fingerprint={"release": RELEASE_PIN, "file": name},
            )
            for key, name in WANTED_FILES.items()
        ] + [
            CacheEntry(
                key=f"lineage_matrix__{self.lineage}",
                filename=f"lineage_dependency__{self.lineage}.npz",
                fingerprint={"release": RELEASE_PIN, "lineage": self.lineage},
            )
        ]

    def _resolve_by_name(self) -> dict[str, dict]:
        """Map wanted file names to a download target and published checksum."""
        req = urllib.request.Request(
            ARTICLE_ENDPOINT, headers={"User-Agent": "cart-platform/1.0"}
        )
        with urllib.request.urlopen(req, timeout=180) as resp:
            listing = json.loads(resp.read())
        by_name = {
            f["name"]: {
                "url": f"{DOWNLOAD_HOST}/files/{f['id']}",
                "md5": f.get("computed_md5") or f.get("supplied_md5"),
                "size": f.get("size"),
            }
            for f in listing
        }
        missing = [n for n in WANTED_FILES.values() if n not in by_name]
        if missing:
            raise FileNotFoundError(
                f"release {RELEASE_PIN} does not list: {', '.join(missing)}"
            )
        return by_name

    def fetch(self) -> None:
        """Download the release files this source is missing."""
        urls: dict[str, str] | None = None
        for entry in self.cache_entries():
            if entry.key.startswith("lineage_matrix__") or self.cache.is_valid(entry):
                continue
            if urls is None:
                urls = self._resolve_by_name()

            def fetcher(tmp: Path, entry: CacheEntry = entry, urls=urls) -> dict:
                """Fetch one artifact into a temporary file."""
                target = urls[entry.filename]
                print(f"  fetching {entry.filename}", flush=True)
                return stream_to_file(
                    target["url"],
                    tmp,
                    progress_label=entry.filename,
                    expected_md5=target["md5"],
                )

            self.cache.ensure(entry, fetcher)

    def path_for(self, key: str) -> Path:
        """The cached path for one key, fetching first if it is absent."""
        entry = next(e for e in self.cache_entries() if e.key == key)
        if not self.cache.is_valid(entry):
            self.fetch()
        return self.cache.path(entry)

    def lineage_models(self) -> dict[str, str]:
        """Model identifier to primary disease, restricted to the lineage."""
        out: dict[str, str] = {}
        with open(self.path_for("model"), encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("OncotreeLineage") == self.lineage:
                    out[row["ModelID"]] = row.get("OncotreePrimaryDisease", "")
        return out

    def build_matrix(self) -> Path:
        """Derive the per-lineage gene-effect matrix from the release files."""
        entry = next(e for e in self.cache_entries()
                     if e.key == f"lineage_matrix__{self.lineage}")

        def fetcher(tmp: Path) -> dict:
            """Fetch one artifact into a temporary file."""
            wanted = self.lineage_models()
            print(f"  {len(wanted)} models in the lineage", flush=True)

            with open(
                self.path_for("gene_effect"), encoding="utf-8", newline=""
            ) as fh:
                reader = csv.reader(fh)
                header = next(reader)

                genes = [c.split(" (")[0] for c in header[1:]]

                models: list[str] = []
                rows: list[np.ndarray] = []
                for row in reader:
                    model_id = row[0]
                    if model_id not in wanted:
                        continue

                    values = np.asarray(
                        [np.nan if c == "" else float(c) for c in row[1:]],
                        dtype=np.float32,
                    )
                    models.append(model_id)
                    rows.append(values)

            matrix = np.vstack(rows).astype(np.float32)
            diseases = np.asarray([wanted[m] for m in models])

            with open(tmp, "wb") as out:
                np.savez_compressed(
                    out,
                    genes=np.asarray(genes),
                    models=np.asarray(models),
                    diseases=diseases,
                    values=matrix,
                    lineage_total=np.asarray([len(wanted)]),
                )
            return {
                "observed_rows": matrix.shape[0],
                "extra": {
                    "lineage_models": len(wanted),
                    "screened_models": matrix.shape[0],
                    "genes": int(matrix.shape[1]),
                },
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> tuple[Dependency, int]:
        """The lineage matrix and the genes it is indexed by."""
        path = self.build_matrix()
        with np.load(path, allow_pickle=False) as data:
            dep = Dependency(
                genes=data["genes"],
                models=data["models"],
                diseases=data["diseases"],
                values=data["values"],
            )
            lineage_total = int(data["lineage_total"][0])
        return dep, lineage_total


def gene_index(dep: Dependency) -> dict[str, int]:
    """Symbol lookup: no stable gene identifier here, so the join is symbol-only."""
    out: dict[str, int] = {}
    for i, g in enumerate(dep.genes):
        out.setdefault(str(g), i)
    return out


if __name__ == "__main__":
    from car_pipeline.data.uniprot import load_surface

    dep, lineage_total = DepMapSource().load()
    idx = gene_index(dep)
    surface, _ = load_surface()

    matched = [r for r in surface if r.gene and r.gene in idx]
    strong = sum(
        1 for r in matched if dep.median_effect(idx[r.gene]) <= DEPENDENCY_THRESHOLD
    )
    adeno = int(np.sum(dep.diseases == "Pancreatic Adenocarcinoma"))
    adenosq = int(np.sum(dep.diseases == "Adenosquamous Carcinoma of the Pancreas"))

    counts = {
        "lineage models": (lineage_total, 68),
        "unscreened": (lineage_total - len(dep.models), 21),
        "screened lines": (len(dep.models), 47),
        "adenocarcinoma": (adeno, 45),
        "adenosquamous": (adenosq, 2),
        "genes": (len(dep.genes), 17916),
        "surface matched": (len(matched), 3169),
        "surface with no row": (len(surface) - len(matched), 327),
        "median effect <= -0.5": (strong, 17),
    }
    for label, (got, exp) in counts.items():
        pct = abs(got - exp) / exp * 100 if exp else 0.0
        print(f"  {label}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    print("\n  gene effect validation")
    checks = [
        ("KRAS", -1.758, 43),
        ("RAN", -4.076, 47),
        ("RPL13A", -1.721, None),
        ("CDK4", -0.681, None),
        ("EGFR", -0.359, None),
        ("BRAF", -0.123, None),
        ("SMAD4", -0.064, None),
        ("TP53", 0.157, None),
    ]
    for gene, exp_val, exp_lines in checks:
        i = idx.get(gene)
        if i is None:
            print(f"    {gene}: absent")
            continue
        val = dep.median_effect(i)
        line = f"    {gene:8s} {val:+7.3f}   expected {exp_val:+.3f}"
        if exp_lines is not None:
            line += (
                f"   dependent in {dep.dependent_lines(i)}/{dep.screened_lines(i)}"
                f"   expected {exp_lines}/47"
            )
        print(line)

    overall = float(np.nanmedian(dep.values))
    print(f"\n    median across all genes {overall:+.3f}   expected -0.034")
