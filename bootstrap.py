"""Get a fresh clone to a working cache."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

ARCHIVE_NAME = "cart-platform-cache.tar.gz"
CHECKSUM_SUFFIX = ".sha256"


# data-v1 is left published so anyone who provisioned from it can tell
# which payload they have. It is missing eight malignant-cell digests
# -- including the one a standard run of each indication asks for, which
# costs a 2.6 GB rebuild on first run -- and the anti-tag binder.
# Replacing it in place would have been silent staleness: same name,
# different bytes, no way for a reader to know which they hold.
RELEASE_TAG = "data-v2"
PREVIOUS_TAGS = ("data-v1",)
RELEASE_REPO = "kakon2002/cart-platform"


HEAVY = ("*.h5ad", "*.h5ad.gz")

# A download that failed part way leaves one of these behind. The
# payload beside it is correct -- an unfinished file is never renamed
# over a finished one -- but a leftover was being counted in the size
# report and would have been packaged into a handover archive, so the
# claim that it is never mistaken for a finished file was true of the
# data and false of the accounting.
PARTIAL = "*.partial"


BUILD_ONLY = {"singlecell": {"archive", "matrix"}}


POOL_DERIVED = {
    "trials": "keyed by the screened antigen list; built during the first run",
}


SHARED_SOURCES = [
    ("uniprot", "the reviewed human proteome, 20,431 entries", "~10 min"),
    ("hpa", "normal tissue, pathology, subcellular, protein atlas", "11 s"),
    ("gtex", "bulk normal medians", "~2 min"),
    ("depmap", "CRISPR gene effect, 432 MB", "~5 min"),
    ("genespan", "gene annotation", "~1 min"),
    ("antibodies", "structure summary and therapeutics", "~1 min"),
    ("domains", "construct part sequences", "~1 min"),
    ("trials", "trial counts per antigen", "during the first run"),
    ("antitag", "the anti-tag binder every adaptor design is built on", "~5 s"),
]


PER_INDICATION = [
    ("tcga", "the tumour cohort, through the GDC API", "~3-21 min"),
    ("singlecell", "derived summaries; the atlas builds them", "minutes to hours"),
    ("depmap", "the per-lineage dependency matrix", "~1 min"),
]

SOURCES = SHARED_SOURCES


def _sha256(path: Path) -> str:
    """One implementation, used by both the writer and the verifier."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _state(name: str) -> tuple[str, int, list[str]]:
    """Whether a source is usable, by checking payloads and not just manifests."""
    directory = DATA / name
    if not directory.exists():
        return "MISSING", 0, []
    size = sum(f.stat().st_size for f in directory.rglob("*")
               if f.is_file() and not f.match(PARTIAL))
    found = sorted(directory.glob("*.manifest.json"))
    if not found:
        return "MISSING", size, []

    build_only = BUILD_ONLY.get(name, set())
    broken = []
    for path in found:
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            broken.append(f"{path.name} is unreadable")
            continue
        if meta.get("key") in build_only:
            continue
        filename = meta.get("filename")
        if filename and not (directory / filename).exists():
            broken.append(f"{filename} is named by a manifest but absent")
    if broken:
        return "BROKEN", size, broken
    return "present", size, []


def _indications():
    """The registry, or an empty list if the package cannot be imported."""
    try:
        from car_pipeline.configs.registry import INDICATIONS
        return sorted(INDICATIONS.values(), key=lambda i: i.key)
    except Exception:
        return []


def _tagged_present(namespace: str, tag: str) -> tuple[bool, int]:
    """Whether any payload in a namespace carries this indication's tag."""
    directory = DATA / namespace
    if not directory.exists() or not tag:
        return False, 0
    hits = [f for f in directory.glob(f"*__{tag}*")
            if f.is_file() and not f.name.endswith(".manifest.json")]
    return bool(hits), sum(f.stat().st_size for f in hits)


def stale_partials() -> list[Path]:
    """Unfinished downloads left by an attempt that did not complete."""
    if not DATA.exists():
        return []
    return [p for p in sorted(DATA.rglob("*"))
            if p.is_file() and p.match(PARTIAL)]


def report() -> int:
    """What is on this machine. Reads the cache only — never the network."""
    leftovers = stale_partials()
    if leftovers:
        print(f"  {len(leftovers)} unfinished download(s), left by an "
              f"attempt that did not complete:")
        for path in leftovers[:4]:
            print(f"    {path.relative_to(ROOT).as_posix()}  "
                  f"{path.stat().st_size / 1e6:.1f} MB")
        print("    They are not counted in the sizes below and are "
              "never packaged. Re-running --from-release replaces "
              "them.")
        print()
    print(f"cache root: {DATA}")
    print()
    usable = 0
    missing = []
    for name, description, cost in SOURCES:
        mark, size, problems = _state(name)
        pending = name in POOL_DERIVED
        if mark == "MISSING" and pending:
            mark = "deferred"
        if mark in ("present", "deferred"):
            usable += 1
        else:
            missing.append(name)
        print(f"  {mark:8s} {name:12s} {size / 1e6:8.1f} MB  {description}")
        for problem in problems:
            print(f"           {'':12s} {'':8s}     {problem}")
        if mark == "deferred":
            print(f"           {'':12s} {'':8s}     {POOL_DERIVED[name]}")
        elif mark == "MISSING":
            print(f"           {'':12s} {'':8s}     rebuild cost {cost}")

    indications = _indications()
    incomplete: list[str] = []
    if indications:
        print()
        print("  per indication:")
        for ind in indications:
            tags = [
                ("tcga", ind.tcga_project),
                ("singlecell", ind.atlas.series if ind.atlas else None),
                ("depmap", ind.depmap_lineage),
            ]
            parts = []
            for namespace, tag in tags:
                if not tag:
                    parts.append(f"{namespace}=declared-absent")
                    continue
                ok, size = _tagged_present(namespace, tag)
                parts.append(f"{namespace}={'ok' if ok else 'MISSING'}")
                if not ok:
                    incomplete.append(f"{ind.cancer_type}: {namespace} ({tag})")
            print(f"    {ind.cancer_type:34s} {'  '.join(parts)}")

    print()
    print(f"  {usable}/{len(SOURCES)} shared sources usable")
    if incomplete:
        print(f"  {len(incomplete)} per-indication cache(s) missing:")
        for row in incomplete:
            print(f"    {row}")
        missing = missing + incomplete
    if BUILD_ONLY:
        print("  the 8.3 GB matrix and its 2.6 GB archive are build-time only "
              "and are not expected here")
    if not missing:
        print("\nThe cache is complete. Nothing to do.")
        return 0
    print(f"\nMissing: {', '.join(missing)}")
    print("Run with --from-release (minutes) or --from-sources (~3 hours).")
    return 1


def _members() -> list[Path]:
    """Every cache file the archive should carry, heavy artifacts excluded."""
    return [p for p in sorted(DATA.rglob("*"))
            if p.is_file() and not any(p.match(pat) for pat in HEAVY)
            and not p.match(PARTIAL)]


def required_digests() -> dict:
    """The malignant-cell digest a standard run of each indication asks for.

    The single-cell malignant-cell cache is keyed by a digest of the exact gene
    set the run requests. That is a content key, so any change upstream of the
    pool silently invalidates a shipped artifact: the lookup misses, the source
    rebuilds from raw, and a 2.5 MB file costs a 2.6 GB download and an 8.3 GB
    expansion. Nothing checked it, which is why a release shipped seven digests
    and none of the one a standard run produces.

    Computing this needs the screening stages, so it is imported here rather
    than at module scope -- bootstrap must stay runnable before dependencies
    are installed.
    """
    from car_pipeline.configs.registry import INDICATIONS, PROJECTS
    from car_pipeline.data.coverage import build_coverage
    from car_pipeline.data.depmap import DepMapSource, gene_index
    from car_pipeline.data.gtex import GTExSource
    from car_pipeline.data.hpa import HPASource, index as atlas_index
    from car_pipeline.data.singlecell import (
        SingleCellSource, _gene_digest, match_surface as cell_match)
    from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
    from car_pipeline.data.uniprot import load_surface
    from car_pipeline.stages import stage3, stage4
    from car_pipeline.stages.stage1 import build_spec

    surface, _ = load_surface()
    by_acc, by_sym = atlas_index(HPASource().load())
    gtex_profiles, gtex_tissues, _ = GTExSource().match_surface(surface, by_acc)

    out = {}
    for key, indication in sorted(INDICATIONS.items()):
        if indication.atlas is None:
            continue
        project = PROJECTS[key]
        spec = build_spec(project)
        ceiling = spec.design_constraints.normal_tissue_risk_ceiling
        overrides = {o: ov.tier for o, ov
                     in spec.inputs.tissue_criticality_overrides.items()}
        cohort = TCGASource(indication.tcga_project).load()
        cohort_join = tcga_match(cohort, surface, by_acc)
        atlas = SingleCellSource(indication.atlas).load()
        cell_index = cell_match(atlas, surface, by_acc)
        dependency, _ = DepMapSource(indication.depmap_lineage).load()
        rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, cohort_join)
        calibration = stage3.calibrate_atlas_levels(
            surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
            stage3.RiskModel(overrides=overrides))
        ranked, _model, _ = stage3.rank(
            rows, {r.accession: r for r in surface}, by_acc, by_sym, atlas,
            cell_index, gtex_profiles, gtex_tissues, cohort, cohort_join,
            dependency, gene_index(dependency), overrides, ceiling, calibration,
            margin_label=indication.gtex_bulk_label)
        pool = stage4.build_pool(ranked)
        genes = sorted({r.gene for r in pool})
        out[indication.cancer_type] = {
            "slug": indication.atlas.slug,
            "digest": _gene_digest(genes),
            "genes": len(genes),
        }
    return out


def check_digests() -> tuple[bool, list[str]]:
    """Whether the cache holds the digest each indication's run will ask for."""
    try:
        wanted = required_digests()
    except Exception as exc:
        return False, [f"could not compute the required digests: "
                       f"{type(exc).__name__}: {exc}"]
    problems = []
    for cancer_type, want in sorted(wanted.items()):
        path = (DATA / "singlecell"
                / f"malignant_cells__{want['slug']}_{want['digest']}.npz")
        if not path.exists():
            have = sorted(
                q.name.split("_")[-1].removesuffix(".npz")
                for q in (DATA / "singlecell").glob(
                    f"malignant_cells__{want['slug']}_*.npz"))
            problems.append(
                f"{cancer_type}: a standard run asks for digest "
                f"{want['digest']} over {want['genes']} genes, and the cache "
                f"holds {len(have)} digest(s) for this atlas, none of them "
                f"that one")
    return not problems, problems


def package(destination: Path) -> int:
    """Build the archive to hand to whoever is deploying next."""
    if not DATA.exists():
        print(f"No cache at {DATA}; nothing to package.")
        return 1
    print("checking the cache's content keys against what a run will ask for")
    ok, problems = check_digests()
    if not ok:
        for detail in problems:
            print(f"  STALE  {detail}")
        print()
        print("  REFUSING to package. A cache whose malignant-cell digest does")
        print("  not match the pool the code produces looks complete and is")
        print("  not: the lookup misses, the source rebuilds from raw, and")
        print("  whoever unpacks this waits for a 2.6 GB download and an")
        print("  8.3 GB expansion to regenerate a file of about 2.5 MB.")
        print()
        print("  Run the pipeline once against this cache so the digest is")
        print("  built, then package again.")
        return 2
    print("  every registered indication's digest is present")

    files = _members()
    raw = sum(f.stat().st_size for f in files)
    print(f"packaging {len(files)} files, {raw / 1e6:.0f} MB uncompressed")
    print(f"  excluding {', '.join(HEAVY)} — build-time inputs, not serve-time")

    started = time.monotonic()
    with tarfile.open(destination, "w:gz", compresslevel=6) as tar:
        for index, path in enumerate(files, 1):
            tar.add(path, arcname=path.relative_to(ROOT).as_posix())
            if index % 10 == 0 or index == len(files):
                print(f"  {index}/{len(files)}", end="\r", flush=True)
    print()

    checksum = Path(str(destination) + CHECKSUM_SUFFIX)
    digest = _sha256(destination)
    checksum.write_text(f"{digest}  {destination.name}\n", encoding="utf-8")
    print(f"  {destination.name}  {destination.stat().st_size / 1e6:.0f} MB  "
          f"in {time.monotonic() - started:.0f}s")
    print(f"  {checksum.name}  {digest}")
    return 0


def from_archive(source: Path, require_checksum: bool = False) -> int:
    """Unpack a prepared cache, refusing one that fails its checksum."""
    if not source.exists():
        print(f"No archive at {source}.")
        print("Get it from whoever handed you this repository, or build the")
        print("cache yourself with --from-sources.")
        return 1

    checksum = Path(str(source) + CHECKSUM_SUFFIX)
    expected = ""
    if checksum.exists():
        parts = checksum.read_text(encoding="utf-8").split()
        expected = parts[0] if parts else ""
        if not expected:
            print(f"{checksum.name} is empty or malformed; treating as a mismatch.")
            return 2
    elif require_checksum:
        print(f"No {checksum.name} beside the archive.")
        print("  Refusing to unpack a downloaded archive unverified.")
        return 2
    else:
        print(f"No {checksum.name} beside the archive; unpacking unverified.")

    if expected:
        actual = _sha256(source)
        if actual != expected:
            print(f"CHECKSUM MISMATCH for {source.name}")
            print(f"  expected {expected}")
            print(f"  got      {actual}")
            print("  The transfer is incomplete or corrupt. Do not use it.")
            return 2
        print(f"checksum ok: {expected[:16]}...")

    started = time.monotonic()
    with tarfile.open(source, "r:gz") as tar:
        tar.extractall(ROOT, filter="data")
    print(f"unpacked in {time.monotonic() - started:.0f}s")
    print()
    return report()


def from_release() -> int:
    """Download the packaged cache from the repository release, then unpack."""
    if shutil.which("gh") is None:
        print("The GitHub CLI is not installed, and the release asset is")
        print("private so it cannot be fetched anonymously.")
        print("  Install:  https://cli.github.com")
        print("  Or download the asset by hand from")
        print(f"    https://github.com/{RELEASE_REPO}/releases/tag/{RELEASE_TAG}")
        print("  then:  python bootstrap.py --from-archive <path>")
        return 1

    # No size is printed. The literal here read 298 MB while the asset
    # was 471.7 MB, and a number nothing checks is how this repository
    # keeps finding stale claims. The download tool reports progress.
    print(f"downloading {ARCHIVE_NAME} from {RELEASE_TAG}")
    result = subprocess.run(
        ["gh", "release", "download", RELEASE_TAG, "-R", RELEASE_REPO,
         "--pattern", f"{ARCHIVE_NAME}*", "--clobber"],
        cwd=ROOT,
    )
    if result.returncode != 0:
        # This asserted an authentication cause it had no way to know,
        # and once sent a reader checking credentials for what was a
        # transient network failure. Report the fact, offer the retry.
        print(f"\nDownload failed: gh exited {result.returncode}. Any error it\n"
              "printed is above this line.")
        print("  --from-release has failed transiently twice in five measured runs.")
        print("  Re-running is safe: an unfinished file is written beside the real")
        print("  one and never renamed over it, so a failed attempt cannot corrupt")
        print("  a good cache. Try the same command again.")
        print("  If it fails twice the same way, check the release exists")
        print(f"    gh release view {RELEASE_TAG} -R {RELEASE_REPO}")
        print("  and that you are authenticated")
        print("    gh auth status")
        return 1

    return from_archive(ROOT / ARCHIVE_NAME, require_checksum=True)


def from_sources() -> int:
    """Rebuild every cache from its origin. Long, and entirely automated."""
    print("Rebuilding from the original sources.")
    print("  Every source has a public programmatic URL. No accounts, no")
    print("  manual downloads, no registration.")
    print()
    print("  Expect roughly 3 hours and ~12 GB of free disk. Almost all of it")
    print("  is one step: deriving the single-cell group means streams the")
    print("  whole 8.3 GB matrix and took 2 h 19 min on the machine that first")
    print("  built this cache. Every download combined was about 15 minutes.")
    print()
    print("  Interrupting is safe. Each artifact is committed with a manifest")
    print("  only once it is complete, so a re-run resumes rather than")
    print("  restarts, and a partial file is never mistaken for a finished one.")
    print()

    try:
        from car_pipeline.data.antibodies import AntibodySource
        from car_pipeline.data.depmap import DepMapSource
        from car_pipeline.data.domains import DomainSource
        from car_pipeline.data.genespan import GeneSpanSource
        from car_pipeline.data.gtex import GTExSource
        from car_pipeline.data.hpa import HPASource
        from car_pipeline.data.singlecell import SingleCellSource
        from car_pipeline.data.tcga import TCGASource
        from car_pipeline.data.uniprot import UniProtSource
    except ImportError as exc:
        print(f"Cannot import the pipeline: {exc}")
        print("  This path needs the project's dependencies, unlike the others.")
        print("  Use the interpreter in .venv:")
        print(r"    .venv\Scripts\python.exe bootstrap.py --from-sources"
              "        (Windows)")
        print("    .venv/bin/python bootstrap.py --from-sources"
              "             (macOS, Linux)")
        return 1

    steps = [
        ("hpa", lambda: HPASource().fetch()),
        ("gtex", lambda: GTExSource().fetch()),
        ("genespan", lambda: GeneSpanSource().fetch()),
        ("domains", lambda: DomainSource().fetch()),
        ("antibodies", lambda: AntibodySource().fetch()),
        ("uniprot", lambda: UniProtSource().fetch()),
    ]

    for ind in _indications():
        label = ind.key
        if ind.tcga_project:
            steps.append((f"tcga/{label}",
                          lambda i=ind: TCGASource(i.tcga_project).build_cohort()))
        if ind.depmap_lineage:
            steps.append((f"depmap/{label}",
                          lambda i=ind: DepMapSource(i.depmap_lineage).build_matrix()))
        if ind.atlas:
            steps.append((f"singlecell/{label}",
                          lambda i=ind: SingleCellSource(i.atlas).build_group_means()))

    overall = time.monotonic()
    for name, run in steps:
        print(f"=== {name} ===", flush=True)
        started = time.monotonic()
        try:
            run()
        except KeyboardInterrupt:
            print(f"\n  interrupted during {name}. Nothing is corrupted: an "
                  "artifact\n  is committed only once complete. Re-run to "
                  "resume.")
            return 130
        except Exception as exc:
            print(f"  FAILED after {time.monotonic() - started:.0f}s: "
                  f"{type(exc).__name__}: {exc}")
            print(traceback.format_exc())
            print("  Re-running resumes from what completed.")
            return 2
        print(f"  done in {time.monotonic() - started:.0f}s", flush=True)
    print(f"\nall sources in {(time.monotonic() - overall) / 60:.0f} min")
    print()
    print("The trial counts and the malignant-cell summaries are keyed by the")
    print("screened pool, and there is one set per indication, so they are built")
    print("by the first run of each rather than here:")
    print("    <the interpreter in .venv> run_all.py --fresh")
    print()
    return report()


def main() -> int:
    """Parse the provisioning mode and run it."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--from-release", action="store_true",
                       help="download the prepared cache and unpack it")
    group.add_argument("--from-archive", nargs="?", const=ARCHIVE_NAME,
                       metavar="PATH", help="unpack a cache you already have")
    group.add_argument("--from-sources", action="store_true",
                       help="rebuild every cache from its origin (~3 hours; "
                            "needs the project interpreter)")
    group.add_argument("--package", nargs="?", const=ARCHIVE_NAME,
                       metavar="PATH", help="build the archive to hand on")
    args = parser.parse_args()

    if args.from_release:
        return from_release()
    if args.package:
        return package(Path(args.package))
    if args.from_archive:
        return from_archive(Path(args.from_archive))
    if args.from_sources:
        return from_sources()
    return report()


if __name__ == "__main__":
    raise SystemExit(main())
