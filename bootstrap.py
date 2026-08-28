"""Get a fresh clone to a working cache.

`data/` is not in git — it is 680 MB of cached sources and none of it is
tracked. A clone has an empty cache, and nothing runs against an empty cache.
This closes that gap two ways.

    python bootstrap.py                 # what is present, what is missing
    python bootstrap.py --from-release  # download a prepared cache  (minutes)
    python bootstrap.py --from-archive  # unpack one you already have
    python bootstrap.py --package       # make the archive to hand on

    .venv/bin/python bootstrap.py --from-sources   # rebuild from origin, ~3 h

**Prefer the release.** Rebuilding is fully automated and every source has a
public programmatic URL, so it needs no accounts and no manual downloads — but
one step in it is unavoidably long. Deriving `group_means.npz` streams the whole
8.3 GB single-cell matrix and took **2 h 19 min** on the machine that first
built this cache, against about 15 minutes for every download combined. The
archive skips that by carrying the 5.7 MB result.

Everything here is standard library **except** `--from-sources`, which imports
the pipeline and therefore needs the project interpreter with its dependencies.
The other paths run under any Python 3, before anything is installed.
"""

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

#: Where the packaged cache lives. The repository is private, so the asset is
#: private with it: anyone who can clone can download, and nobody else can.
RELEASE_TAG = "data-v1"
RELEASE_REPO = "kakon2002/cart-platform"

#: Excluded from the archive, and the single authority on that. `_members()`
#: iterates this rather than repeating the patterns, so what is printed as
#: skipped is what is actually skipped.
HEAVY = ("*.h5ad", "*.h5ad.gz")

#: Cache keys that exist only to build something else. Their payloads are
#: deliberately absent from a deployed cache — a served run never opens either,
#: which is measured rather than assumed — so a missing payload here is the
#: intended state and is reported as such instead of as a gap.
BUILD_ONLY = {"singlecell": {"archive", "matrix"}}

#: Sources that cannot be fetched ahead of a run, with the reason. The trial
#: cache is fingerprinted by the antigen list it covers, so it has no meaning
#: until a pool exists; it is built during the first screen. Listing it as
#: merely missing would report a complete cache as incomplete forever.
POOL_DERIVED = {
    "trials": "keyed by the screened antigen list; built during the first run",
}

#: Every source, with what it holds and roughly what a cold rebuild costs.
#: The times are measured from the manifests of the original build rather than
#: estimated: each manifest records `retrieved_at`, and these are the gaps
#: between consecutive entries.
#: Sources that describe the human body. One copy, whatever is being screened.
SHARED_SOURCES = [
    ("uniprot", "the reviewed human proteome, 20,431 entries", "~10 min"),
    ("hpa", "normal tissue, pathology, subcellular, protein atlas", "11 s"),
    ("gtex", "bulk normal medians", "~2 min"),
    ("depmap", "CRISPR gene effect, 432 MB", "~5 min"),
    ("genespan", "gene annotation", "~1 min"),
    ("antibodies", "structure summary and therapeutics", "~1 min"),
    ("domains", "construct part sequences", "~1 min"),
    ("trials", "trial counts per antigen", "during the first run"),
]

#: Sources that describe a tumour. One copy PER INDICATION, which is why the
#: report and the rebuild both iterate the registry rather than a fixed list.
#: The old flat list named "tcga" and "singlecell" as though each were a single
#: thing, so a clone provisioned the reference indication and reported itself
#: complete while a second indication had nothing.
PER_INDICATION = [
    ("tcga", "the tumour cohort, through the GDC API", "~3-21 min"),
    ("singlecell", "derived summaries; the atlas builds them", "minutes to hours"),
    ("depmap", "the per-lineage dependency matrix", "~1 min"),
]

SOURCES = SHARED_SOURCES


def _sha256(path: Path) -> str:
    """One implementation, used by both the writer and the verifier.

    Two copies would let the two sides diverge, and a divergence there produces
    a mismatch on a perfectly good archive with no way to tell which is wrong.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _state(name: str) -> tuple[str, int, list[str]]:
    """Whether a source is usable, by checking payloads and not just manifests.

    A manifest is the marker that a fetch finished, but its presence alone does
    not mean the file it describes is still there. Globbing for manifests would
    report a cache as complete while the payload beside it had been deleted,
    which is the failure this is most likely to be asked about.
    """
    directory = DATA / name
    if not directory.exists():
        return "MISSING", 0, []
    size = sum(f.stat().st_size for f in directory.rglob("*") if f.is_file())
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
            continue                      # absent on purpose; see BUILD_ONLY
        filename = meta.get("filename")
        if filename and not (directory / filename).exists():
            broken.append(f"{filename} is named by a manifest but absent")
    if broken:
        return "BROKEN", size, broken
    return "present", size, []


def _indications():
    """The registry, or an empty list if the package cannot be imported.

    bootstrap runs before dependencies are installed on a fresh clone, so it
    must degrade to the shared-source report rather than failing on an import.
    """
    try:
        from car_pipeline.configs.registry import INDICATIONS
        return sorted(INDICATIONS.values(), key=lambda i: i.key)
    except Exception:                                  # noqa: BLE001
        return []


def _tagged_present(namespace: str, tag: str) -> tuple[bool, int]:
    """Whether any payload in a namespace carries this indication's tag."""
    directory = DATA / namespace
    if not directory.exists() or not tag:
        return False, 0
    hits = [f for f in directory.glob(f"*__{tag}*")
            if f.is_file() and not f.name.endswith(".manifest.json")]
    return bool(hits), sum(f.stat().st_size for f in hits)


def report() -> int:
    """What is on this machine. Reads the cache only — never the network."""
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
    # Per indication. A clone that has the shared sources and one indication's
    # tumour caches is not ready for a registry that declares two, and saying
    # "10/10 sources usable" would tell it that it is.
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
    return [p for p in sorted(DATA.rglob("*"))
            if p.is_file() and not any(p.match(pat) for pat in HEAVY)]


def package(destination: Path) -> int:
    """Build the archive to hand to whoever is deploying next."""
    if not DATA.exists():
        print(f"No cache at {DATA}; nothing to package.")
        return 1
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
        # An unparseable sidecar is treated as a mismatch, not as a crash. A
        # truncated checksum file is exactly the state this machinery exists to
        # survive, so it must not be the thing that raises.
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
            # Refused rather than unpacked. A truncated transfer produces a
            # cache that reads as present and answers with the wrong data,
            # which is worse than having no cache at all.
            print(f"CHECKSUM MISMATCH for {source.name}")
            print(f"  expected {expected}")
            print(f"  got      {actual}")
            print("  The transfer is incomplete or corrupt. Do not use it.")
            return 2
        print(f"checksum ok: {expected[:16]}...")

    started = time.monotonic()
    with tarfile.open(source, "r:gz") as tar:
        # filter="data" refuses absolute paths and parent traversal in member
        # names. No fallback: this project runs on 3.13, and an unpacker that
        # silently drops the guard on older interpreters is worse than one that
        # refuses to run there.
        tar.extractall(ROOT, filter="data")
    print(f"unpacked in {time.monotonic() - started:.0f}s")
    print()
    return report()


def from_release() -> int:
    """Download the packaged cache from the repository release, then unpack.

    Uses the GitHub CLI because the repository is private and the asset
    inherits that: an unauthenticated fetch gets a 404, which is the right
    answer but an unhelpful one, so the credential the clone already needed is
    reused rather than a token being invented.
    """
    if shutil.which("gh") is None:
        print("The GitHub CLI is not installed, and the release asset is")
        print("private so it cannot be fetched anonymously.")
        print("  Install:  https://cli.github.com")
        print("  Or download the asset by hand from")
        print(f"    https://github.com/{RELEASE_REPO}/releases/tag/{RELEASE_TAG}")
        print("  then:  python bootstrap.py --from-archive <path>")
        return 1

    print(f"downloading {ARCHIVE_NAME} (298 MB) from {RELEASE_TAG}")
    result = subprocess.run(
        ["gh", "release", "download", RELEASE_TAG, "-R", RELEASE_REPO,
         "--pattern", f"{ARCHIVE_NAME}*", "--clobber"],
        cwd=ROOT,
    )
    if result.returncode != 0:
        print("\nDownload failed. If this is an authentication error, run:")
        print("    gh auth login")
        return 1
    # require_checksum: a download is the one case where the sidecar must be
    # there. Unpacking 298 MB off the network unverified is the failure mode the
    # checksum exists for.
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

    # `fetch()` is the method that populates a cache; several sources have no
    # `load()` at all (GTEx exposes only match_surface, which needs the surface
    # proteome). The two that derive an artifact from a download — TCGA's cohort
    # and the single-cell group means — are named explicitly, because the
    # download alone leaves the expensive part undone.
    #
    # Ordered cheapest-first so a broken network fails in seconds rather than
    # after the 2.6 GB download.
    steps = [
        ("hpa", lambda: HPASource().fetch()),
        ("gtex", lambda: GTExSource().fetch()),
        ("genespan", lambda: GeneSpanSource().fetch()),
        ("domains", lambda: DomainSource().fetch()),
        ("antibodies", lambda: AntibodySource().fetch()),
        ("uniprot", lambda: UniProtSource().fetch()),
    ]

    # Then the tumour-side caches, once per registered indication. This used to
    # be three unparameterised calls, so a clone provisioned the reference
    # indication only and the multi-indication stage failed on a fresh machine
    # while every other stage passed.
    for ind in _indications():
        label = ind.key
        if ind.tcga_project:
            steps.append((f"tcga/{label}",
                          lambda i=ind: TCGASource(i.tcga_project).build_cohort()))
        if ind.depmap_lineage:
            steps.append((f"depmap/{label}",
                          lambda i=ind: DepMapSource(i.depmap_lineage).build_matrix()))
        if ind.atlas:
            # For the reference indication this downloads 2.6 GB, expands it to
            # 8.3 GB and streams the whole thing: the ~2 h 19 min step. A
            # CELLxGENE export needs no expansion and takes seconds.
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
        except Exception as exc:                       # noqa: BLE001
            # The traceback, not just the message. A failure two hours into the
            # last step with one line of context means re-running two hours to
            # find out where it was.
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
