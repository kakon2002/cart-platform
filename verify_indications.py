"""Multi-indication, against the criteria in specs/multi-indication.md."""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.api import pipeline
from car_pipeline.configs.registry import INDICATIONS, registered, resolve

DATA = Path(__file__).resolve().parent / "data"


SHARED = ("uniprot", "gtex", "hpa", "genespan", "antibodies", "domains")


MODE_A_PIN = "CD19"


def digests(prefixes: tuple[str, ...]) -> dict[str, str]:
    """sha256 of every cache payload whose name carries an indication tag."""
    out: dict[str, str] = {}
    for path in sorted(DATA.rglob("*")):
        if not path.is_file() or "__" not in path.name:
            continue
        if not any(p in path.name for p in prefixes):
            continue
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        out[f"{path.parent.name}/{path.name}"] = h.hexdigest()
    return out


def main() -> int:
    """Run the multi-indication criteria."""
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    print("registered indications:")
    for row in registered():
        print(f"  {row['cancer_type']:36s} cohort={row['cohort']:10s} "
              f"atlas={row['atlas']} lineage={row['dependency_lineage']}")

    tags = tuple(
        [i.tcga_project for i in INDICATIONS.values() if i.tcga_project]
        + [i.atlas.series for i in INDICATIONS.values() if i.atlas]
        + [i.depmap_lineage for i in INDICATIONS.values() if i.depmap_lineage]
    )

    print("\nrunning both indications", flush=True)
    before = digests(tags)
    results = {}
    for ind in sorted(INDICATIONS.values(), key=lambda i: i.key):
        print(f"  {ind.cancer_type}", flush=True)
        results[ind.key] = pipeline.run(ind.cancer_type,
                                        progress=lambda s, n="": None)
    after = digests(tags)

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)

    moved = sorted(k for k in before if before[k] != after.get(k))
    vanished = sorted(set(before) - set(after))
    criterion("M1", bool(moved or vanished),
              f"{len(before)} indication-tagged artifacts; {len(moved)} changed, "
              f"{len(vanished)} disappeared after running both"
              + (f" -> {moved + vanished}" if (moved or vanished) else ""))

    dupes = [p.name for ns in SHARED for p in (DATA / ns).glob("*")
             if p.is_file() and "__" in p.name]
    criterion("M2", bool(dupes),
              f"shared sources carry no per-indication copy ({len(dupes)} found)")

    import car_pipeline.data.singlecell as sc
    import car_pipeline.data.tcga as tc
    import car_pipeline.data.depmap as dm
    leftovers = [n for n, mod in (("SERIES", sc), ("ARCHIVE", sc),
                                  ("MALIGNANT_LEVEL1", sc), ("LEVEL1", sc),
                                  ("PROJECT", tc), ("LINEAGE", dm))
                 if hasattr(mod, n)]
    criterion("M3", bool(leftovers),
              f"indication-specific module constants remaining: "
              f"{leftovers or 'none'}")

    from car_pipeline.configs.indication import Indication
    bare = Indication(
        key="atlasless", cancer_type="Atlas-less Control",
        tcga_project="TCGA-PAAD", depmap_lineage="Pancreas",
        gtex_bulk_label="Pancreas", atlas=None)
    from car_pipeline.configs import registry as _reg
    _reg.INDICATIONS[bare.cancer_type.lower()] = bare
    _reg.PROJECTS[bare.cancer_type.lower()] = _reg.PROJECTS["pancreatic ductal adenocarcinoma"]
    try:
        bare_run = pipeline.run(bare.cancer_type, progress=lambda s, n="": None)
        crashed = None
    except Exception as exc:
        bare_run, crashed = None, f"{type(exc).__name__}: {exc}"
    finally:
        _reg.INDICATIONS.pop(bare.cancer_type.lower(), None)
        _reg.PROJECTS.pop(bare.cancer_type.lower(), None)

    criterion("M4", crashed is not None
              or bare_run.get("usability") != pipeline.NOT_USABLE
              or bare_run.get("final"),
              f"an atlas-less indication returns {(bare_run or {}).get('usability')} "
              f"with no ranking"
              + (f" -- it raised {crashed}" if crashed else ""))

    reasons_text = " ".join((bare_run or {}).get("reasons", []))
    criterion("M5", "malignant_vs_stroma" not in reasons_text,
              "the refusal names malignant_vs_stroma as the missing "
              "discriminator, not just a lost weight")

    verdict = pipeline.validate("pdac", MODE_A_PIN)
    criterion("M6", not verdict.get("verdict"),
              f"Mode A on {MODE_A_PIN} returns {verdict.get('verdict')} "
              f"with {len(verdict.get('reasons', []))} reasons")

    rank = verdict.get("rank")

    criterion("M7", rank is None or rank <= 20,
              f"{MODE_A_PIN} ranks {rank} of {verdict.get('of')} -- outside the "
              "top 20, so the verdict is not self-agreement")

    screen = results["pdac"]
    same = next((x for x in screen["ranked"] if x.gene == MODE_A_PIN), None)
    agree = (same is not None
             and same.risk == verdict.get("risk")
             and same.composite == verdict.get("composite"))
    criterion("M8", not agree,
              f"Mode A and Mode B report the same evidence for {MODE_A_PIN} "
              f"(risk {verdict.get('risk')}, composite {verdict.get('composite')})")

    ref = results["pdac"]
    ref_scored = sorted((x for x in ref["ranked"]
                         if x.composite is not None and x.gene),
                        key=lambda x: (-x.composite, x.gene, x.accession))
    from collections import Counter
    outcomes = Counter(d["outcome"] for d in ref["decisions"])
    expected_top = ["CEACAM5", "TMC5", "MUCL3"]
    expected_hash = "a91c696f2e1318f7"
    actual_top = [x.gene for x in ref_scored[:3]]
    criterion("M9",
              actual_top != expected_top
              or ref["stage3_hash"] != expected_hash
              or len(ref["pool"]) != 200,
              f"reference unchanged: top3 {actual_top}, pool {len(ref['pool'])}, "
              f"hash {ref['stage3_hash']}, outcomes {dict(outcomes)}")

    degraded = Indication(
        key="nodep", cancer_type="No-Dependency Control",
        tcga_project="TCGA-PAAD", depmap_lineage="NoSuchLineage",
        gtex_bulk_label="Pancreas", atlas=INDICATIONS["pancreatic ductal adenocarcinoma"].atlas)
    _reg.INDICATIONS[degraded.cancer_type.lower()] = degraded
    _reg.PROJECTS[degraded.cancer_type.lower()] = _reg.PROJECTS["pancreatic ductal adenocarcinoma"]
    try:
        deg = pipeline.run(degraded.cancer_type, progress=lambda s, n="": None)
        deg_unavailable = deg.get("unavailable", [])
    except Exception as exc:
        deg_unavailable = [f"RAISED {type(exc).__name__}: {exc}"]
    finally:
        _reg.INDICATIONS.pop(degraded.cancer_type.lower(), None)
        _reg.PROJECTS.pop(degraded.cancer_type.lower(), None)

    named = (bool(deg_unavailable)
             and all(":" in u for u in deg_unavailable)
             and not any(u.startswith("RAISED") for u in deg_unavailable)
             and all(":" in u for u in ref.get("unavailable", [])))
    criterion("M10", not named,
              f"a degraded indication names its missing source: "
              f"{deg_unavailable}")

    print("=" * 72)
    print(f"  {10 - len(tripped)}/10 criteria clear")

    print()
    print("=" * 72)
    print("WHAT EACH INDICATION PRODUCED")
    print("=" * 72)
    for key, r in sorted(results.items()):
        sc_ = sorted((x for x in r["ranked"] if x.composite is not None and x.gene),
                     key=lambda x: (-x.composite, x.gene, x.accession))
        oc = Counter(d["outcome"] for d in r["decisions"])
        print(f"  {r['indication'].cancer_type}")
        print(f"    usability {r['usability']}  status {r['status']}")
        print(f"    scored {len(sc_):,}  pool {len(r['pool'])}  "
              f"top {[x.gene for x in sc_[:5]]}")
        print(f"    outcomes {dict(oc)}")
        if r.get("unavailable"):
            for u in r["unavailable"]:
                print(f"    unavailable: {u}")
        print()

    if tripped:
        print(f"  STOPPING: {', '.join(tripped)} tripped.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
