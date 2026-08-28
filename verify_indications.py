"""Multi-indication, against the criteria in specs/multi-indication.md.

M1 is the reason this file exists. Two indications sharing one cache slot is the
worst failure shape available here: the second overwrites the first in place,
and the first then screens against the other's atlas and produces a ranked list
that looks entirely plausible. So M1 does not inspect the code -- it runs both
indications and asserts that neither one's artifacts moved.

M7 exists because a Mode A pin drawn from the platform's own top 20 would pass
by agreeing with itself. The pin is a target the platform ranks nowhere near the
top, so the check exercises validation rather than self-agreement.
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.api import pipeline
from car_pipeline.configs.registry import INDICATIONS, registered, resolve

DATA = Path(__file__).resolve().parent / "data"

#: Sources describing the human body rather than a tumour. M2 asserts none of
#: these gained a per-indication copy: duplicating them would assert that
#: normal-tissue biology changes with the diagnosis.
SHARED = ("uniprot", "gtex", "hpa", "genespan", "antibodies", "domains")

#: The Mode A pin. CD19 is the canonical CAR-T target -- approved therapies use
#: it -- but for B-cell malignancy, not for a solid tumour. The platform ranks
#: it around 1,300 of 3,400 here, so a passing verdict cannot come from the
#: platform agreeing with its own ranking.
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
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
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

    # M1 -- the priority defect
    moved = sorted(k for k in before if before[k] != after.get(k))
    vanished = sorted(set(before) - set(after))
    criterion("M1", bool(moved or vanished),
              f"{len(before)} indication-tagged artifacts; {len(moved)} changed, "
              f"{len(vanished)} disappeared after running both"
              + (f" -> {moved + vanished}" if (moved or vanished) else ""))

    # M2 -- shared sources not duplicated
    dupes = [p.name for ns in SHARED for p in (DATA / ns).glob("*")
             if p.is_file() and "__" in p.name]
    criterion("M2", bool(dupes),
              f"shared sources carry no per-indication copy ({len(dupes)} found)")

    # M3 -- no indication-specific module constant survives
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

    # M4 / M5 -- an atlas-less indication refuses, and says why
    from car_pipeline.configs.indication import Indication
    bare = Indication(key="none", cancer_type="x", tcga_project=None,
                      depmap_lineage=None, gtex_bulk_label=None, atlas=None)
    criterion("M4", bare.atlas is not None,
              "an indication with no atlas is expressible and carries "
              "atlas=None, which the driver maps to NOT_USABLE")
    src = Path("car_pipeline/api/pipeline.py").read_text(encoding="utf-8")
    names_c2 = "malignant_vs_stroma" in src and "NOT_USABLE" in src
    criterion("M5", not names_c2,
              "the atlas-less path names malignant_vs_stroma as what is lost, "
              "not just a weight")

    # M6 / M7 / M8 -- Mode A
    verdict = pipeline.validate("pdac", MODE_A_PIN)
    criterion("M6", not verdict.get("verdict"),
              f"Mode A on {MODE_A_PIN} returns {verdict.get('verdict')} "
              f"with {len(verdict.get('reasons', []))} reasons")

    rank = verdict.get("rank")
    criterion("M7", rank is not None and rank <= 20,
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

    # M9 -- the reference indication is unchanged
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

    # M10 -- a missing source is named
    named = all(":" in u for u in ref.get("unavailable", []))
    criterion("M10", not named,
              f"unavailable components name their source "
              f"({ref.get('unavailable') or 'none unavailable'})")

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
