"""Stage 4a against the criteria fixed in specs/stage4a-architecture-routing.md.

Criteria before biology. A5 and A6 are positive pins — the shape that caught the
Stage 5 join bug, where every check was a negative and all 200 passed against a
dead route. A9 and A10 exist because §3 fixes a number this pipeline cannot
measure: A9 makes its effect visible, A10 makes moving it after seeing output a
tripped criterion rather than an edit.
"""

from __future__ import annotations

import os
import sys

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.stages import routing, stage4
from car_pipeline.stages.stage1 import build_spec

#: The value recorded in the spec. A10 compares the config against this, so
#: moving the ceiling without moving the spec is a tripped criterion.
SPEC_TERMINABLE_CEILING = 0.35

#: A9 reports the admitted count across this range so the choice of ceiling is
#: visible rather than argued.
SWEEP = (0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70)


def main() -> int:
    spec = build_spec(PDAC_PROJECT)
    persistent = spec.design_constraints.normal_tissue_risk_ceiling
    terminable = spec.design_constraints.terminable_risk_ceiling
    tol = routing.Tolerances(persistent=persistent, terminable=terminable)
    print(f"  persistent ceiling {persistent}   terminable ceiling {terminable}")

    print("\nloading the ranked pool", flush=True)
    from car_pipeline.api import pipeline
    run = pipeline.run("Pancreatic Ductal Adenocarcinoma",
                       progress=lambda s, n="": None)
    pool = run["pool"]
    decisions = {d["gene"]: d for d in run["decisions"]}
    routed = {
        r.gene: routing.route(r.gene, r.risk, r.risk_organ, tol)
        for r in pool
    }
    print(f"  {len(pool)} pool members")

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    # A1 order independence
    shuffled = list(reversed(pool))
    again = {r.gene: routing.route(r.gene, r.risk, r.risk_organ, tol).architecture
             for r in shuffled}
    drift = [g for g in routed if routed[g].architecture != again[g]]
    criterion("A1", bool(drift),
              f"architecture is order-independent; {len(drift)} differ on a "
              "reversed pool")

    # A2 simplest admitting architecture wins
    wrong = [g for g, r in routed.items()
             if r.architecture == routing.ADAPTOR
             and r.risk is not None and r.risk <= persistent]
    criterion("A2", bool(wrong),
              f"{len(wrong)} targets routed ADAPTOR that CONVENTIONAL would "
              "have admitted")

    # A3 routing never rewrites a risk number
    by_gene = {r.gene: r.risk for r in pool}
    moved = [g for g, r in routed.items() if r.risk != by_gene[g]]
    criterion("A3", bool(moved),
              f"{len(moved)} routed risks differ from the Stage 3 risk")

    # A4 the persistent ceiling still binds
    over = [g for g, r in routed.items()
            if r.architecture == routing.CONVENTIONAL
            and (r.risk is None or r.risk > persistent)]
    criterion("A4", bool(over),
              f"{len(over)} CONVENTIONAL targets sit above the persistent "
              "ceiling")

    # A5 positive pin: NPSR1
    npsr1 = routed.get("NPSR1")
    criterion("A5", npsr1 is None or npsr1.architecture != routing.CONVENTIONAL,
              f"NPSR1 (risk {npsr1.risk:.4f}) routes "
              f"{npsr1.architecture}" if npsr1 else "NPSR1 not in the pool")

    # A6 positive pin: MSLN
    msln = routed.get("MSLN")
    criterion("A6", msln is None or msln.architecture != routing.ADAPTOR,
              f"MSLN (risk {msln.risk:.4f}, {msln.risk_organ}) routes "
              f"{msln.architecture}" if msln else "MSLN not in the pool")

    # A7 unbuilt rows are named, not silently dropped
    silent = [g for g, r in routed.items()
              if r.architecture == routing.NO_ARCHITECTURE and not r.reason]
    criterion("A7", bool(silent),
              f"{len(silent)} targets resolve to NO_ARCHITECTURE with no reason")

    # A8 the terminable ceiling is declared, not defaulted
    bare = routing.Tolerances(persistent=persistent)
    leak = [r.gene for r in pool
            if routing.route(r.gene, r.risk, r.risk_organ,
                             bare).architecture == routing.ADAPTOR]
    criterion("A8", bool(leak),
              f"with no declared terminable ceiling, {len(leak)} targets still "
              "route ADAPTOR")

    # A9 sensitivity is reported
    counts = routing.sweep(
        [(r.gene, r.risk, r.risk_organ) for r in pool], tol, SWEEP)
    criterion("A9", not counts,
              "adaptor admissions across the ceiling sweep: "
              + ", ".join(f"{k}->{v}" for k, v in counts.items()))

    # A10 not tuned to an outcome
    criterion("A10", terminable != SPEC_TERMINABLE_CEILING,
              f"declared ceiling {terminable} matches the spec value "
              f"{SPEC_TERMINABLE_CEILING}")

    # A11 an unsupplied binder is never reported as a finished sequence
    adaptors = [c for c in run["constructs"] if c.outcome == stage4.ADAPTOR]
    lying = [c for c in adaptors if c.amino_acid_sequence and not c.binder_supplied]
    criterion("A11", bool(lying),
              f"{len(adaptors)} adaptor constructs; {len(lying)} emit a "
              "sequence despite an unsupplied binder")

    print("=" * 72)
    print(f"  {11 - len(tripped)}/11 criteria clear")

    print()
    print("=" * 72)
    print("WHAT ROUTING CHANGED")
    print("=" * 72)
    from collections import Counter
    arch = Counter(r.architecture for r in routed.values())
    for name, n in arch.most_common():
        print(f"    {name:18s} {n:4d}")
    print()
    outcomes = Counter(d["outcome"] for d in decisions.values())
    print(f"    stage 4 outcomes: {dict(outcomes)}")
    for c in adaptors[:5]:
        print(f"    {c.gene:10s} {c.total_bp:5d} bp  headroom {c.headroom_bp:5d}"
              f"  {c.verdict}  sequence_supplied={c.binder_supplied}")

    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
