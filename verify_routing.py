"""Stage 4a against the criteria fixed in specs/stage4a-architecture-routing.md."""

from __future__ import annotations

import os
import sys

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.stages import routing, stage4, stage9
from car_pipeline.stages.stage1 import build_spec


SPEC_TERMINABLE_CEILING = 0.35


SWEEP = (0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70)


def main() -> int:
    """Run the architecture-routing criteria."""
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

    checked: list[str] = []
    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    shuffled = list(reversed(pool))
    again = {r.gene: routing.route(r.gene, r.risk, r.risk_organ, tol).architecture
             for r in shuffled}
    drift = [g for g in routed if routed[g].architecture != again[g]]
    criterion("A1", bool(drift),
              f"architecture is order-independent; {len(drift)} differ on a "
              "reversed pool")

    wrong = [g for g, r in routed.items()
             if r.architecture == routing.ADAPTOR
             and r.risk is not None and r.risk <= persistent]
    criterion("A2", bool(wrong),
              f"{len(wrong)} targets routed ADAPTOR that CONVENTIONAL would "
              "have admitted")

    by_gene = {r.gene: r.risk for r in pool}
    moved = [g for g, r in routed.items() if r.risk != by_gene[g]]
    criterion("A3", bool(moved),
              f"{len(moved)} routed risks differ from the Stage 3 risk")

    over = [g for g, r in routed.items()
            if r.architecture == routing.CONVENTIONAL
            and (r.risk is None or r.risk > persistent)]
    criterion("A4", bool(over),
              f"{len(over)} CONVENTIONAL targets sit above the persistent "
              "ceiling")

    npsr1 = routed.get("NPSR1")
    criterion("A5", npsr1 is None or npsr1.architecture != routing.CONVENTIONAL,
              f"NPSR1 (risk {npsr1.risk:.4f}) routes "
              f"{npsr1.architecture}" if npsr1 else "NPSR1 not in the pool")

    msln = routed.get("MSLN")
    criterion("A6", msln is None or msln.architecture != routing.ADAPTOR,
              f"MSLN (risk {msln.risk:.4f}, {msln.risk_organ}) routes "
              f"{msln.architecture}" if msln else "MSLN not in the pool")

    silent = [g for g, r in routed.items()
              if r.architecture == routing.NO_ARCHITECTURE and not r.reason]
    criterion("A7", bool(silent),
              f"{len(silent)} targets resolve to NO_ARCHITECTURE with no reason")

    bare = routing.Tolerances(persistent=persistent)
    leak = [r.gene for r in pool
            if routing.route(r.gene, r.risk, r.risk_organ,
                             bare).architecture == routing.ADAPTOR]
    criterion("A8", bool(leak),
              f"with no declared terminable ceiling, {len(leak)} targets still "
              "route ADAPTOR")

    counts = routing.sweep(
        [(r.gene, r.risk, r.risk_organ) for r in pool], tol, SWEEP)
    criterion("A9", not counts,
              "adaptor admissions across the ceiling sweep: "
              + ", ".join(f"{k}->{v}" for k, v in counts.items()))

    criterion("A10", terminable != SPEC_TERMINABLE_CEILING,
              f"declared ceiling {terminable} matches the spec value "
              f"{SPEC_TERMINABLE_CEILING}")

    adaptors = [c for c in run["constructs"] if c.outcome == stage4.ADAPTOR]
    lying = [c for c in adaptors if c.amino_acid_sequence and not c.binder_supplied]
    criterion("A11", bool(lying),
              f"{len(adaptors)} adaptor constructs; {len(lying)} emit a "
              "sequence despite an unsupplied binder")

    gated = {g.gene: g for g in run["gated"]}
    blind = []
    connected = []
    for c in adaptors:
        record = gated.get(c.gene)
        if record is None:
            blind.append(f"{c.gene}=no safety record")
            continue
        if (record.binder_origin == "human"
                or not record.binder_source_organism
                or not record.binder_structure_accession):
            blind.append(
                f"{c.gene}={record.binder_origin or 'unset'}/"
                f"{record.binder_source_organism or 'no organism'}")
        if record.epitope_immunogenicity != stage9.NOT_CONNECTED:
            connected.append(c.gene)
    criterion(
        "A12", bool(blind) or bool(connected),
        (f"{len(adaptors) - len(blind)} of {len(adaptors)} adaptor constructs "
         f"carry a structure-derived binder the origin check can see"
         + (f"; blind on {blind[:3]}" if blind else "")
         + (f"; epitope immunogenicity unexpectedly connected on {connected[:3]}"
            if connected else
            "; epitope immunogenicity stays NOT_CONNECTED on all of them, so "
            "the species gap and the immunogenicity gap remain separate"))
    )

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")

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
