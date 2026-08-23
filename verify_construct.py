"""Runs the construct stage and tests it against the criteria fixed in the spec.

Criteria before biology. K1 and K2 are the two that would have caught Stage 5's
dead route: one checks the assembly did not corrupt the sequence, the other that
the binder arrived at all. Both are positive pins.
"""

from __future__ import annotations

import sys

from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.domains import PROTEOME, SYNTHETIC
from car_pipeline.stages import stage4, stage5, stage6

#: Pinned before the run. The construct for each of these must contain the named
#: therapeutic's variable regions verbatim. This is the Stage 5 to Stage 6 join,
#: and the equivalent join one stage earlier is the one that silently returned
#: nothing for all 200 targets.
#: Of the fourteen recommendations, only two duals have a binder on both arms.
#: They are the pins because they are the constructs the stage actually builds:
#: MSLN and CLDN18 both carry binders but Stage 4 returned NO_DESIGN for the
#: first and a partner without a binder for the second, so neither yields a
#: construct and neither can test the join.
PINNED_ASSEMBLED = ["MUC16", "MUC17"]


def main() -> int:
    print("loading upstream", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    source = AntibodySource()
    records = stage5.retrieve(decisions, source, progress=False)
    binders = {r.gene: r for r in records}
    stage5_hash = stage5.configuration_hash(
        manifest["stage4_hash"], [r.gene for r in records])
    constructs = stage6.build(decisions, binders)
    by_gene = {c.gene: c for c in constructs}
    built = [c for c in constructs if c.amino_acid_sequence]
    print(f"  {len(decisions)} decisions, {len(built)} assembled")

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    bad_round_trip = [c.gene for c in built
                      if stage6.translate(c.dna) != c.amino_acid_sequence]
    criterion("K1", bool(bad_round_trip),
              f"{len(built)} constructs translate back to their own sequence"
              if not bad_round_trip else f"round trip fails for {bad_round_trip[:5]}")

    k2_bad = []
    # A stage that assembled nothing would pass a check phrased only over what it
    # assembled. The pins make that impossible.
    for gene in PINNED_ASSEMBLED:
        c = by_gene.get(gene)
        if c is None or not c.amino_acid_sequence:
            k2_bad.append(f"{gene}: expected a construct, got none")
    if not built:
        k2_bad.append("nothing was assembled at all")
    for c in built:
        record = binders.get(c.gene)
        picked = next((t for t in record.sequence if t.name == c.binder_name), None)
        if picked is None:
            k2_bad.append(f"{c.gene}: chosen binder {c.binder_name} not in Stage 5")
            continue
        if picked.heavy_sequence not in c.amino_acid_sequence:
            k2_bad.append(f"{c.gene}: VH of {c.binder_name} absent")
        if picked.light_sequence not in c.amino_acid_sequence:
            k2_bad.append(f"{c.gene}: VL of {c.binder_name} absent")
        if c.partner_binder_name:
            pr = binders.get(c.partner)
            pp = next((t for t in pr.sequence if t.name == c.partner_binder_name),
                      None)
            if pp is None or pp.heavy_sequence not in c.amino_acid_sequence:
                k2_bad.append(f"{c.gene}: partner VH absent")
    criterion("K2", bool(k2_bad),
              f"{len(built)} constructs carry their binders verbatim, including "
              f"the pinned {', '.join(PINNED_ASSEMBLED)}"
              if not k2_bad else "; ".join(k2_bad[:4]))

    k3_bad = []
    for c in built:
        pos = 0
        for seg in c.segments:
            if seg.aa_start != pos or seg.aa_end > len(c.amino_acid_sequence):
                k3_bad.append(c.gene)
                break
            pos = seg.aa_end
        else:
            if pos != len(c.amino_acid_sequence):
                k3_bad.append(c.gene)
    criterion("K3", bool(k3_bad),
              "domain boundaries partition every construct exactly"
              if not k3_bad else f"gaps or overlaps in {k3_bad[:5]}")

    undescribed = [
        (c.gene, s.name) for c in built for s in c.segments
        if not (s.provenance == SYNTHETIC
                or (s.provenance == PROTEOME and s.accession and s.start_residue)
                or s.provenance == "stage5")
    ]
    criterion("K4", bool(undescribed),
              f"every part of every construct names its source "
              f"({len(built[0].segments) if built else 0} parts in the first construct)"
              if not undescribed else f"{len(undescribed)} parts without a source")

    k5_bad = []
    for c in built:
        # Both terms recomputed from the segments rather than from the construct's
        # own properties. Comparing `headroom_bp` against its own definition is a
        # tautology and cannot fail.
        summed = sum(s.residues for s in c.segments) * 3 + len(stage6.STOP)
        expected_headroom = stage6.BUDGET_BP - summed
        if summed != c.total_bp or expected_headroom != c.headroom_bp:
            k5_bad.append(c.gene)
    criterion("K5", bool(k5_bad),
              "part costs sum to the printed total for every construct"
              if not k5_bad else f"arithmetic disagrees for {k5_bad[:5]}")

    no_switch = [c.gene for c in built
                 if c.verdict == stage6.BUILDABLE and not c.has_switch]
    criterion("K6", bool(no_switch),
              "every buildable construct carries the mandatory safety switch"
              if not no_switch else f"{no_switch[:5]} buildable without a switch")

    # Per spec 5.1: a binder is necessary and not sufficient. A construct is owed
    # only where the target also carries a recommendation, and for a dual, where
    # the partner has a binder too.
    def usable(gene):
        r = binders.get(gene)
        return bool(r and [t for t in r.sequence
                           if t.heavy_sequence and t.light_sequence
                           and not stage6.assemblable(
                               t.heavy_sequence + t.light_sequence)])

    k7_bad, withheld_by_outcome = [], 0
    for c in constructs:
        owed = (
            usable(c.gene)
            and c.outcome in stage6.BUILDABLE_OUTCOMES
            and (c.outcome != "DUAL" or usable(c.partner))
        )
        if owed and not c.amino_acid_sequence:
            k7_bad.append(f"{c.gene}: owed a construct, got none")
        if not usable(c.gene) and c.amino_acid_sequence:
            k7_bad.append(f"{c.gene}: construct without a usable binder")
        if usable(c.gene) and c.outcome not in stage6.BUILDABLE_OUTCOMES:
            withheld_by_outcome += 1
    criterion("K7", bool(k7_bad),
              f"every owed construct was built and none was built without a binder; "
              f"{withheld_by_outcome} targets have a binder but no recommendation "
              f"(§5.1)" if not k7_bad else "; ".join(k7_bad[:3]))

    in_genes = {d["gene"] for d in decisions}
    out_genes = {c.gene for c in constructs}
    criterion("K8", len(constructs) != len(decisions) or in_genes != out_genes,
              f"{len(constructs)} rows out of {len(decisions)}; "
              f"{len(in_genes - out_genes)} dropped, {len(out_genes - in_genes)} added")

    print("=" * 72)
    print(f"  {8 - len(tripped)}/8 criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    # ---------------- the biology ---------------------------------------
    print()
    print("=" * 72)
    print("WHAT CAN BE BUILT")
    print("=" * 72)
    counts: dict[str, int] = {}
    for c in constructs:
        counts[c.verdict] = counts.get(c.verdict, 0) + 1
    for verdict in (stage6.BUILDABLE, stage6.BUDGET_EXCEEDED, stage6.NO_CONSTRUCT):
        n = counts.get(verdict, 0)
        print(f"    {verdict:18s} {n:4d}  ({n / len(constructs):.0%})")

    singles = [c for c in built if c.outcome != "DUAL"]
    duals = [c for c in built if c.outcome == "DUAL"]
    print()
    print(f"    single-antigen assembled {len(singles)}, "
          f"{sum(1 for c in singles if c.verdict == stage6.BUILDABLE)} within budget")
    print(f"    dual assembled           {len(duals)}, "
          f"{sum(1 for c in duals if c.verdict == stage6.BUILDABLE)} within budget")

    print()
    print("  The budget, itemised")
    print("  " + "-" * 68)
    example = next((c for c in built if c.outcome == "SINGLE"),
                   built[0] if built else None)
    if example:
        print(f"    {example.gene} — {example.architecture}, binder "
              f"{example.binder_name}")
        for s in example.segments:
            where = (f"{s.accession} {s.start_residue}-{s.end_residue}"
                     if s.provenance == PROTEOME else s.provenance)
            print(f"      {s.name:24s} {s.residues:4d} aa  {s.bp_end - s.bp_start:5d} bp"
                  f"   {where}")
        print(f"      {'stop':24s} {'':4s}     {len(stage6.STOP):5d} bp")
        print(f"      {'TOTAL':24s} {'':4s}     {example.total_bp:5d} bp"
              f"   headroom {example.headroom_bp} of {stage6.BUDGET_BP}")

    if duals:
        d = duals[0]
        print()
        print(f"  And one dual: {d.gene} + {d.partner} — {d.total_bp} bp, "
              f"{d.verdict}")
        if d.verdict == stage6.BUDGET_EXCEEDED:
            print(f"      {d.reason}")
            print("      This is the expected result. Two single-chain binders plus")
            print("      the mandatory switch do not fit, and Stage 5 retrieved one")
            print("      single-domain binder in 720 candidates, so the smaller")
            print("      format that would fit has no inventory to draw on.")

    print()
    print(f"  configuration hash "
          f"{stage6.configuration_hash(stage5_hash, [c.gene for c in constructs])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
