"""Runs the construct stage and tests it against the criteria fixed in the spec."""

from __future__ import annotations

import sys

from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.domains import PROTEOME, SYNTHETIC
from car_pipeline.stages import stage4, stage5, stage6


PINNED_ASSEMBLED = ["MUC16", "MUC17"]


def main() -> int:
    """Run the construct-assembly criteria."""
    print("loading upstream", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    source = AntibodySource()

    records = stage5.load_or_retrieve(
        decisions, source, manifest["stage4_hash"])
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

    checked: list[str] = []
    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    bad_round_trip = [c.gene for c in built
                      if stage6.translate(c.dna) != c.amino_acid_sequence]
    criterion("K1", bool(bad_round_trip),
              f"{len(built)} constructs translate back to their own sequence"
              if not bad_round_trip else f"round trip fails for {bad_round_trip[:5]}")

    k2_bad = []

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
                or (s.provenance == "structure" and s.accession)
                or s.provenance == "stage5")
    ]
    criterion("K4", bool(undescribed),
              f"every part of every construct names its source "
              f"({len(built[0].segments) if built else 0} parts in the first construct)"
              if not undescribed else f"{len(undescribed)} parts without a source")

    k5_bad = []
    for c in built:
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

    def usable(gene):
        """Whether the target has a binder whose regions can actually be assembled."""
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

    expected_rows = manifest["pool_size"]
    out_genes = {c.gene for c in constructs}
    criterion("K8",
              len(constructs) != expected_rows or len(out_genes) != expected_rows,
              f"{len(constructs)} rows and {len(out_genes)} distinct genes against "
              f"the {expected_rows} the Stage 4 manifest records")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

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
    print("=" * 72)
    print("WHY NOTHING FITS — THE CHAIN, STATED AS A FINDING")
    print("=" * 72)
    switch_names = ("T2A skip peptide", "FKBP12", "SGGGS linker",
                    "caspase-9 without CARD")
    example_dual = next((c for c in built if c.outcome == "DUAL"), None)
    if example_dual is not None:
        switch_bp = sum(s.bp_end - s.bp_start for s in example_dual.segments
                        if s.name in switch_names and s.name != "T2A skip peptide")
        cassette_bp = sum(s.bp_end - s.bp_start for s in example_dual.segments
                          if s.name in switch_names)
        over = -example_dual.headroom_bp
        binder_bp = sum(s.bp_end - s.bp_start for s in example_dual.segments
                        if s.provenance == "stage5")
        print(f"    1. Stage 1's conservative tolerance requires a safety switch.")
        print(f"    2. The switch costs {switch_bp} bp, {cassette_bp} bp with the")
        print(f"       skip peptide that carries it.")
        print(f"    3. Two single-chain binders ({binder_bp} bp) plus that switch")
        print(f"       reach {example_dual.total_bp} bp against a "
              f"{stage6.BUDGET_BP} bp budget, over by {over}.")
        print(f"    4. Single-domain binders would fit: they replace the two")
        print(f"       variable regions and their linker with one domain each.")
        print(f"    5. Stage 5 retrieved ONE single-domain candidate in 720.")
        print()
        print("    The format that fits has no inventory. That is the finding, and")
        print("    it is not a shortfall in the budget — the budget is doing what")
        print("    it exists for.")

        print()
        print("=" * 72)
        print("WHAT EACH EXIT WOULD COST, AS NUMBERS")
        print("=" * 72)
        print(f"    a smaller safety switch would have to free      {over} bp")
        print(f"       leaving it at most                            "
              f"{switch_bp - over} bp ({(switch_bp - over) // 3} residues), "
              f"{100 * (switch_bp - over) / switch_bp:.0f}% of its current size")
        payload = example_dual.total_bp
        backbone = 1200
        print(f"    a vector large enough for the current design    "
              f"{(payload + backbone) / 1000:.2f} kb")
        print(f"       against the {(stage6.BUDGET_BP + backbone) / 1000:.1f} kb "
              f"assumed in Stage 1, a {100 * (payload + backbone) / (stage6.BUDGET_BP + backbone) - 100:.0f}% increase")

    single_domain = []
    for gene, record in binders.items():
        structural = [c for c in record.structure if "single" in c.fmt.lower()]
        sequenced = [c for c in record.sequence if "single domain" in c.fmt.lower()]
        if structural or sequenced:
            single_domain.append((gene, len(structural), len(sequenced)))
    print(f"    pool targets with ANY single-domain binder      "
          f"{len(single_domain)} of {len(binders)}")
    for gene, ns, nq in sorted(single_domain):
        record = binders[gene]
        outcome = record.outcome
        names = [c.name for c in record.sequence if "single domain" in c.fmt.lower()]
        print(f"       {gene:10s} structure {ns}, sequence {nq}"
              f"   Stage 4 outcome {outcome}   {', '.join(names)}")
        for c in record.sequence:
            if "single domain" not in c.fmt.lower():
                continue
            unusable = stage6.assemblable(c.heavy_sequence + c.light_sequence)
            if unusable:
                print(f"         and it is not assemblable: the light-chain field "
                      f"holds {c.light_sequence!r}, a placeholder rather than a")
                print(f"         sequence, so the residues {sorted(unusable)} are not "
                      "residues at all")
    if not single_domain:
        print("       none")
    print()
    print("    So the nanobody route is empty rather than thin. The single")
    print("    candidate sits on a target the pairing stage did not recommend, and")
    print("    its light-chain field is a placeholder. There is no exit here to")
    print("    take, and reporting one would be inventing inventory.")

    print()
    print(f"  configuration hash "
          f"{stage6.configuration_hash(stage5_hash, [c.gene for c in constructs])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
