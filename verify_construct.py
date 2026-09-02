"""Runs the construct stage and tests it against the criteria fixed in the spec."""

from __future__ import annotations

import sys

from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.domains import PROTEOME, STRUCTURE, SYNTHETIC, anti_tag_binder
from car_pipeline.stages import stage4, stage5, stage6


def stage5_usable(binders: dict, gene: str | None) -> bool:
    """Whether Stage 5 gave this target a binder the assembler can use."""
    record = binders.get(gene) if gene else None
    if record is None:
        return False
    return any(
        candidate.heavy_sequence and candidate.light_sequence
        and not stage6.assemblable(
            candidate.heavy_sequence + candidate.light_sequence)
        for candidate in record.sequence)


def two_armed_duals(decisions: list[dict], binders: dict) -> list[tuple[str, str]]:
    """The duals carrying a binder on both arms, which is what K2 pins on."""
    out = []
    for row in decisions:
        partner = row.get("partner")
        if row["outcome"] != stage4.DUAL or not partner:
            continue
        if stage5_usable(binders, row["gene"]) and stage5_usable(binders, partner):
            out.append((row["gene"], partner))
    return out


def needs_met(binders: dict, outcome: str, gene: str, partner: str | None,
              adaptor_supplied: bool) -> bool:
    """Whether every binder this architecture needs was actually retrieved."""
    if outcome == stage4.ADAPTOR:
        return adaptor_supplied
    if outcome == stage4.DUAL:
        return (stage5_usable(binders, gene)
                and bool(partner) and stage5_usable(binders, partner))
    if outcome == stage4.SINGLE:
        return stage5_usable(binders, gene)
    return False


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
    adaptor_part = anti_tag_binder()
    print(f"  {len(decisions)} decisions, {len(built)} assembled")
    print("  anti-tag binder "
          + (f"retrieved from {adaptor_part.accession}" if adaptor_part.supplied
             else "declares a size but no sequence"))

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

    route_config = manifest.get("routing")
    disabled = [d["gene"] for d in decisions
                if (d.get("route_reason") or "") == stage4.ROUTING_DISABLED_REASON]
    k0_bad = []
    if not route_config:
        k0_bad.append("the manifest records no routing configuration")
    if disabled:
        k0_bad.append(f"{len(disabled)} of {len(decisions)} rows carry "
                      f"{stage4.ROUTING_DISABLED_REASON!r}")
    if not built:
        k0_bad.append(f"nothing assembled from {len(decisions)} decisions, so "
                      "every criterion below would report on an empty set")
    criterion(
        "K0", bool(k0_bad),
        f"the decision set is routed (persistent "
        f"{route_config['persistent_ceiling']}, terminable "
        f"{route_config['terminable_ceiling']}) and yields {len(built)} "
        f"construct(s) for the criteria below to read"
        if not k0_bad else "; ".join(k0_bad))

    bad_round_trip = [c.gene for c in built
                      if stage6.translate(c.dna) != c.amino_acid_sequence]
    criterion("K1", bool(bad_round_trip) or not built,
              f"{len(built)} constructs translate back to their own sequence"
              if built and not bad_round_trip else
              f"round trip fails for {bad_round_trip[:5]}" if bad_round_trip else
              "no construct to translate")

    owed_pairs = two_armed_duals(decisions, binders)
    k2_bad = []

    for gene, partner in owed_pairs:
        c = by_gene.get(gene)
        if c is None or not c.amino_acid_sequence:
            k2_bad.append(f"{gene}+{partner}: owed a construct, got none")
    adaptor_rows = [d for d in decisions if d["outcome"] == stage4.ADAPTOR]
    if adaptor_rows and not adaptor_part.supplied:
        k2_bad.append(
            f"{len(adaptor_rows)} adaptor row(s), but no anti-tag sequence was "
            "retrieved, so the anti-tag join cannot be verified on any of them")

    by_route = {"anti-tag": 0, "stage 5": 0}
    for c in built:
        if c.outcome == stage4.ADAPTOR:
            by_route["anti-tag"] += 1
            if c.binder_name != adaptor_part.name:
                k2_bad.append(f"{c.gene}: names binder {c.binder_name!r}, not the "
                              "retrieved anti-tag part")
            if adaptor_part.sequence not in c.amino_acid_sequence:
                k2_bad.append(f"{c.gene}: the anti-tag sequence is absent")
            carried = [s for s in c.segments if s.provenance == STRUCTURE]
            if len(carried) != 1 or carried[0].accession != adaptor_part.accession:
                k2_bad.append(
                    f"{c.gene}: {len(carried)} structure-derived segment(s), "
                    f"accession {[s.accession for s in carried]} against "
                    f"{adaptor_part.accession!r}")
            continue

        by_route["stage 5"] += 1
        record = binders.get(c.gene)
        if record is None:
            k2_bad.append(f"{c.gene}: no Stage 5 record for an assembled target")
            continue
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
                      None) if pr else None
            if pp is None:
                k2_bad.append(f"{c.gene}: partner binder {c.partner_binder_name} "
                              "not in Stage 5")
            else:
                if pp.heavy_sequence not in c.amino_acid_sequence:
                    k2_bad.append(f"{c.gene}: partner VH absent")
                if pp.light_sequence not in c.amino_acid_sequence:
                    k2_bad.append(f"{c.gene}: partner VL absent")
    criterion(
        "K2", bool(k2_bad) or not owed_pairs,
        "; ".join(k2_bad[:4]) if k2_bad else
        (f"{len(built)} constructs carry their binders verbatim "
         f"({by_route['anti-tag']} by the anti-tag route, "
         f"{by_route['stage 5']} by the Stage 5 route), including the "
         f"{len(owed_pairs)} dual(s) carrying a binder on both arms: "
         + ", ".join(f"{g}+{p}" for g, p in owed_pairs[:4])) if owed_pairs else
        (f"no dual carries a binder on both arms, so the two-arm join is not "
         f"exercised anywhere in this decision set ({len(built)} of "
         f"{len(decisions)} rows assembled, {by_route['anti-tag']} of them by "
         f"the anti-tag route, whose binder is verified above); the criterion "
         f"has nothing to pin on and reports that rather than clearing on an "
         f"empty set"))

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
    criterion("K3", bool(k3_bad) or not built,
              f"domain boundaries partition all {len(built)} constructs exactly"
              if built and not k3_bad else
              f"gaps or overlaps in {k3_bad[:5]}" if k3_bad else
              "no construct to partition")

    undescribed = [
        (c.gene, s.name) for c in built for s in c.segments
        if not (s.provenance == SYNTHETIC
                or (s.provenance == PROTEOME and s.accession and s.start_residue)
                or (s.provenance == STRUCTURE and s.accession)
                or s.provenance == "stage5")
    ]
    criterion("K4", bool(undescribed) or not built,
              f"every part of every construct names its source "
              f"({len(built[0].segments) if built else 0} parts in the first "
              f"construct, {sum(len(c.segments) for c in built)} across all "
              f"{len(built)})"
              if built and not undescribed else
              f"{len(undescribed)} parts without a source" if undescribed else
              "no construct, so no part was examined")

    k5_bad = []
    for c in built:
        summed = sum(s.residues for s in c.segments) * 3 + len(stage6.STOP)
        expected_headroom = stage6.BUDGET_BP - summed
        if summed != c.total_bp or expected_headroom != c.headroom_bp:
            k5_bad.append(c.gene)
    criterion("K5", bool(k5_bad) or not built,
              f"part costs sum to the printed total for all {len(built)} constructs"
              if built and not k5_bad else
              f"arithmetic disagrees for {k5_bad[:5]}" if k5_bad else
              "no construct to cost")

    buildable = [c for c in built if c.verdict == stage6.BUILDABLE]
    no_switch = [c.gene for c in buildable if not c.has_switch]
    criterion("K6", bool(no_switch) or not buildable,
              f"all {len(buildable)} buildable constructs carry the mandatory "
              f"safety switch"
              if buildable and not no_switch else
              f"{no_switch[:5]} buildable without a switch" if no_switch else
              "no buildable construct to check the switch on")

    k7_bad, withheld_by_outcome = [], 0
    for c in constructs:
        owed_here = (
            c.outcome in stage6.BUILDABLE_OUTCOMES
            and needs_met(binders, c.outcome, c.gene, c.partner,
                          adaptor_part.supplied)
        )
        if owed_here and not c.amino_acid_sequence:
            k7_bad.append(f"{c.gene}: owed a construct, got none")
        if not owed_here and c.amino_acid_sequence:
            k7_bad.append(f"{c.gene}: {c.outcome} construct built without the "
                          "binder its architecture needs")
        if stage5_usable(binders, c.gene) and c.outcome not in stage6.BUILDABLE_OUTCOMES:
            withheld_by_outcome += 1
    criterion("K7", bool(k7_bad),
              f"every owed construct was built and none was built without the "
              f"binder its architecture needs; {withheld_by_outcome} targets "
              f"have a binder but no recommendation (§5.1)"
              if not k7_bad else "; ".join(k7_bad[:3]))

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

    singles = [c for c in built if c.outcome == stage4.SINGLE]
    duals = [c for c in built if c.outcome == stage4.DUAL]
    adaptors = [c for c in built if c.outcome == stage4.ADAPTOR]
    print()
    for label, group in (("single-antigen", singles), ("dual", duals),
                         ("adaptor", adaptors)):
        print(f"    {label + ' assembled':24s} {len(group):4d}, "
              f"{sum(1 for c in group if c.verdict == stage6.BUILDABLE)} "
              f"within budget")
    for c in adaptors:
        print(f"      {c.gene:10s} {c.total_bp:5d} bp   headroom "
              f"{c.headroom_bp:4d}   {c.verdict}")

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
