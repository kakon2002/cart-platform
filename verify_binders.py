"""Runs the binder stage and tests it against the criteria fixed in the spec.

Order is deliberate: the known-answer checks first, then the criteria, then the
biology. A retrieval that silently returns nothing looks exactly like a target
with no binders, and only a target whose answer is known separates the two — so
nothing else is worth reading until those have passed.
"""

from __future__ import annotations

import sys

from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.stages import stage4, stage5

#: The five known targets for this indication, taken from the ranking stage's own
#: list rather than restated, so the two cannot drift apart.
from verify_ranking import KNOWN_TARGETS

#: Pinned before the run. Every entry verified present by accession-anchored
#: search while the specification was written; every negative verified absent.
#: Without the negatives the check is one-sided and a stage that returned
#: everything would pass it.
EXPECTED_SEQUENCE = {
    "MSLN": {"Amatuximab", "Anetumab"},
    "CLDN18": {"Zolbetuximab"},
    "CEACAM5": {"Tusamitamab"},
    "CEACAM6": {"Tinurilimab"},
    "MUC1": {"Cantuzumab"},
}
#: Entries that MUST appear on the structure route. Every one verified present by
#: accession-anchored search and confirmed antibody-containing in the curated
#: summary. This half of the check is not optional: the first run of this stage
#: returned zero structure candidates for all 200 targets because the two sources
#: identify entries differently, and every structure-route check in place at the
#: time was a negative, so all of them passed while the route was dead.
EXPECTED_STRUCTURE = {
    "MSLN": {"4f3f", "7ued", "8cxc", "8cz8"},
    "CLDN18": {"9v32"},
}

#: Targets that must return no structure-route candidate. CEACAM6 has entries but
#: none containing an antibody; NRG3 has no entries at all. A stage that emitted
#: a candidate for either is not filtering.
EXPECTED_NO_STRUCTURE = ["CEACAM6", "NRG3"]

#: CEACAM5 has six accession-anchored entries and exactly one containing an
#: antibody. It is the sharpest single check in the file: a stage that echoed its
#: entry count would return six, and one that had lost the join would return
#: none, so only the correct answer sits at one.
CEACAM5_STRUCTURE_ENTRIES = 1


def main() -> int:
    print("loading stage 4 decisions", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    print(f"  {len(decisions)} decisions, stage4 hash {manifest['stage4_hash']}, "
          f"usable as a result: {'yes' if manifest['usable_as_result'] else 'no'}")
    if not manifest["usable_as_result"]:
        print(f"  upstream tripped {', '.join(manifest['criteria_tripped'])} — this "
              "stage annotates them and does not treat them as settled")

    source = AntibodySource()
    print("retrieving binders", flush=True)
    records = stage5.retrieve(decisions, source)
    by_gene = {r.gene: r for r in records}

    # ---------------- known answers ------------------------------------
    print()
    print("=" * 72)
    print("KNOWN-ANSWER CHECKS")
    print("=" * 72)
    failures: list[str] = []

    for gene, expected in EXPECTED_SEQUENCE.items():
        record = by_gene.get(gene)
        if record is None:
            print(f"  FAIL  {gene:9s} absent from the pool")
            failures.append(gene)
            continue
        found = {c.name for c in record.sequence}
        missing = expected - found
        ok = not missing
        if not ok:
            failures.append(gene)
        print(f"  {'ok  ' if ok else 'FAIL'}  {gene:9s} "
              f"{len(record.sequence)} sequence binder(s); expected "
              f"{sorted(expected)} {'' if ok else '-> MISSING ' + str(sorted(missing))}")

    print()
    for gene, expected in EXPECTED_STRUCTURE.items():
        record = by_gene.get(gene)
        if record is None:
            print(f"  FAIL  {gene:9s} absent from the pool")
            failures.append(gene)
            continue
        found = {c.identifier.split(":")[0].lower() for c in record.structure}
        missing = expected - found
        ok = not missing
        if not ok:
            failures.append(gene)
        print(f"  {'ok  ' if ok else 'FAIL'}  {gene:9s} "
              f"{len(record.structure)} structure candidate(s) over "
              f"{len(found)} entries; expected {sorted(expected)}"
              f"{'' if ok else ' -> MISSING ' + str(sorted(missing))}")

    print()
    for gene in EXPECTED_NO_STRUCTURE:
        record = by_gene.get(gene)
        if record is None:
            print(f"  ----  {gene:9s} not in the pool, cannot be checked")
            continue
        ok = not record.structure
        if not ok:
            failures.append(gene)
        print(f"  {'ok  ' if ok else 'FAIL'}  {gene:9s} "
              f"{len(record.entries)} entries, {len(record.structure)} structure "
              f"candidate(s); expected 0")

    if failures:
        print(f"\n  STOP: the known answers fail for {', '.join(sorted(set(failures)))}. "
              "Nothing below this line is worth reading.")
        return 1

    # ---------------- criteria -----------------------------------------
    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    # B1' — the filters must subtract. A stage that echoed the entry count would
    # return a candidate for every entry and never for none.
    partial = [r for r in records if r.entries and not r.structure]
    echoed = [r for r in records
              if r.entries and len(r.structure) == len(r.entries)]
    # Both halves gate. A stage that returned a candidate per entry would show
    # zero partials and every target echoing; one that had lost the join would
    # show every target partial and none echoing. The correct answer is a mix.
    with_entries = [r for r in records if r.entries]
    criterion("B1",
              not partial
              or len(echoed) == len(with_entries)
              or len(partial) == len(with_entries),
              f"{len(partial)} targets have entries but no antibody among them "
              f"(an inert stage would have none); {len(echoed)} of "
              f"{len([r for r in records if r.entries])} with entries echo their count")

    # Recomputed rather than reading `failures`, which the early return above
    # guarantees is empty by the time this runs — a criterion that cannot fail is
    # the thing this project keeps deleting.
    ceacam5 = by_gene.get("CEACAM5")
    ceacam5_entries = len(ceacam5.entries) if ceacam5 else 0
    ceacam5_candidates = len({c.identifier.split(":")[0].lower()
                              for c in ceacam5.structure} if ceacam5 else set())
    b3_bad = (
        ceacam5 is None
        or ceacam5_candidates != CEACAM5_STRUCTURE_ENTRIES
        or any(
            not (EXPECTED_STRUCTURE[g] <= {
                c.identifier.split(":")[0].lower()
                for c in by_gene[g].structure})
            for g in EXPECTED_STRUCTURE if g in by_gene
        )
        or any(
            not (EXPECTED_SEQUENCE[g] <= {c.name for c in by_gene[g].sequence})
            for g in EXPECTED_SEQUENCE if g in by_gene
        )
        or any(by_gene[g].structure for g in EXPECTED_NO_STRUCTURE if g in by_gene)
    )
    criterion("B3", b3_bad,
              f"known answers hold: {len(EXPECTED_SEQUENCE)} sequence targets, "
              f"{len(EXPECTED_STRUCTURE)} pinned structure targets, "
              f"{len(EXPECTED_NO_STRUCTURE)} documented negatives, and CEACAM5 at "
              f"{ceacam5_candidates} structure entries of {ceacam5_entries} "
              f"(expected {CEACAM5_STRUCTURE_ENTRIES})")

    with_any = [r for r in records if r.verdict != stage5.NO_BINDER]
    knowns = [g for g in KNOWN_TARGETS if by_gene.get(g)
              and by_gene[g].verdict != stage5.NO_BINDER]
    criterion("B8", len(knowns) < len(KNOWN_TARGETS),
              f"{len(knowns)} of {len(KNOWN_TARGETS)} known targets return a binder "
              f"on some route")

    numeric_affinity = [c for r in records for c in (r.structure + r.sequence)
                        if c.affinity != stage5.NOT_CONNECTED]
    criterion("B5", bool(numeric_affinity),
              f"{len(numeric_affinity)} candidates carry an affinity value "
              "(the source does not have one)")

    ordered = [r.pool_index for r in records]
    criterion("B7", ordered != sorted(ordered),
              "pool order carried through unchanged" if ordered == sorted(ordered)
              else "the stage reordered its input")

    in_genes = {d["gene"] for d in decisions}
    out_genes = {r.gene for r in records}
    missing = in_genes - out_genes
    extra = out_genes - in_genes
    criterion("B13", bool(missing or extra) or len(records) != len(decisions),
              f"{len(records)} records out of {len(decisions)} decisions in; "
              f"{len(missing)} genes dropped, {len(extra)} added")

    unresolved = [c for r in records for c in (r.structure + r.sequence)
                  if c.isoform != stage5.ISOFORM_UNRESOLVED]
    criterion("B11", bool(unresolved),
              f"{len(unresolved)} candidates claim a resolved isoform "
              "(neither route can determine one)")

    print("=" * 72)
    print(f"  {7 - len(tripped)}/7 criteria clear")

    # ---------------- the biology ---------------------------------------
    print()
    print("=" * 72)
    print("WHAT THE LITERATURE HOLDS FOR THIS POOL")
    print("=" * 72)
    counts: dict[str, int] = {}
    for r in records:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    for verdict in (stage5.STRUCTURE_AND_SEQUENCE, stage5.BINDER_STRUCTURE_ONLY,
                    stage5.BINDER_SEQUENCE_ONLY, stage5.NO_BINDER):
        n = counts.get(verdict, 0)
        print(f"    {verdict:24s} {n:4d}  ({n / len(records):.0%})")
    print()
    print(f"    targets with entries but no antibody in any of them: "
          f"{len(partial)}")
    print("    That gap is the point: a protein having been crystallised is not")
    print("    the same as a binder existing against it, and only the curated")
    print("    chain annotation separates the two.")

    print()
    print("  The two validation targets in full")
    print("  " + "-" * 68)
    for gene in ("MSLN", "CLDN18"):
        r = by_gene[gene]
        print(f"    {gene} ({r.accession}) — {r.verdict}, Stage 4 outcome {r.outcome}")
        print(f"      {len(r.entries)} entries, {len(r.structure)} antibody complexes")
        for c in r.sequence:
            # Built as one string rather than a conditional expression over a
            # concatenation: written the latter way the condition binds to the
            # whole thing, so a therapeutic with only one variable region printed
            # a blank line and lost its name, stage and status with it.
            size = f"  scFv {c.car_bp} bp" if c.car_bp else "  size NOT_COMPUTABLE"
            print(f"      {c.name:16s} {c.clinical_stage:12s} {c.status:14s} "
                  f"VH {len(c.heavy_sequence):3d} VL {len(c.light_sequence):3d} aa"
                  + size)
        print()

    print("  Known targets, all routes")
    print("  " + "-" * 68)
    for gene in KNOWN_TARGETS:
        r = by_gene.get(gene)
        if r is None:
            continue
        print(f"    {gene:9s} {r.verdict:24s} entries {len(r.entries):3d}  "
              f"structure {len(r.structure):3d}  sequence {len(r.sequence):3d}")

    written = stage5.write_binders(records, manifest["stage4_hash"])
    print()
    print(f"  binders written to {written}")
    print(f"  configuration hash {stage5.configuration_hash(manifest['stage4_hash'], [r.gene for r in records])}")
    print("  affinity: NOT_CONNECTED for every candidate, and that is measured")

    if tripped:
        print()
        print(f"  STOPPING: {', '.join(tripped)} tripped.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
