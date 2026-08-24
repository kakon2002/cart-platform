"""Runs the developability stage and tests it against its criteria.

The pins are synthetic controls with hand-checkable answers, so the arithmetic is
verified independently of whatever the data happens to contain. A stage whose
only checks run over its own output can be uniformly wrong and still pass.
"""

from __future__ import annotations

import sys

from car_pipeline.stages import stage4, stage5, stage10


def main() -> int:
    print("loading upstream", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)
    # From the persisted binder artifact, not the network. Re-querying two
    # hundred accessions for every downstream stage is slow and loses the whole
    # run to one dropped connection.
    records, binder_manifest = stage5.read_binders()
    binders = {r.gene: r for r in records}
    rows, status = stage10.assess(binders)
    with_sequence = sum(1 for r in records
                        for c in r.sequence if c.heavy_sequence)
    print(f"  {len(rows)} binder sequences scored, status {status}")

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        if is_tripped:
            tripped.append(cid)

    basic = stage10.isoelectric_point("KKKKKKKKKK")
    acidic = stage10.isoelectric_point("EEEEEEEEEE")
    criterion("D1", not (basic > 9.0 and acidic < 5.0),
              f"poly-K pI {basic}, poly-E pI {acidic} — the charge model orders "
              "the two known answers correctly")

    yes = stage10.glycosylation_sequons("AAANSTAAA")
    no = stage10.glycosylation_sequons("AAANPTAAA")
    criterion("D2", not (len(yes) == 1 and len(no) == 0),
              f"NST yields {len(yes)} sequon, NPT yields {len(no)} — the proline "
              "exclusion holds")

    three = stage10.score("CACAC", "control", "control")
    criterion("D3", three.cysteine_parity != "odd",
              f"a 3-cysteine control reports parity {three.cysteine_parity!r}, "
              "and never 'unpaired: 0'")

    bad = [r.gene for r in rows
           if not (1.0 <= r.isoelectric_point <= 14.0) or not (0 <= r.flag_count <= 5)]
    criterion("D4", bool(bad),
              "every scored binder has a pI in 1..14 and 0..5 flags"
              if not bad else f"out of range: {bad[:5]}")

    criterion("D5", len(rows) != with_sequence,
              f"{len(rows)} rows against {with_sequence} binders carrying a sequence")

    summed = any(hasattr(r, "score") or hasattr(r, "developability_score")
                 for r in rows)
    criterion("D6", summed,
              "no liability is summed into a single score; flags are counted and "
              "listed")

    print("=" * 72)
    print(f"  {6 - len(tripped)}/6 criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHAT THE SEQUENCES CARRY")
    print("=" * 72)
    if status == stage10.NOTHING_TO_SCORE:
        print("    NOTHING_TO_SCORE — no binder carries a sequence.")
        print("    Reported as a status rather than as an empty table, which would")
        print("    read as 'nothing had liabilities'.")
        return 0

    counts: dict[str, int] = {}
    for r in rows:
        for kind in r.kinds:
            counts[kind] = counts.get(kind, 0) + 1
    print(f"    binder sequences scored          {len(rows)}")
    print(f"    carrying at least one flag       "
          f"{sum(1 for r in rows if r.flags)}")
    print(f"    clean on all five               "
          f"{sum(1 for r in rows if not r.flags)}")
    print()
    print("    liability by kind:")
    for kind, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"      {kind:44s} {n:4d}")

    universal = [k for k, n in counts.items() if n == len(rows)]
    if universal:
        print()
        print("    A FLAG THAT FIRES ON EVERY INPUT CARRIES NO INFORMATION.")
        for kind in universal:
            print(f"      {kind!r} fired on all {len(rows)} sequences.")
        print("      The threshold was fixed before the run and is NOT adjusted")
        print("      now that its output is visible — that is the tuning this")
        print("      project forbids. It is recorded as uninformative at its")
        print("      stated setting, which is a result about the threshold rather")
        print("      than about the sequences. Every antibody variable region has")
        print("      hydrophobic windows; discriminating between them needs a")
        print("      measure of surface exposure, which sequence alone cannot give.")

    print()
    print("  The binders for the two targets under discussion")
    print("  " + "-" * 68)
    for gene in ("MSLN", "CLDN18"):
        for r in rows:
            if r.gene != gene:
                continue
            print(f"    {gene:8s} {r.binder:16s} {r.residues:4d} aa  "
                  f"pI {r.isoelectric_point:6.3f}  charge {r.net_charge:+7.3f}  "
                  f"Cys {r.cysteines} ({r.cysteine_parity})  "
                  f"sequons {len(r.sequons)}  APR {len(r.apr_starts)}")
            for kind, detail in r.flags:
                print(f"        - {kind}: {detail}")

    print()
    print(f"  configuration hash "
          f"{stage10.configuration_hash(manifest['stage4_hash'], [r.gene for r in rows])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
