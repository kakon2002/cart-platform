"""Tests the MSLN binder benchmark against its criteria.

Runs no comparison against the literature panel. That happens only after the
blind output is frozen and committed, and this file must never import or read
the answers.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

from benchmark import blind
from car_pipeline.api import pipeline, server
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import binder_check

ROOT = Path(__file__).resolve().parent
SABDAB = ROOT / "data/antibodies/sabdab_summary_all.csv"

BENCHMARK_TARGET = "MSLN"


def sabdab_rows() -> dict[str, list[dict]]:
    """Every antibody instance, grouped by entry."""
    grouped: dict[str, list[dict]] = {}
    with SABDAB.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            grouped.setdefault((row["PDB"] or "").lower(), []).append(row)
    return grouped


def main() -> int:
    """Run the benchmark criteria that do not need the answers."""
    print("loading sources", flush=True)
    surface, _ = load_surface()
    by_gene = {r.gene: r for r in surface if r.gene}
    entries = sabdab_rows()

    def entry(code: str, piped: bool = False) -> dict:
        """One instance of an entry. An entry has several and they differ."""
        rows = entries[f"pdb_0000{code.lower()}"]
        if piped:
            return next(r for r in rows if "|" in (r["antigen_name"] or ""))
        return rows[0]

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []
    checked: list[str] = []

    def criterion(cid, is_tripped, detail):
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    # ---------------- the blind rule, before anything else ----------------
    lit = blind.literal_violations()
    criterion("BM1", bool(lit),
              f"no held-out value appears in any of the "
              f"{len(set().union(*(blind.reachable_modules(e) for e in blind.ENTRY_POINTS)))} "
              f"modules reachable from the algorithm"
              if not lit else f"{lit}")

    caught = [held for held in ("SS1", "M11", "26.9")
              if blind.literal_violations(
                  extra_source={"car_pipeline.stages.binder_check":
                                f'X = "{held}"\n'})]
    ordinary = blind.literal_violations(
        extra_source={"car_pipeline.stages.binder_check": 'X = "Mesothelin"\n'})
    criterion("BM1b", len(caught) != 3 or bool(ordinary),
              "blinded: three held-out literals injected into a reachable "
              "module are each caught, and an ordinary literal is not"
              if len(caught) == 3 and not ordinary else
              f"caught {caught}, ordinary {ordinary}")

    reach = blind.answers_reachable()
    criterion("BM2", bool(reach),
              "the answers are unreachable from the algorithm's import graph"
              if not reach else f"reachable from {reach}")

    probes = [blind.answers_reachable(
                  extra_source={"car_pipeline.stages.binder_check": snippet})
              for snippet in ('P = "benchmark/answers/msln-literature.json"\n',
                              "from benchmark.answers import msln\n")]
    criterion("BM2b", not all(probes),
              "blinded: a reachable module naming the answers path or package "
              "is caught either way"
              if all(probes) else f"probes {probes}")

    a = {"binders": [{"binder_id": "b_1", "target_match": "PASS"}]}
    b = {"binders": [{"binder_id": "b_1", "target_match": "FAIL"}]}
    reordered = {"binders": [{"target_match": "PASS", "binder_id": "b_1"}]}
    moves = blind.fingerprint(a) != blind.fingerprint(b)
    stable = blind.fingerprint(reordered) == blind.fingerprint(a)
    criterion("BM3", not (moves and stable),
              "blinded: the fingerprint moves when the output changes and "
              "holds when only key order does"
              if moves and stable else f"moves {moves}, stable {stable}")

    refused = False
    try:
        blind.ANSWERS.parent.mkdir(parents=True, exist_ok=True)
        blind.ANSWERS.write_text('{"probe": true}', encoding="utf-8")
        try:
            blind.freeze({"probe": True})
        except blind.BlindViolation:
            refused = True
    finally:
        if blind.ANSWERS.exists():
            blind.ANSWERS.unlink()
    criterion("BM4", not refused,
              "blinded: freezing while the answers file exists is refused, so "
              "an output cannot be frozen after the answers are on disk"
              if refused else "freezing succeeded with the answers present")

    # BM5 -- the commit ordering. This can only be violated once the answers
    # exist, so until then it is satisfied by their absence rather than by a
    # check that cannot fail. The first draft was written `criterion("BM5",
    # False, ...)`, which is the shape this repository keeps finding: a line
    # that reports rather than tests.
    order = blind.committed_before_answers()
    answers_present = blind.ANSWERS.exists() or bool(order.get("answers_commit"))
    criterion("BM5", answers_present and not order["ordered"],
              (f"the answers are present and {order['reason']}"
               if answers_present else
               "the answers are not on disk and have never been committed, so "
               "no output could have been produced after reading them")
              + f" [frozen output committed: "
                f"{bool(order.get('frozen_commit'))}]")

    msln_names = binder_check.name_set(by_gene[BENCHMARK_TARGET])

    # BM6 -- the target-match verdicts, pinned before the check was written.
    # The first four are the specification's; the last two are cases the
    # retrieval run surfaced and are pinned for the same reason.
    PINS = [
        ("4F3F", BENCHMARK_TARGET, binder_check.PASS, False),
        ("7U8C", BENCHMARK_TARGET, binder_check.PASS, False),
        ("8H8J", "GPR35", binder_check.FAIL, False),
        ("1P4B", BENCHMARK_TARGET, binder_check.FAIL, False),
        ("7UED", BENCHMARK_TARGET, binder_check.UNKNOWN, False),
        ("8CZ8", BENCHMARK_TARGET, binder_check.PASS, True),
    ]
    wrong = []
    for code, gene, expected, piped in PINS:
        row = entry(code, piped)
        names = binder_check.name_set(by_gene[gene]) if gene in by_gene else set()
        got = binder_check.evaluate(
            row["antigen_name"], names,
            row["antigen_type"], row["antigen_species"])["target_match"]
        if got != expected:
            wrong.append(f"{code} expected {expected}, got {got}")
    criterion("BM6", bool(wrong),
              f"all {len(PINS)} pinned verdicts hold: "
              + ", ".join(f"{c}={e}" for c, _g, e, _p in PINS)
              if not wrong else "; ".join(wrong))

    # BM6a -- absent annotation reads UNKNOWN, never PASS and never FAIL.
    absent = [p for p in ("NA", "", "none", "-", "n/a")
              if binder_check.evaluate(p, msln_names)["target_match"]
              != binder_check.UNKNOWN]
    criterion("BM6a", bool(absent),
              "absent annotation reads UNKNOWN for every marker the source "
              "uses" if not absent else f"{absent} did not read UNKNOWN")

    # BM6b -- matching is per element. Blinded against the whole-string test,
    # which is what a naive check would do and what this case must defeat.
    piped_row = entry("8CZ8", piped=True)
    recorded = piped_row["antigen_name"]
    whole = recorded.strip().lower() in {n.lower() for n in msln_names}
    per_element = binder_check.evaluate(
        recorded, msln_names)["target_match"] == binder_check.PASS
    criterion("BM6b", whole or not per_element,
              f"per-element matching is load-bearing: {recorded!r} fails a "
              f"whole-string test and passes per element"
              if per_element and not whole else
              f"whole-string {whole}, per-element {per_element}")

    # BM6c -- a wrong antigen raises the flag rather than passing quietly.
    v = binder_check.evaluate(
        "Guanine nucleotide-binding protein subunit alpha-13", msln_names)
    criterion("BM6c",
              v["target_match"] != binder_check.FAIL or not v["wrong_antigen_flag"],
              "an unrelated antigen against the target's name set reads FAIL "
              "and raises wrong_antigen_flag"
              if v["wrong_antigen_flag"] else
              f"read {v['target_match']}, flag {v['wrong_antigen_flag']}")

    # BM7 -- species is carried and never gates.
    na_entries = [c for c, _g, _e, _p in PINS
                  if entry(c)["antigen_species"].strip().upper() == "NA"]
    u8c = binder_check.evaluate(
        entry("7U8C")["antigen_name"], msln_names,
        entry("7U8C")["antigen_type"],
        entry("7U8C")["antigen_species"])["target_match"]
    criterion("BM7", u8c != binder_check.PASS,
              f"{len(na_entries)} pinned entries record species NA and 7U8C, a "
              f"genuine {BENCHMARK_TARGET} entry among them, still passes"
              if u8c == binder_check.PASS else
              f"7U8C with species NA read {u8c}, so species is gating")

    # ---------------- the live retrieval ----------------
    print("\n  running the pipeline", flush=True)
    run = pipeline.run("Pancreatic Ductal Adenocarcinoma",
                       progress=lambda *a: None)
    rows = server._checked_rows(run, server._binder_rows(run))
    struct = [b for b in rows if b["route"] == "structure"]
    seq = [b for b in rows if b["route"] == "sequence"]

    # BM9 -- structural evidence only where target match passed.
    leaked = [b for b in rows
              if b["target_match"] != binder_check.PASS
              and (b.get("method") or b.get("antigen_chain"))
              and b.get("structural_evidence")]
    criterion("BM9", bool(leaked),
              "no binder carries structural evidence without a passing "
              "target match; the field is gated rather than filtered later"
              if not leaked else f"{len(leaked)} carry it regardless")

    # BM12 -- affinity is UNKNOWN everywhere while no source is connected.
    values = {b["affinity"] for b in rows}
    criterion("BM12", values != {"NOT_CONNECTED"},
              f"affinity is NOT_CONNECTED on all {len(rows)} retrieved "
              f"candidates; no connected release carries the column"
              if values == {"NOT_CONNECTED"} else f"affinity values {values}")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    # ---------------- what the check found ----------------
    print()
    print("=" * 72)
    print("SOURCE COVERAGE, BEFORE ANY RESULT")
    print("=" * 72)
    print(f"    affinity values present: 0 of {len(rows)}")
    print("    No connected evidence release carries an affinity, KD or")
    print("    free-energy column. Fold error against a literature value and")
    print("    rank correlation against a literature ordering are therefore")
    print("    structurally uncomputable rather than computed and failed.")

    print()
    print("=" * 72)
    print("B3 ACROSS THE POOL")
    print("=" * 72)
    print(f"    {len(rows)} retrieved binder(s): {len(struct)} structure "
          f"route, {len(seq)} sequence route")
    print(f"    verdicts        {dict(Counter(b['target_match'] for b in rows))}")
    print(f"    structure route {dict(Counter(b['target_match'] for b in struct))}")
    flagged = [b for b in struct if b.get("wrong_antigen_flag")]
    print(f"    wrong_antigen_flag raised on {len(flagged)} of {len(struct)} "
          f"structure-route rows")

    print()
    print(f"  {BENCHMARK_TARGET}, the benchmark target")
    msln = [b for b in rows if b["target_id"] == BENCHMARK_TARGET]
    print(f"    {len(msln)} row(s) "
          f"{dict(Counter(b['target_match'] for b in msln))}")
    for b in msln:
        print(f"      {b['binder_id']}  {b['route']:<9} "
              f"{str(b['identifier'])[:14]:<14} {b['target_match']:<8} "
              f"{str(b['recorded_antigens'])[:44]}")

    print()
    print("=" * 72)
    print("SEQUENCE COVERAGE, WHICH BOUNDS B4 AND B5")
    print("=" * 72)
    withseq = [b for b in msln if b["sequence_available"]]
    print(f"    {BENCHMARK_TARGET}: {len(withseq)} of {len(msln)} rows carry a "
          f"binder sequence")
    print(f"      structure route "
          f"{sum(1 for b in msln if b['route'] == 'structure' and b['sequence_available'])}"
          f" of {sum(1 for b in msln if b['route'] == 'structure')}")
    print(f"      sequence route  "
          f"{sum(1 for b in msln if b['route'] == 'sequence' and b['sequence_available'])}"
          f" of {sum(1 for b in msln if b['route'] == 'sequence')}")
    print(f"    pool-wide: {sum(1 for b in rows if b['sequence_available'])} "
          f"of {len(rows)}")
    print()
    print("    Sanitation and antibody annotation therefore apply to the")
    print("    sequence-route rows only. The structure route carries an")
    print("    identifier and an annotation and no residues, which is what")
    print("    the target-match and structural-evidence checks read and is")
    print("    nothing for a sequence step to work on.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
