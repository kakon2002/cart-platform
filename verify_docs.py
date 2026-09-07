"""Checks the numbers written into the documents against the code.

A figure in prose is a claim with no criterion. Nothing re-evaluates a
sentence, so it accretes authority by remaining on the page while the code
underneath it moves. This project shipped three such figures in one week: a
criteria total that went stale when two idioms were converted, a cache count
that went stale when a source was tracked, and a download size that had been
wrong since it was typed.

These criteria give the prose a criterion. Every count asserted in the
documents is read out of the file and compared against the live value.
Timings are not checked -- they are dated measurements rather than claims, and
they drift with the machine.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
    "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12,
    "thirteen": 13, "fourteen": 14, "fifteen": 15, "sixteen": 16,
    "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
}


def spelled(text: str) -> int | None:
    """A number written as a word, if that is what this is."""
    return WORDS.get(text.strip().lower())


def main() -> int:
    """Run the documentation criteria."""
    runbook = (ROOT / "RUNBOOK.md").read_text(encoding="utf-8")
    handover = (ROOT / "HANDOVER.md").read_text(encoding="utf-8")
    orientation = (ROOT / "ORIENTATION.md").read_text(encoding="utf-8")

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

    import bootstrap
    from car_pipeline.api import pipeline
    from car_pipeline.configs.registry import INDICATIONS

    # --- V1: the count of shared sources, wherever it is asserted ----------
    shared = len(bootstrap.SHARED_SOURCES)
    asserted = set(re.findall(r"(\d+)/(\d+) shared sources usable", runbook))
    wrong = sorted({b for _a, b in asserted if int(b) != shared})
    criterion("V1", bool(wrong),
              f"every 'n/n shared sources usable' in the runbook uses "
              f"{shared}, matching bootstrap.SHARED_SOURCES"
              if not wrong else
              f"the runbook says out of {wrong} where the code has {shared}")

    # --- V2: the MISSING count a first bootstrap prints --------------------
    # Shared sources minus those deferred, plus one line per per-indication
    # cache per indication that declares an atlas.
    deferred = sum(1 for _n, _d, cost in bootstrap.SHARED_SOURCES
                   if "first run" in cost)
    with_atlas = sum(1 for i in INDICATIONS.values() if i.atlas is not None)
    expected_missing = (shared - deferred) + len(bootstrap.PER_INDICATION) * with_atlas
    m = re.search(r"Expected on a first run: (\w+) `MISSING` lines",
                  runbook) or re.search(
                      r"Expected on a first run: (\w+) `MISSING` lines", runbook)
    stated = spelled(m.group(1)) if m else None
    criterion("V2", stated != expected_missing,
              f"the runbook says {m.group(1) if m else '?'} MISSING lines and "
              f"the code yields {expected_missing} "
              f"({shared} shared - {deferred} deferred + "
              f"{len(bootstrap.PER_INDICATION)}x{with_atlas} per-indication)"
              if stated == expected_missing else
              f"the runbook says {stated}, the code yields {expected_missing}")

    # --- V3: the dependency count ------------------------------------------
    reqs = [l for l in (ROOT / "requirements.txt").read_text(
        encoding="utf-8").splitlines() if l.strip() and not l.startswith("#")]
    m = re.search(r"(\w+) direct dependencies", runbook)
    stated = spelled(m.group(1)) if m else None
    criterion("V3", stated != len(reqs),
              f"the runbook says {m.group(1)} direct dependencies and "
              f"requirements.txt pins {len(reqs)}"
              if stated == len(reqs) else
              f"the runbook says {stated}, requirements.txt pins {len(reqs)}")

    # --- V4: the suite total, against the last recorded run ----------------
    run = json.loads((ROOT / "reports/full-run.json").read_text(encoding="utf-8"))
    m = re.search(r"\*\*(\d+) of (\d+)\s*\n?criteria clear\*\*", runbook)
    if not m:
        m = re.search(r"\*\*(\d+) of (\d+)", runbook)
    stated = (int(m.group(1)), int(m.group(2))) if m else None
    actual = (run["clear"], run["total"])
    criterion("V4", stated != actual,
              f"the runbook's {stated[0]} of {stated[1]} matches the last "
              f"recorded run" if stated == actual else
              f"the runbook says {stated}, the last run was {actual}")

    # --- V5: the expected trips, by name -----------------------------------
    unexpected = {f"{u['stage']}/{u['id']}" for u in run["unexpected"]}
    named = set(re.findall(r"`(\d+[a-z]?/[A-Z]+\d+)`", runbook)) or \
        set(re.findall(r"`(\d+[a-z]?/[A-Z]+\d+)`", runbook))
    criterion("V5", not unexpected <= named,
              f"every open trip {sorted(unexpected)} is named in the runbook"
              if unexpected <= named else
              f"the run reports {sorted(unexpected - named)} which the runbook "
              f"does not name")

    # --- V6: the accepted exemptions, by name ------------------------------
    accepted = {c for _n, c in bootstrap and () } if False else set()
    import run_all
    accepted = {cid for _stage, cid in run_all.ACCEPTED}
    named_all = set(re.findall(r"`([A-Z]+\d+)`", runbook))
    missing = sorted(accepted - named_all)
    criterion("V6", bool(missing),
              f"all {len(accepted)} accepted exemptions are named in the runbook"
              if not missing else f"the runbook does not name {missing}")

    # --- V7: the stage count in the progress example -----------------------
    stages = len(pipeline.STAGES)
    m = re.search(r'"stages_complete": \d+,\s*\n\s*"stages_total": (\d+)', runbook)
    if not m:
        m = re.search(r"stages_total\D{0,4}(\d+)", runbook)
    stated = int(m.group(1)) if m else None
    criterion("V7", stated != stages,
              f"the runbook's stages_total {stated} matches pipeline.STAGES"
              if stated == stages else
              f"the runbook says {stated}, pipeline.STAGES has {stages}")

    # --- V8: no document states a bare suite total in prose ----------------
    # ORIENTATION deliberately points at the report rather than naming one.
    bare = re.findall(r"\*\*\d{3} criteria\*\*", orientation + handover)
    criterion("V8", bool(bare),
              "neither ORIENTATION nor HANDOVER writes a suite total into "
              "prose; both point at reports/full-run.md"
              if not bare else f"a bare total appears: {bare}")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
