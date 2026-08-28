"""One command for the whole platform: every stage, in order, from scratch.

Runs each stage's verifier as its own process and folds the results into a
single report. Separate processes on purpose — a stage that dies takes its own
interpreter with it and the run carries on to the next one, so a late failure
still reports every earlier stage rather than losing the lot.

``--fresh`` deletes the derived artifacts (the pairing decisions and the binder
cache) so they are rebuilt rather than read. The raw source caches under data/
are never touched: 11 GB of matrices and atlases are the input to this run, not
part of what it recomputes.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

#: Stage order is the pipeline order, and it is also a dependency order: the
#: pairing verifier reads what ranking wrote. Running these out of order would
#: read a stale artifact and call the agreement a result.
#:
#: Stage 5 in particular must precede 6 and 9. It is the only verifier that
#: exercises the binder retrieval route live; the two after it read the cache it
#: blessed, which is the path the API takes. Run standalone against a warm cache
#: those two cannot see a dead retrieval route — Stage 5 is what catches that,
#: and it is why it retrieves for real even though it costs five minutes.
STAGES = [
    ("1", "Design spec", "verify_schema.py"),
    ("2", "Surface proteome", "verify_surface.py"),
    ("3", "Target discovery", "verify_ranking.py"),
    ("4", "Target pairing", "verify_pairing.py"),
    ("4a", "Architecture routing", "verify_routing.py"),
    ("5", "Binder discovery", "verify_binders.py"),
    ("6", "Construct assembly", "verify_construct.py"),
    ("9", "Safety gate", "verify_safety.py"),
    ("10", "Developability", "verify_developability.py"),
    ("11", "Final ranking", "verify_ranking_final.py"),
    ("API", "HTTP surface", "verify_api.py"),
]

#: Derived artifacts. Everything else under data/ is a raw source cache.
DERIVED = ["stage4", "stage5"]

#: Criteria that are known to trip, each with the decision that accepted it.
#:
#: This exists so the exit code carries a signal. Five criteria trip on every
#: run, so a runner that simply failed on any trip would be red forever and a
#: genuine regression would land in a report nobody reads as changed.
#:
#: It is deliberately two-sided. A criterion tripping that is not listed here is
#: a regression. A criterion listed here that *stops* tripping is also
#: reported — a stale entry silently grants an exemption nothing needs, which
#: is how an accepted-limitations list turns into a place to hide things.
ACCEPTED = {
    ("3", "R13"): "Withdrawn. The two populations differ by construction "
                  "(breadth 51 vs 7) and the max-over-sources gate rewards "
                  "being unmeasured, so no combination rule reaches the 5x "
                  "limit. Replaced by R13-prime, which clears.",
    ("4", "P4"): "Coverage is span-confounded: f_AB tracks genomic span "
                 "(+0.68) more than expression (+0.20). Reported beside a "
                 "span-matched percentile and removed from partner selection.",
    ("4", "P7"): "One cleared pair contains HLA-A. Recorded rather than "
                 "filtered, because the pool is not curated by hand.",
    ("4", "P8"): "48.6% of cleared pairs stop clearing if an unmeasured "
                 "antigen saturates its organ. This is the cost of treating "
                 "missing as a third state instead of imputing it.",
    ("4", "P15"): "Partner choice is unstable under pool halving (71.4%). "
                  "The pairing stage is complete-with-limitations by decision.",
    ("4a", "A6"): "A positive pin written before the run. It expected MSLN to "
                  "route to an adaptor because it matches that row's condition "
                  "in words: serious normal-tissue expression. It does not, "
                  "because its measured risk 0.6366 is nearly twice the "
                  "declared terminable ceiling of 0.35. Admitting it needs a "
                  "ceiling near 0.65, which also admits about 120 others - a "
                  "clinical policy decision, not a code change. The ceiling "
                  "stays where the spec pinned it and A9 reports the whole "
                  "sweep so the trade is visible.",
}

# Both spacings are in use across the verifiers; matching only one would
# silently report zero criteria for half the stages.
_CRITERION = re.compile(r"^\s{2}(TRIPPED|clear)\s+([A-Za-z]+\d+[a-z']?):\s*(.*)$")
_SUMMARY = re.compile(r"^\s*(\d+)/(\d+) criteria clear")
_SCHEMA = re.compile(r"^checks passed: (\d+)/(\d+)")
# The spec verifier reports named checks rather than numbered criteria. Without
# this its 31 checks showed a count in the summary and nothing underneath it,
# which is the one shape a verification report must not have.
_CHECK = re.compile(r"^\s{2}(ok|FAIL)\s+(.+?): got (.+?)\s{2,}expected (.+)$")
# re.M matters: without it these anchors only match at the start of the whole
# captured output, so the surface stage parsed as zero criteria while exiting 0
# — a stage that passed and reported nothing.
_SETS = re.compile(r"^validation sets: (pass|FAIL)", re.M)
_DRIFT = re.compile(r"^filter decisions within .*: (yes|NO)", re.M)


class Stage:
    def __init__(self, number: str, name: str, script: str):
        self.number, self.name, self.script = number, name, script
        self.criteria: list[tuple[str, bool, str]] = []
        self.clear = 0
        self.total = 0
        self.code: int | None = None
        self.seconds = 0.0
        self.output = ""

    @property
    def ok(self) -> bool:
        return self.code == 0 and self.total > 0 and self.clear == self.total

    def parse(self) -> None:
        """Read the counts out of the verifier's own report.

        Deliberately not re-derived here. The verifier is the authority on its
        own criteria, and a second count living in this file would be free to
        disagree with it.
        """
        for line in self.output.splitlines():
            hit = _CRITERION.match(line)
            if hit:
                self.criteria.append(
                    (hit.group(2), hit.group(1) == "TRIPPED", hit.group(3).strip())
                )
                continue
            hit = _SUMMARY.match(line)
            if hit:
                self.clear, self.total = int(hit.group(1)), int(hit.group(2))
                continue
            hit = _CHECK.match(line)
            if hit:
                self.criteria.append((
                    f"check {len(self.criteria) + 1}",
                    hit.group(1) == "FAIL",
                    f"{hit.group(2)} — got {hit.group(3).strip()}, "
                    f"expected {hit.group(4).strip()}",
                ))
                continue
            hit = _SCHEMA.match(line)
            if hit:
                self.clear, self.total = int(hit.group(1)), int(hit.group(2))
        # The surface verifier reports two validation sets rather than numbered
        # criteria, so its pass/fail is folded into the same shape.
        sets = _SETS.search(self.output)
        drift = _DRIFT.search(self.output)
        if sets and drift:
            for label, value, good in (
                ("validation sets", sets.group(1), "pass"),
                ("count drift", drift.group(1), "yes"),
            ):
                self.criteria.append((label, value != good, value))
            self.total = 2
            self.clear = sum(1 for c in self.criteria if not c[1])

    def run(self, logs: Path, timeout: float) -> None:
        started = time.monotonic()
        try:
            proc = subprocess.run(
                [sys.executable, str(ROOT / self.script)],
                cwd=ROOT, capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=timeout,
            )
            self.code = proc.returncode
            self.output = (proc.stdout or "") + (proc.stderr or "")
        except subprocess.TimeoutExpired as exc:
            # A stalled socket in one verifier must not hang the whole run.
            # Recorded as a stage failure so the remaining stages still run and
            # still report, which is the reason these are separate processes.
            self.code = 124
            self.output = ((exc.stdout or "") if isinstance(exc.stdout, str)
                           else (exc.stdout or b"").decode("utf-8", "replace"))
            self.output += f"\n\nTIMEOUT after {timeout:.0f}s; stage abandoned.\n"
        self.seconds = time.monotonic() - started
        self.parse()
        # Written now, not at the end of the run. Buffering every transcript
        # until all ten stages finish means an interrupt during stage 5's five
        # minutes discards the nine that already succeeded.
        (logs / f"{self.script}.log").write_text(self.output, encoding="utf-8")


def api_result(output: str) -> list[str]:
    """The API's own closing block, quoted rather than recomputed."""
    lines = output.splitlines()
    for index, line in enumerate(lines):
        if "WHAT THE API RETURNS" in line:
            return [x.rstrip() for x in lines[index + 2:] if x.strip()]
    return []


def render(stages: list[Stage], elapsed: float, fresh: bool) -> str:
    total = sum(s.total for s in stages)
    clear = sum(s.clear for s in stages)
    provenance = ("Derived artifacts deleted and rebuilt" if fresh
                  else "Derived artifacts reused")
    lines = [
        "# Full run",
        "",
        f"{provenance}; raw sources read from `data/` unchanged.",
        "",
        f"**{clear}/{total} criteria clear across {len(stages)} stages**, "
        f"{elapsed / 60:.1f} minutes.",
        "",
        "| Stage | | Criteria | | Time |",
        "| --- | --- | --- | --- | --- |",
    ]
    for s in stages:
        mark = "clear" if s.ok else "**TRIPPED**"
        lines.append(
            f"| {s.number} | {s.name} | {s.clear}/{s.total} | {mark} "
            f"| {s.seconds:.0f}s |"
        )

    lines += ["", "## Every criterion", ""]
    for s in stages:
        lines += [f"### Stage {s.number} — {s.name} ({s.clear}/{s.total})", ""]
        if not s.criteria:
            lines += [f"No criteria parsed; the verifier exited {s.code}.", ""]
            continue
        for cid, tripped, detail in s.criteria:
            mark = "**TRIPPED**" if tripped else "clear"
            lines.append(f"- {mark} `{cid}` — {detail}")
        lines.append("")

    api = next((s for s in stages if s.number == "API"), None)
    block = api_result(api.output) if api else []
    if block:
        lines += ["## What the platform returns for this indication", "", "```"]
        lines += block
        lines += ["```", ""]

    tripped = [(s, c) for s in stages for c in s.criteria if c[1]]
    if tripped:
        lines += ["## Tripped", ""]
        for s, (cid, _t, detail) in tripped:
            note = ACCEPTED.get((s.number, cid))
            lines.append(f"- Stage {s.number} `{cid}` — {detail}")
            if note:
                lines.append(f"  - Accepted: {note}")
            else:
                lines.append("  - **Not on the accepted list. This is new.**")
    else:
        lines.append("Nothing tripped.")

    stale = stale_exemptions(stages)
    if stale:
        lines += ["", "## Stale exemptions", "",
                  "Listed as accepted, but did not trip on this run. Remove "
                  "them from the list or find out what changed.", ""]
        for number, cid in stale:
            lines.append(f"- Stage {number} `{cid}`")
    return "\n".join(lines) + "\n"


def unexpected(stages: list[Stage]) -> list[tuple[str, str, str]]:
    """Criteria that tripped without a decision on record."""
    return [(s.number, c[0], c[2]) for s in stages for c in s.criteria
            if c[1] and (s.number, c[0]) not in ACCEPTED]


def stale_exemptions(stages: list[Stage]) -> list[tuple[str, str]]:
    """Accepted entries that no longer trip, including ones that vanished.

    Reported for the same reason as a regression. An exemption that has stopped
    applying is an exemption nobody is checking, and it will still be sitting
    there covering something else the next time that identifier is used.

    The vanished case is the one that matters most here. R13 is *scheduled* to
    disappear — it was withdrawn and replaced. Checking only identifiers the run
    still emits would let its entry sit in this table forever, silently
    pre-approving whatever later reuses the name. A stage that reported nothing
    at all is excluded: it cannot distinguish "gone" from "never ran".
    """
    reported = {s.number for s in stages if s.total > 0}
    tripping = {(s.number, c[0]) for s in stages for c in s.criteria if c[1]}
    return sorted(key for key in ACCEPTED
                  if key[0] in reported and key not in tripping)


def project_interpreter() -> Path:
    """The venv interpreter, at whichever path this platform puts it."""
    for suffix in (("Scripts", "python.exe"), ("bin", "python")):
        candidate = ROOT.joinpath(".venv", *suffix)
        if candidate.exists():
            return candidate
    return ROOT / ".venv"


def preflight() -> str | None:
    """The interpreter this was launched with must be the one with the deps.

    Every stage runs under ``sys.executable``. Launched with the wrong Python
    all ten fail identically on the first import, and the report would say
    "no criteria parsed" ten times instead of naming the one cause.
    """
    # Read from requirements.txt rather than a list written here, which would
    # be free to drift from what the project actually declares. Only lines that
    # name a distribution: an option line is not a module.
    declared = []
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-")):
            continue
        name = re.split(r"[=<>!~\[;]", line)[0].strip()
        if name:
            declared.append(name)
    # find_spec rather than __import__: this process only shells out, and
    # importing numpy and h5py in full to ask whether they exist costs a second
    # and a hundred megabytes for an answer that needs neither.
    missing = [m for m in declared if importlib.util.find_spec(m) is None]
    if not missing:
        return None

    venv = project_interpreter()
    # Two different faults with two different remedies. Telling someone already
    # running the venv interpreter to run the venv interpreter is not advice.
    if Path(sys.executable).resolve() == venv.resolve():
        return (f"{sys.executable}\n  is the project interpreter but is "
                f"missing: {', '.join(missing)}\n"
                f"  Install them:\n    {venv} -m pip install -r requirements.txt")
    return (f"{sys.executable}\n  is missing: {', '.join(missing)}\n"
            f"  Run this with the project interpreter instead:\n"
            f"    {venv} run_all.py --fresh")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run every stage end to end.")
    parser.add_argument("--fresh", action="store_true",
                        help="delete derived artifacts so they are rebuilt")
    parser.add_argument("--report", default="reports/full-run.md")
    # Generous against the slowest stage (binder retrieval, ~320 s) but finite,
    # so a socket that stalls without closing fails one stage instead of the run.
    parser.add_argument("--timeout", type=float, default=1800.0,
                        help="seconds allowed per stage before it is abandoned")
    args = parser.parse_args()

    print("=" * 78)
    print("CAR-T DESIGN PLATFORM - FULL RUN")
    print("=" * 78)

    wrong = preflight()
    if wrong:
        print(f"\n  STOPPING before any stage runs. Wrong interpreter:\n  {wrong}")
        return 3
    print(f"  interpreter {sys.executable}")

    if args.fresh:
        for name in DERIVED:
            path = ROOT / "data" / name
            if path.exists():
                shutil.rmtree(path)
                print(f"  cleared derived artifact   data/{name}")
            else:
                print(f"  derived artifact absent    data/{name}")
    else:
        print("  reusing derived artifacts (--fresh rebuilds them)")
    raw = sorted(p.name for p in (ROOT / "data").iterdir()
                 if p.is_dir() and p.name not in DERIVED)
    print(f"  raw caches read, not rebuilt: {', '.join(raw)}")
    print()

    logs = ROOT / "reports" / "run-logs"
    logs.mkdir(parents=True, exist_ok=True)

    stages = [Stage(*s) for s in STAGES]
    started = time.monotonic()
    for stage in stages:
        label = "API" if stage.number == "API" else f"Stage {stage.number}"
        print(f"  {label:9s} {stage.name:20s} ", end="", flush=True)
        stage.run(logs, args.timeout)
        if stage.ok:
            verdict = "clear"
        elif stage.total:
            verdict = f"{stage.total - stage.clear} TRIPPED"
        else:
            verdict = f"no criteria, exit {stage.code}"
        print(f"{stage.clear:3d}/{stage.total:<3d} {verdict:16s} "
              f"{stage.seconds:6.1f}s")
    elapsed = time.monotonic() - started

    report = render(stages, elapsed, fresh=args.fresh)
    print()
    print(report)

    clear = sum(s.clear for s in stages)
    total = sum(s.total for s in stages)
    out = ROOT / args.report
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")

    # The structured run, so the page generator reads data rather than
    # re-parsing prose. A regex pinned to render()'s f-strings would break on a
    # formatting tweak and, worse, could half-match and render empty sections.
    (out.parent / (out.stem + ".json")).write_text(json.dumps({
        "clear": clear, "total": total, "minutes": round(elapsed / 60, 1),
        "fresh": args.fresh, "report": args.report,
        "accepted": {f"{n}/{c}": why for (n, c), why in ACCEPTED.items()},
        "unexpected": [{"stage": n, "id": c, "detail": d}
                       for n, c, d in unexpected(stages)],
        "stale": [{"stage": n, "id": c} for n, c in stale_exemptions(stages)],
        "api": api_result(next((s.output for s in stages
                                if s.number == "API"), "")),
        "stages": [{
            "number": s.number, "name": s.name, "clear": s.clear,
            "total": s.total, "ok": s.ok, "seconds": round(s.seconds, 1),
            "exit": s.code,
            "criteria": [{"id": i, "tripped": t, "detail": d}
                         for i, t, d in s.criteria],
        } for s in stages],
    }, indent=2), encoding="utf-8")

    print(f"written to {args.report} and {out.stem}.json")
    print("full transcripts in reports/run-logs/")

    # Two distinct failures, kept apart because conflating them produces a
    # message that is simply untrue. A verifier exits non-zero *because* a
    # criterion tripped, so a non-zero exit alongside tripped criteria is the
    # expected shape, not a fault.
    silent = [s for s in stages if s.total == 0]
    crashed = [s for s in stages
               if s.code != 0 and s.total > 0
               and not any(c[1] for c in s.criteria)]
    new = unexpected(stages)
    stale = stale_exemptions(stages)

    print()
    for s in silent:
        print(f"  NO CRITERIA  Stage {s.number} reported none at all "
              f"(exit {s.code}); an empty verification is not a passing one")
    for s in crashed:
        print(f"  UNEXPLAINED  Stage {s.number} exited {s.code} with "
              f"{s.clear}/{s.total} clear and nothing tripped")
    for number, cid, detail in new:
        print(f"  REGRESSION  Stage {number} {cid}: {detail}")
    for number, cid in stale:
        print(f"  STALE       Stage {number} {cid} is on the accepted list "
              "but did not trip")
    bad = silent or crashed or new or stale
    if not bad:
        accepted = sum(1 for s in stages for c in s.criteria if c[1])
        print(f"  {sum(s.clear for s in stages)}/{sum(s.total for s in stages)} "
              f"clear; {accepted} tripped, all {accepted} accepted on record.")
    return 0 if not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
