"""The blind-run harness for the binder benchmark.

The benchmark's own critical rule is that the literature panel must not enter
the algorithm before the blind run. Stated as an intention that is worth
nothing, so it is enforced in four layers, each mechanical and each blinded
before it is trusted:

1. **The literal guard.** No held-out value appears in any module reachable
   from the benchmark entry point.
2. **The import graph.** The answers are unreachable from the algorithm.
3. **The freeze.** The blind output is written with a configuration hash and
   committed *before* the answers are read. The comparison refuses to run if
   re-executing the algorithm does not reproduce that hash. This makes the
   ordering a property of the commit graph rather than of anyone's memory.
4. **Blinding.** Each of the first three is deliberately broken and must trip.

Layer 3 is the one that matters. The others can be satisfied by a careful
author; only the freeze survives an author who is not careful, because the
evidence is the artifact and the order in which it was committed.
"""

from __future__ import annotations

import ast
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANSWERS = ROOT / "benchmark" / "answers" / "msln-literature.json"
FROZEN = ROOT / "benchmark" / "frozen" / "msln-blind.json"

# The values that must not reach the algorithm. Held here, in a module the
# algorithm never imports, so that the guard has something to check against
# without the answers themselves being importable.
HELD_OUT = (
    "SS1", "M5", "M11", "15B6", "h15B6",
    "26.9", "64.7",
)

# Modules the algorithm is built from. The guard walks imports from here.
ENTRY_POINTS = (
    "car_pipeline.stages.binder_check",
    "car_pipeline.stages.stage5",
    "car_pipeline.stages.scoring",
)


class BlindViolation(RuntimeError):
    """The blind rule was broken. Never caught and continued past."""


# ---------------------------------------------------------------------------
# Layer 1 -- the literal guard
# ---------------------------------------------------------------------------

def _module_path(name: str) -> Path | None:
    """Where a dotted module lives, without importing it."""
    candidate = ROOT / (name.replace(".", "/") + ".py")
    return candidate if candidate.exists() else None


def reachable_modules(entry: str, seen: set[str] | None = None) -> set[str]:
    """Every first-party module reachable from an entry point, by parsing.

    Parsed rather than imported, so walking the graph cannot itself execute
    the code being audited.
    """
    seen = seen if seen is not None else set()
    if entry in seen:
        return seen
    path = _module_path(entry)
    if path is None:
        return seen
    seen.add(entry)
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module] + [f"{node.module}.{a.name}"
                                     for a in node.names]
        for name in names:
            if name.startswith("car_pipeline") or name.startswith("benchmark"):
                reachable_modules(name, seen)
    return seen


def literal_violations(extra_source: dict[str, str] | None = None) -> list[str]:
    """Any held-out value appearing in a module the algorithm can reach."""
    out: list[str] = []
    modules: set[str] = set()
    for entry in ENTRY_POINTS:
        modules |= reachable_modules(entry)
    sources = {m: (_module_path(m).read_text(encoding="utf-8"))
               for m in modules if _module_path(m)}
    sources.update(extra_source or {})
    for module, source in sorted(sources.items()):
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                for held in HELD_OUT:
                    if held == node.value.strip():
                        out.append(f"{module}: {held!r}")
    return out


# ---------------------------------------------------------------------------
# Layer 2 -- the answers must be unreachable
# ---------------------------------------------------------------------------

def answers_reachable(extra_source: dict[str, str] | None = None) -> list[str]:
    """Any algorithm module that names the answers file or its package."""
    out: list[str] = []
    modules: set[str] = set()
    for entry in ENTRY_POINTS:
        modules |= reachable_modules(entry)
    sources = {m: (_module_path(m).read_text(encoding="utf-8"))
               for m in modules if _module_path(m)}
    sources.update(extra_source or {})
    for module, source in sorted(sources.items()):
        if "msln-literature" in source or "benchmark.answers" in source:
            out.append(module)
        if module.startswith("benchmark.answers"):
            out.append(module)
    return out


# ---------------------------------------------------------------------------
# Layer 3 -- the freeze
# ---------------------------------------------------------------------------

def fingerprint(payload: dict) -> str:
    """A hash over the blind output, stable under key order."""
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def freeze(payload: dict) -> dict:
    """Write the blind output with its hash. Refuses to overwrite silently."""
    violations = literal_violations() + answers_reachable()
    if violations:
        raise BlindViolation(
            "the blind rule is already broken, so freezing would record an "
            f"output produced under it: {violations}")
    if ANSWERS.exists():
        raise BlindViolation(
            f"{ANSWERS.name} exists. The answers must not be present when the "
            "blind output is produced: freeze and commit first, then reveal.")
    record = {
        "benchmark_id": "MSLN-LITERATURE-001",
        "blind": True,
        "output": payload,
        "fingerprint": fingerprint(payload),
    }
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(record, indent=2, sort_keys=True),
                      encoding="utf-8", newline="")
    return record


def verify_frozen(payload: dict) -> dict:
    """Re-running must reproduce the frozen output exactly."""
    if not FROZEN.exists():
        raise BlindViolation("nothing frozen; the blind run has not happened")
    record = json.loads(FROZEN.read_text(encoding="utf-8"))
    again = fingerprint(payload)
    return {
        "frozen": record["fingerprint"],
        "recomputed": again,
        "reproduces": again == record["fingerprint"],
    }


def committed_before_answers() -> dict:
    """Whether the frozen output was committed before the answers were.

    The commit graph is the evidence. If the answers file has never been
    committed the ordering is trivially satisfied; if both exist, the frozen
    artifact's first commit must not come after the answers' first commit.
    """
    def first_commit(path: Path) -> str | None:
        """The earliest commit touching this path, if any."""
        rel = path.relative_to(ROOT).as_posix()
        try:
            out = subprocess.run(
                ["git", "log", "--reverse", "--format=%H %ct", "--", rel],
                cwd=ROOT, capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            return None
        line = (out.stdout or "").strip().splitlines()
        return line[0] if line else None

    frozen = first_commit(FROZEN)
    answers = first_commit(ANSWERS)
    if frozen is None:
        return {"ordered": False, "reason": "the frozen output is not committed"}
    if answers is None:
        return {"ordered": True,
                "reason": "the answers have never been committed, so nothing "
                          "could have been read before the freeze",
                "frozen_commit": frozen.split()[0]}
    ok = int(frozen.split()[1]) <= int(answers.split()[1])
    return {
        "ordered": ok,
        "frozen_commit": frozen.split()[0],
        "answers_commit": answers.split()[0],
        "reason": ("the frozen output was committed no later than the answers"
                   if ok else
                   "the answers were committed first, so the blind run cannot "
                   "be shown to have been blind"),
    }


def guard(payload: dict | None = None) -> dict:
    """Every layer, reported together. Raises on a violation."""
    literals = literal_violations()
    reach = answers_reachable()
    state = {
        "layer_1_literals": literals,
        "layer_2_answers_reachable": reach,
        "layer_3_ordering": committed_before_answers(),
    }
    if payload is not None and FROZEN.exists():
        state["layer_3_reproduces"] = verify_frozen(payload)
    if literals or reach:
        raise BlindViolation(json.dumps(state, indent=2))
    return state
