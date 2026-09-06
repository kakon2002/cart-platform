"""Tests the /v1/cart surface against its criteria.

D1 and D2 are the two that matter: they are what makes "the adapter translates
and never decides" a property rather than an intention. Both are blinded here,
because a guard that has never failed is not known to work.
"""

from __future__ import annotations

import ast
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

from car_pipeline.api import adapter
from car_pipeline.api.server import Handler
from car_pipeline.configs.registry import registered, resolve
from car_pipeline.stages import stage11, validation

HOST, PORT = "127.0.0.1", 8171
BASE = f"http://{HOST}:{PORT}"
ADAPTER_SOURCE = Path(adapter.__file__)

DEFAULT_INDICATION = "Pancreatic Ductal Adenocarcinoma"

# Modules that decide. The adapter may reach none of them directly.
FORBIDDEN_MODULES = {
    "stage1", "stage3", "stage4", "stage5", "stage6", "stage9", "stage10",
    "stage11", "stage12", "scoring", "routing", "validation", "construct_safety",
}


def call(method, path, body=None):
    """One HTTP call, returning status and decoded body."""
    data = json.dumps(body).encode() if body is not None else None
    request = urllib.request.Request(
        BASE + path, data=data, method=method,
        headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=1800) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read() or b"{}")


def direct_imports(path: Path) -> set[str]:
    """Every module this file imports by name, without importing it."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            names.add(module)
            names |= {f"{module}.{a.name}" for a in node.names}
    return names


def forbidden_in(names: set[str]) -> list[str]:
    """Which forbidden modules appear among these import names."""
    return sorted({m for name in names for m in FORBIDDEN_MODULES
                   if name.split(".")[-1] == m})


def vocabulary() -> set[str]:
    """Every literal the adapter is forbidden to spell."""
    return (set(stage11.DECISIONS) | set(stage11.GATE_STATUS.values())
            | {stage11.PASSED_ALL_GATES, validation.CONSERVATIVE,
               validation.INNOVATIVE})


def literals_in(source: str) -> set[str]:
    """Which forbidden values appear as string literals in this source."""
    tree = ast.parse(source)
    found = set()
    banned = vocabulary()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if node.value in banned:
                found.add(node.value)
    return found


def main() -> int:
    """Run the adapter criteria against a live server."""
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

    # ---------------- D1: the import graph ----------------
    names = direct_imports(ADAPTER_SOURCE)
    bad = forbidden_in(names)
    criterion("D1", bool(bad),
              f"the adapter imports {sorted(names)} directly and reaches no "
              f"deciding module"
              if not bad else f"the adapter imports {bad} directly")

    blinded = forbidden_in(names | {"car_pipeline.stages.stage11"})
    criterion("D1b", not blinded,
              "blinded: adding an import of stage11 trips D1"
              if blinded else "D1 does not trip when a stage import is added")

    # ---------------- D2: the vocabulary ----------------
    source = ADAPTER_SOURCE.read_text(encoding="utf-8")
    spelt = literals_in(source)
    criterion("D2", bool(spelt),
              f"none of the {len(vocabulary())} decision, gate or design-class "
              f"values appears as a literal in adapter source"
              if not spelt else f"the adapter spells {sorted(spelt)}")

    probe = source + '\nSMUGGLED = "ADVANCE"\n'
    criterion("D2b", not literals_in(probe),
              "blinded: adding the literal ADVANCE trips D2"
              if literals_in(probe) else "D2 does not trip on a decision literal")

    # ---------------- the live surface ----------------
    server = ThreadingHTTPServer((HOST, PORT), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()

    indication, _ = resolve(DEFAULT_INDICATION)
    cancer_type = indication.cancer_type

    # D4/D5: the input contract, through the adapter
    d4 = []
    for field, payload in (
            ("objective", {"indication": cancer_type, "objective": "rank"}),
            ("delivery_mode", {"indication": cancer_type,
                               "delivery_mode": "AUTO"})):
        st, b = call("POST", "/v1/cart/projects", payload)
        if st != 400 or b.get("field") != field:
            d4.append(f"{field} -> {st} naming {b.get('field')!r}")
    criterion("D4", bool(d4),
              "objective and delivery_mode are still refused by name through "
              "the adapter" if not d4 else "; ".join(d4))

    d5 = []
    st, b = call("POST", "/v1/cart/projects",
                 {"indication": cancer_type, "binder_mode": "RETRIEVAL_FIRST"})
    if st != 201:
        d5.append(f"RETRIEVAL_FIRST refused with {st}: {b.get('error','')[:60]}")
    for value in ("DE_NOVO", "OPTIMIZE", "ANYTHING"):
        st, b = call("POST", "/v1/cart/projects",
                     {"indication": cancer_type, "binder_mode": value})
        if st != 400 or b.get("field") != "binder_mode":
            d5.append(f"{value} -> {st} naming {b.get('field')!r}")
    criterion("D5", bool(d5),
              "binder_mode RETRIEVAL_FIRST is honoured and three other values "
              "are refused by name" if not d5 else "; ".join(d5))

    # D12: consulting the adapter must not consume the request body. The
    # handler reads the stream once and passes it on; reading it twice returns
    # empty the second time, which would hand the native endpoint an empty
    # request and refuse it for a field the caller did supply.
    st, native = call("POST", "/projects", {"cancer_type": cancer_type})
    criterion("D12", st != 201 or native.get("cancer_type") != cancer_type,
              "a native POST still receives its body after the adapter is "
              f"consulted ({st}, cancer_type {native.get('cancer_type')!r})"
              if st == 201 else
              f"the native path lost its body: {st} {native.get('error','')[:60]}")

    # D7: the two alias refusals
    d7 = []
    st, b = call("POST", "/v1/cart/projects",
                 {"indication": cancer_type, "project_id": "abcdef123456"})
    if st != 400 or b.get("field") != "project_id":
        d7.append(f"a canonical-shaped reference was accepted ({st})")
    st, first = call("POST", "/v1/cart/projects",
                     {"indication": cancer_type, "project_id": "WM-DUP-001"})
    st2, b2 = call("POST", "/v1/cart/projects",
                   {"indication": cancer_type, "project_id": "WM-DUP-001"})
    if st2 != 400 or not b2.get("existing_project_id"):
        d7.append(f"a duplicate reference was accepted ({st2})")
    criterion("D7", bool(d7),
              "a reference shaped like a canonical id and a duplicate "
              "reference are both refused by name" if not d7 else "; ".join(d7))

    # The run everything else reads
    st, project = call("POST", "/v1/cart/projects",
                       {"project_id": "WM-CART-PDAC-001",
                        "indication": cancer_type,
                        "target_mode": "DISCOVER",
                        "architecture_mode": "AUTO",
                        "binder_mode": "RETRIEVAL_FIRST",
                        "max_final_candidates": 5})
    reference = project["project_id"]
    canonical = project["canonical_project_id"]
    print(f"\n  project {reference} = {canonical}")

    st, run = call("POST", "/v1/cart/runs", {"project_id": reference})
    run_id = run.get("run_id")
    print(f"  run {run_id} {run.get('status')}", flush=True)
    while True:
        st, state = call("GET", f"/v1/cart/runs/{run_id}")
        if state.get("status") in ("COMPLETE", "FAILED"):
            break
        print(f"    {state.get('status'):<9} {state.get('stage')}  "
              f"{state.get('stages_complete')}/{state.get('stages_total')}",
              flush=True)
        time.sleep(20)
    print(f"  finished {state.get('status')}")

    # D3: progress is null with a reason, not a fraction
    criterion("D3", state.get("progress") is not None or not state.get("reasons"),
              f"progress is null with a reason; stage {state.get('stage')} at "
              f"{state.get('stages_complete')} of {state.get('stages_total')}"
              if state.get("progress") is None
              else f"progress reported as {state.get('progress')!r}")

    # D6: both names return identical bodies
    d6 = []
    for path in ("/candidates/rank", "/targets/discover"):
        _s1, a = call("POST", f"/v1/cart{path}", {"project_id": reference})
        _s2, b = call("POST", f"/v1/cart{path}", {"project_id": canonical})
        if json.dumps(a, sort_keys=True) != json.dumps(b, sort_keys=True):
            d6.append(path)
    criterion("D6", bool(d6),
              "an alias and its canonical id return byte-identical bodies from "
              "every view tried" if not d6 else f"bodies differ on {d6}")

    # D9: every documented path answers
    missing = []
    for path in adapter.PATHS:
        if "{run_id}" in path:
            st, _ = call("GET", f"/v1/cart/runs/{run_id}")
        elif "{binder_id}" in path:
            _s, blist = call("POST", "/v1/cart/binders/retrieve",
                             {"project_id": reference})
            bid = (blist.get("binders") or [{}])[0].get("binder_id")
            st, _ = call("GET", f"/v1/cart/binders/{bid}?project_id={reference}")
        else:
            st, _ = call("POST", path, {"project_id": reference,
                                        "target_ids": ["MSLN"],
                                        "candidate_ids": []})
        if st == 404:
            missing.append(path)
    criterion("D9", bool(missing),
              f"all {len(adapter.PATHS)} documented paths answer"
              if not missing else f"{len(missing)} return 404: {missing[:3]}")

    # D8: the three constants, through the adapter
    d8 = []
    for path, field, expected in (
            ("/v1/cart/structure/evaluate", "structure_status", "UNKNOWN"),
            ("/v1/cart/function/predict", "status",
             "INSUFFICIENT_VALIDATED_DATA"),
            ("/v1/cart/learning/calibrate", "status",
             "BLOCKED_INSUFFICIENT_DATA")):
        st, b = call("POST", path, {"construct_id": "C001"})
        if st != 200 or b.get(field) != expected:
            d8.append(f"{path} -> {st} {field}={b.get(field)!r}")
    criterion("D8", bool(d8),
              "the three constants return their declared status through the "
              "adapter" if not d8 else "; ".join(d8))

    # logic OR refused by name
    st, b = call("POST", "/v1/cart/targets/pair",
                 {"project_id": reference, "logic_modes": ["AND", "OR"]})
    criterion("D10", st != 400 or b.get("field") != "logic_modes",
              f"logic mode OR is refused by name: {b.get('error','')[:70]}"
              if st == 400 else f"OR accepted with {st}")

    # failed_gates carries at most one element
    st, gates = call("POST", "/v1/cart/gates/evaluate", {"project_id": reference})
    over = [c for c in gates.get("candidates", []) if len(c["failed_gates"]) > 1]
    criterion("D11", bool(over),
              f"{len(gates.get('candidates', []))} candidate(s), every "
              f"failed_gates carries at most one element"
              if not over else f"{len(over)} carry more than one")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        server.shutdown()
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHAT THE DASHBOARD RECEIVES")
    print("=" * 72)
    st, ranked = call("POST", "/v1/cart/candidates/rank",
                      {"project_id": reference})
    print(f"    run_status {ranked['run_status']}, "
          f"{ranked['eligible_candidate_count']} eligible, "
          f"{ranked['returned_candidate_count']} returned")
    for c in ranked["candidates"]:
        print(f"      {c['candidate_id']}  {c['targets']}  {c['architecture']}"
              f"  {c['decision']}  overall {c['scores'].get('overall')}")
    st, binders = call("POST", "/v1/cart/binders/retrieve",
                       {"project_id": reference, "target_id": "MSLN"})
    print(f"    MSLN binders retrieved: {binders['total_binders_found']}")
    for b in binders["binders"][:3]:
        print(f"      {b['binder_id']}  {b['route']:<9} {b['identifier']}"
              f"  match {b['target_match']}  affinity {b['affinity']}")

    st, made = call("POST", "/v1/cart/constructs/generate",
                    {"project_id": reference, "target_id": "FER1L6"})
    print(f"    constructs/generate for FER1L6 -> {st}, "
          f"{len(made.get('constructs') or [])} construct(s)")
    for c in (made.get("constructs") or [])[:2]:
        print(f"      {c.get('gene')}  {c.get('architecture')}  "
              f"{c.get('total_bp')} bp  {c.get('verdict')}")

    server.shutdown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
