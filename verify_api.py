"""Exercises the API end to end against a live server."""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

from car_pipeline.api.server import Handler
from car_pipeline.configs.registry import registered, resolve
from car_pipeline.stages import stage3

HOST, PORT = "127.0.0.1", 8137
BASE = f"http://{HOST}:{PORT}"


DEFAULT_INDICATION = "Pancreatic Ductal Adenocarcinoma"


def call(method: str, path: str, body: dict | None = None):
    """Make one HTTP call and return the status and decoded body."""
    data = json.dumps(body).encode() if body is not None else None
    request = urllib.request.Request(
        BASE + path, data=data, method=method,
        headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read() or b"{}")


def main() -> int:
    """Run the API criteria against a live server."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--indication", default=DEFAULT_INDICATION,
        help="cancer type to exercise; resolved through the registry, so an "
             "alias works. Configured: "
             + "; ".join(row["cancer_type"] for row in registered()))
    args = parser.parse_args()

    indication, _project = resolve(args.indication)
    cancer_type = indication.cancer_type
    print(f"  indication {cancer_type}")

    server = ThreadingHTTPServer((HOST, PORT), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    print(f"  server up on {BASE}")

    tripped: list[str] = []

    checked: list[str] = []
    def criterion(cid, is_tripped, detail):
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)

    status, project = call("POST", "/projects",
                           {"cancer_type": cancer_type})
    pid = project.get("project_id", "")
    criterion("A1", status != 201 or project.get("target_antigen") is not None,
              f"project created ({status}), target_antigen "
              f"{project.get('target_antigen')} and discovery mode "
              f"{project.get('discovery_mode')}")

    status, early = call("GET", f"/projects/{pid}/constructs")
    criterion("A2", status != 409 or early.get("status") != "RUN_NOT_COMPLETE",
              f"a view before any run answers {status} {early.get('status')} "
              "with instructions, not an empty list")

    status, job = call("POST", f"/projects/{pid}/runs")
    jid = job.get("job_id", "")
    criterion("A3", status != 202 or not jid,
              f"a run returns {status} with job {jid} rather than blocking")

    print()
    print("  polling", flush=True)
    seen, last = [], None
    deadline = time.time() + 3600
    while time.time() < deadline:
        _s, state = call("GET", f"/jobs/{jid}")
        if state.get("stage") and state["stage"] != last:
            last = state["stage"]
            seen.append(last)
            print(f"    {state['status']:9s} {last}  {state.get('note','')}")
        if state.get("status") in ("complete", "failed"):
            break
        time.sleep(5)

    criterion("A4", state.get("status") != "complete",
              f"job finished {state.get('status')} after stages {seen}"
              + (f" — {state.get('error')}" if state.get("error") else ""))
    if state.get("status") != "complete":
        print(state.get("trace", ""))
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    status, constructs = call("GET", f"/projects/{pid}/constructs")

    states = {row.get("state") for row in constructs.get("constructs", [])}
    counted = constructs.get("complete", 0) + constructs.get("awaiting_binder", 0)
    criterion("A5",
              status != 200
              or not constructs.get("reasons")
              or counted != constructs.get("buildable")
              or (constructs.get("awaiting_binder")
                  and "AWAITING_BINDER" not in states),
              f"{status} {constructs.get('status')}: "
              f"{constructs.get('buildable')} buildable = "
              f"{constructs.get('complete')} complete + "
              f"{constructs.get('awaiting_binder')} awaiting a binder; "
              f"{constructs.get('over_budget')} over budget, "
              f"{len(constructs.get('reasons', []))} reasons")

    status, result = call("GET", f"/projects/{pid}/result")
    total = sum(step["dropped"] for step in result.get("attrition", []))

    reached = result.get("reached_the_end", 0)
    criterion("A6",
              status != 200
              or total + reached != result.get("pool_size")
              or result.get("complete", 0) + result.get("awaiting_binder", 0)
                 != reached,
              f"end state {result.get('status')}, attrition accounts for "
              f"{total} + {reached} of {result.get('pool_size')}; {reached} "
              f"reached = {result.get('complete')} complete + "
              f"{result.get('awaiting_binder')} awaiting")

    status, targets = call("GET", f"/projects/{pid}/targets")

    rows = targets.get("targets") or []
    top = rows[0] if rows else {}
    breakdown = top.get("breakdown") or {}
    components = set(stage3.WEIGHTS)
    criterion("A7",
              status != 200 or not rows or not top.get("gene")
              or set(breakdown) != components,
              f"top target {top.get('gene')} ranked {top.get('rank')} carries "
              f"all {len(components)} scoring components"
              if rows and set(breakdown) == components else
              "no ranked target returned" if not rows else
              f"breakdown is {sorted(breakdown)}, expected {sorted(components)}")

    status, pairs = call("GET", f"/projects/{pid}/pairs")
    prows = pairs.get("pairs") or [{}]
    first = prows[0]
    criterion("A8",
              status != 200 or first.get("coverage_span_percentile") is None,
              "pairs carry the span percentile beside the raw fraction "
              f"({first.get('coverage_f_ab')} at percentile "
              f"{first.get('coverage_span_percentile')})")

    s_unknown, unknown = call("GET", "/projects/ffffffffffff/result")
    _s, fresh = call("POST", "/projects",
                     {"cancer_type": cancer_type})
    s_unrun, unrun = call("GET", f"/projects/{fresh.get('project_id')}/result")
    criterion("A10",
              s_unknown != 404
              or unknown.get("status") != "NOT_FOUND"
              or not unknown.get("reasons")
              or s_unrun != 409
              or unrun.get("status") != "RUN_NOT_COMPLETE",
              f"an unknown project answers {s_unknown} "
              f"{unknown.get('status')} and one that exists without a finished "
              f"run answers {s_unrun} {unrun.get('status')}: a client can tell "
              "a bad id from a run in progress")

    status, evidence = call("GET", f"/projects/{pid}/evidence/MSLN")
    stages = [k for k in evidence if k.startswith("stage")]
    criterion("A9",
              status != 200 or len(stages) < 6,
              f"evidence trail for MSLN spans {len(stages)} stages: "
              f"{', '.join(s.split('_')[0] for s in stages)}")

    # W10: every field the document names is either honoured or refused by
    # name. Both directions are exercised -- a field wrongly refused fails this
    # as surely as one wrongly accepted.
    must_refuse = [
        ("objective", {"cancer_type": cancer_type, "objective": "rank safely"}),
        ("delivery_mode", {"cancer_type": cancer_type, "delivery_mode": "AUTO"}),
        ("architecture_mode", {"cancer_type": cancer_type,
                               "architecture_mode": "OR"}),
        ("architecture_mode", {"cancer_type": cancer_type,
                               "architecture_mode": "AND-NOT"}),
        ("max_final_candidates", {"cancer_type": cancer_type,
                                  "max_final_candidates": 0}),
    ]
    must_accept = [
        {"cancer_type": cancer_type, "max_final_candidates": 3},
        {"cancer_type": cancer_type, "architecture_mode": "AUTO"},
        {"cancer_type": cancer_type, "architecture_mode": "ADAPTOR"},
        {"cancer_type": cancer_type, "target_mode": "DISCOVER"},
        {"cancer_type": cancer_type, "project_id": "WM-CART-001"},
    ]
    w10_bad = []
    for field, payload in must_refuse:
        st, b = call("POST", "/projects", payload)
        if st != 400:
            w10_bad.append(f"{field} accepted with {st}, not refused")
        elif b.get("field") != field:
            w10_bad.append(f"{field} refused but the body names "
                           f"{b.get('field')!r}")
        elif not b.get("error"):
            w10_bad.append(f"{field} refused with no message")
    for payload in must_accept:
        st, b = call("POST", "/projects", payload)
        if st != 201:
            sent = [k for k in payload if k != "cancer_type"]
            w10_bad.append(f"{sent} refused with {st}: {b.get('error', '')[:60]}")
    # And the document's own example, which must refuse on a named field
    # rather than on a bare 400.
    st, b = call("POST", "/projects", {
        "project_id": "WM-CART-PDAC-001", "indication": cancer_type,
        "objective": "discover and rank safe CAR-T candidates",
        "target_mode": "DISCOVER", "architecture_mode": "AUTO",
        "delivery_mode": "AUTO", "max_final_candidates": 5})
    if st != 400 or not b.get("field") or not b.get("remove"):
        w10_bad.append("the document's example JSON is not refused by name")
    criterion("W10", bool(w10_bad),
              f"{len(must_refuse)} unhonoured field(s) refused by name and "
              f"{len(must_accept)} honoured field(s) accepted; the document's "
              f"own example is refused naming {b.get('field')!r} and saying "
              f"to remove {b.get('remove')}"
              if not w10_bad else "; ".join(w10_bad[:3]))

    # W12: the output contract carries every field section 16 names, and the
    # ones the platform cannot fill are named absent rather than emitted empty.
    st, contract = call("GET", f"/projects/{pid}/contract")
    required = ("run_status", "eligible_candidate_count", "candidates",
                "excluded_candidates", "next_best_experiments", "audit_id")
    w12_bad = [k for k in required if k not in contract] if st == 200 else \
        [f"GET /contract returned {st}"]
    if st == 200:
        if contract.get("next_best_experiments") == []:
            w12_bad.append("next_best_experiments is an empty list, which says "
                           "no experiment is recommended rather than that none "
                           "was computed")
        if not contract.get("next_best_experiments_state"):
            w12_bad.append("next_best_experiments is absent without saying why")
        if not contract.get("audit_id"):
            w12_bad.append("audit_id is empty")
        for c in contract.get("candidates", []):
            missing = [k for k in ("candidate_id", "gate_status", "scores",
                                   "uncertainty", "decision") if k not in c]
            if missing:
                w12_bad.append(f"{c.get('candidate_id')} lacks {missing}")
            zeroed = [k for k, v in (c.get("scores") or {}).items() if v == 0.0]
            if zeroed:
                w12_bad.append(f"{c.get('candidate_id')} reports {zeroed} as "
                               "0.00; an unmeasured component must be null")
    criterion("W12", bool(w12_bad),
              f"the contract carries all {len(required)} named fields for "
              f"{len(contract.get('candidates', []))} candidate(s) and "
              f"{contract.get('excluded_candidate_count')} excluded; "
              f"next_best_experiments is null with a stated reason and "
              f"audit_id is {contract.get('audit_id')}"
              if not w12_bad else "; ".join(w12_bad[:3]))

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHAT THE API RETURNS FOR THIS INDICATION")
    print("=" * 72)
    print(f"    GET /constructs -> {constructs['status']}")
    for reason in constructs["reasons"]:
        print(f"      - {reason}")
    print()
    print(f"    GET /result     -> {result['status']}")
    for step in result["attrition"]:
        print(f"      {step['gate']:34s} -{step['dropped']:4d}  "
              f"{step['remaining']:4d} remain")
    return 0


if __name__ == "__main__":
    sys.exit(main())
