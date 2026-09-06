"""The /v1/cart surface, translated onto the pipeline's own.

**This module translates and never decides.** It renames fields, re-routes
paths, and fans a batch request out over a per-item view. It applies no gate,
computes no ranking, and carries no vocabulary of its own. Anything it cannot do
by translation is a gap in the pipeline and is built there.

Stated as an intention that would be worth nothing, so it is enforced:

* **D1** parses this module's own import statements and trips if any names a
  stage module, the scoring frame, routing or validation. Note the check is on
  *direct* imports rather than the transitive closure: `server` legitimately
  imports every stage, and reaching those stages *through* its published views
  is exactly the intended design. What must not happen is this module reaching
  them directly and forming its own opinion.
* **D2** scans this source for any value from the decision, gate-status or
  design-class vocabularies. If this module needed to spell a decision it would
  be deciding; passing one through never requires naming it.

Both are blinded before they are trusted.
"""

from __future__ import annotations

from car_pipeline.api import constants, server

PREFIX = "/v1/cart"

# Job statuses are renamed, not interpreted: a lookup with no branch on meaning.
RUN_STATUS = {"queued": "QUEUED", "running": "RUNNING",
              "complete": "COMPLETE", "failed": "FAILED"}

# Logic modes the pair request may ask for. The refused one has no implemented
# architecture behind it, and accepting it then returning something else would
# be a silent substitution.
LOGIC_SUPPORTED = {"AND"}
LOGIC_REFUSED = {
    "OR": "no OR-gate architecture is implemented; routing has no such row",
    "AND-NOT": "the inhibitory architecture is recorded as not implemented, "
               "with its reason, so it cannot be selected",
}


class Refusal(Exception):
    """A translation that cannot be made, named rather than guessed at."""

    def __init__(self, status: int, payload: dict):
        super().__init__(payload.get("error", "refused"))
        self.status = status
        self.payload = payload


def _refuse(field: str, error: str, **extra) -> Refusal:
    """A 400 naming the field, in the same shape the pipeline already uses."""
    return Refusal(400, {"status": "UNSUPPORTED_INPUT", "field": field,
                         "error": error, "remove": [field], **extra})


def _project(body: dict) -> str:
    """The canonical project id from whichever name the caller used."""
    reference = body.get("project_id") or body.get("run_id")
    if not reference:
        raise _refuse("project_id", "project_id is required")
    try:
        return server.resolve_project(str(reference).strip())
    except KeyError:
        raise Refusal(404, {"status": "NOT_FOUND", "error": str(reference),
                            "reasons": ["No project with this id or client "
                                        "reference."]})


def _limit(body: dict, *names: str, default: int = 50) -> int:
    """Whichever cap the caller named, bounded the way the pipeline bounds it."""
    for name in names:
        if body.get(name) is not None:
            try:
                return max(1, min(500, int(body[name])))
            except (TypeError, ValueError):
                raise _refuse(name, f"{name} must be a whole number")
    return default


# --------------------------------------------------------------------------
# The eight that translate
# --------------------------------------------------------------------------

def projects(body: dict) -> dict:
    """Create a project. Field names differ; the contract does not."""
    project = server.create_project(body)
    reference = project.get("client_reference")
    return {
        "project_id": reference or project["project_id"],
        "canonical_project_id": project["project_id"],
        "status": "CREATED",
        "indication": project["cancer_type"],
        "target_mode": project["target_mode"],
        "architecture_mode": project["architecture_mode"],
        "binder_mode": project["binder_mode"],
        "max_final_candidates": project["max_final_candidates"],
        "reasons": [
            "project_id echoes the reference you supplied; "
            "canonical_project_id is the server-generated id. Both address the "
            "same project, because ids must be unique and server-controlled "
            "while a client reference must stay yours.",
        ],
    }


def runs(body: dict) -> dict:
    """Start a run. The project moves from the body into the path."""
    job = server.start_run(_project(body))
    return {"run_id": job["job_id"],
            "project_id": body.get("project_id"),
            "status": RUN_STATUS.get(job["status"], job["status"].upper())}


def run_status(run_id: str) -> dict:
    """Poll a run. Progress is reported as stages, because that is what exists."""
    job = server.job_status(run_id)
    if not job:
        raise Refusal(404, {"status": "NOT_FOUND", "error": run_id,
                            "reasons": ["No run with this id. Runs live in "
                                        "memory and do not survive a restart."]})
    stages = list(job.get("stages") or [])
    stage = job.get("stage")
    reached = (stages.index(stage) + 1) if stage in stages else 0
    return {
        "run_id": run_id,
        "project_id": job.get("project_id"),
        "status": RUN_STATUS.get(job["status"], job["status"].upper()),
        "stage": stage,
        "stages_complete": reached,
        "stages_total": len(stages),
        "progress": None,
        "error": job.get("error"),
        "reasons": [
            "progress is null rather than a fraction. A fraction implies a "
            "measured proportion of the work done; what exists is an ordered "
            "list of stages and the one currently reached. stages_complete "
            "and stages_total carry that without inventing a proportion.",
        ],
    }


def evidence_build(body: dict) -> dict:
    """Normalised evidence per target. One view, fanned out over the list."""
    pid = _project(body)
    targets = body.get("target_ids") or []
    if not targets:
        raise _refuse("target_ids", "target_ids must name at least one target")
    out, missing = [], []
    for target in targets:
        try:
            out.append(server.evidence_view(pid, str(target)))
        except KeyError:
            missing.append(str(target))
    return {"evidence": out, "requested": len(targets),
            "returned": len(out), "not_found": missing,
            "reasons": ["A target with no evidence row is named in not_found "
                        "rather than omitted, so the counts reconcile."]}


def targets_discover(body: dict) -> dict:
    """Ranked targets."""
    return server.targets_view(_project(body), _limit(body, "max_targets"))


def targets_pair(body: dict) -> dict:
    """Evaluated pairs. A logic mode with no implemented row is refused."""
    for mode in body.get("logic_modes") or []:
        key = str(mode).strip().upper()
        if key in LOGIC_REFUSED:
            raise _refuse("logic_modes",
                          f"logic mode {key!r} is not implemented: "
                          f"{LOGIC_REFUSED[key]}",
                          supported_values=sorted(LOGIC_SUPPORTED),
                          rejected_values=sorted(LOGIC_REFUSED))
        if key not in LOGIC_SUPPORTED:
            raise _refuse("logic_modes",
                          f"logic mode {key!r} is not a value this platform "
                          "recognises",
                          supported_values=sorted(LOGIC_SUPPORTED),
                          rejected_values=sorted(LOGIC_REFUSED))
    return server.pairs_view(_project(body), _limit(body, "max_pairs"))


def gates_evaluate(body: dict) -> dict:
    """Gate status per candidate, selected from the run's own attribution."""
    contract = server.contract_view(_project(body))
    wanted = {str(c) for c in (body.get("candidate_ids") or [])}
    rows = contract["candidates"] + contract["excluded_candidates"]
    if wanted:
        rows = [c for c in rows
                if c.get("candidate_id") in wanted or c.get("gene") in wanted]
    return {
        "candidates": [
            {"candidate_id": c.get("candidate_id"),
             "target_id": c.get("gene") or (c.get("targets") or [None])[0],
             "gate_status": c.get("gate_status"),
             "failed_gates": ([] if c.get("candidate_id") and c in contract["candidates"]
                              else [c.get("gate_status")])}
            for c in rows
        ],
        "reasons": [
            "failed_gates carries at most one element. Every candidate is "
            "attributed to the first gate it failed and evaluation stops "
            "there, which is what makes the attrition chain sum to the pool. "
            "A longer list would imply every gate was evaluated on every "
            "candidate, and it was not.",
        ],
    }


def candidates_rank(body: dict) -> dict:
    """The final ranked list, which the pipeline already serves in this shape."""
    return server.contract_view(_project(body))


# --------------------------------------------------------------------------
# The four that needed a shape change, now served from pipeline views
# --------------------------------------------------------------------------

def binders_retrieve(body: dict) -> dict:
    """Retrieved binders for one target."""
    return server.binders_view(_project(body), body.get("target_id"))


def binder_scorecard(project_id: str, binder_id: str) -> dict:
    """One binder's record and provenance."""
    try:
        return server.binder_view(server.resolve_project(project_id), binder_id)
    except KeyError as exc:
        raise Refusal(404, {"status": "NOT_FOUND", "error": str(exc),
                            "reasons": ["No binder with this id in this run."]})


def constructs_generate(body: dict) -> dict:
    """The construct assembled for a target, selected from the completed run."""
    pid = _project(body)
    view = server.constructs_view(pid)
    target = (body.get("target_id") or "").strip().upper()
    rows = [c for c in view["constructs"]
            if not target or str(c.get("gene", "")).upper() == target]
    return {
        **{k: v for k, v in view.items() if k != "constructs"},
        "target_id": target or None,
        "constructs": rows,
        "reasons": list(view.get("reasons") or []) + [
            "Constructs are assembled during a run from the routed decision "
            "for each target, not built on demand from a supplied binder. "
            "This returns the construct the run produced for that target; a "
            "binder_id in the request does not select a different one.",
        ],
    }


def manufacturability_evaluate(body: dict) -> dict:
    """Payload and packaging, with the score's definition declared."""
    pid = _project(body)
    view = server.constructs_view(pid)
    contract = server.contract_view(pid)
    target = (body.get("target_id") or "").strip().upper()

    scores = {}
    for candidate in contract["candidates"]:
        card = candidate.get("scorecard") or {}
        for component in card.get("components") or []:
            if component.get("component") == "manufacturability":
                gene = (candidate.get("targets") or [None])[0]
                scores[gene] = component
    rows = []
    for row in view["constructs"]:
        gene = row.get("gene")
        if target and str(gene).upper() != target:
            continue
        component = scores.get(gene) or {}
        rows.append({
            "target_id": gene,
            "payload_bp": row.get("total_bp"),
            "payload_budget_bp": view.get("budget_bp"),
            "packaging_status": row.get("verdict"),
            "manufacturability_score": component.get("value"),
            "manufacturability_state": component.get("state"),
            "score_source": component.get("source"),
        })
    return {
        "constructs": rows,
        "score_definition": (
            "manufacturability is payload headroom as a fraction of the "
            "budget, measured after the hard packaging gate. It is not the "
            "same quantity as any illustrative figure in the interface "
            "specification, and it is reported under its own definition "
            "rather than rescaled to match one."),
        "reasons": list(view.get("reasons") or []),
    }


def experiments_recommend(body: dict) -> dict:
    """The validation plan per candidate, and the ranking that does not exist."""
    pid = _project(body)
    contract = server.contract_view(pid)
    plan = server.validation_view(pid)
    return {
        "experiments": None,
        "experiments_state": contract["next_best_experiments_state"],
        "validation_plan": plan,
        "reasons": [
            "experiments is null rather than an empty list. An empty list "
            "would say no experiment is recommended; the true statement is "
            "that nothing ranked one, because the stage that would has no "
            "producer. The per-candidate validation plan is carried instead, "
            "and it is a template of what to measure rather than a ranked "
            "recommendation of what is most informative.",
        ],
    }


def experiments_results(body: dict) -> dict:
    """Ingestion, which has no store behind it."""
    return {
        "status": "NOT_ACCEPTED",
        "calibration_ready": False,
        "ingested": 0,
        "reasons": [
            "No experimental result store exists, so nothing can be ingested "
            "and nothing is silently discarded. Returning ACCEPTED would tell "
            "a client their measurements were kept when they were not.",
            "This is the first link of the feedback path: until it exists, "
            "calibration stays blocked for want of data rather than for want "
            "of a model.",
        ],
    }


def binders_optimize(body: dict) -> dict:
    """Variant search, which has no producer."""
    return {
        "status": "NOT_IMPLEMENTED",
        "job_id": None,
        "reasons": [
            "Binder discovery retrieves and stops. Nothing generates a "
            "binder, models a complex, optimises an interface or germlines a "
            "framework, so there is no variant search to queue.",
            "A queued job id would promise a result that will never arrive.",
        ],
    }


def binders_rank(body: dict) -> dict:
    """Binder ranking, which is specified and not yet built."""
    return {
        "status": "NOT_IMPLEMENTED",
        "ranked_binders": None,
        "reasons": [
            "Binder ranking is specified and not yet built. It is null rather "
            "than an empty list, because an empty list would say no binder "
            "ranked and the true statement is that nothing ranked them.",
        ],
    }


# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------

_POST = {
    "/projects": projects,
    "/runs": runs,
    "/evidence/build": evidence_build,
    "/targets/discover": targets_discover,
    "/targets/pair": targets_pair,
    "/binders/retrieve": binders_retrieve,
    "/binders/optimize": binders_optimize,
    "/binders/rank": binders_rank,
    "/constructs/generate": constructs_generate,
    "/manufacturability/evaluate": manufacturability_evaluate,
    "/gates/evaluate": gates_evaluate,
    "/candidates/rank": candidates_rank,
    "/experiments/recommend": experiments_recommend,
    "/experiments/results": experiments_results,
}

_CONSTANTS = {
    "/structure/evaluate": lambda b: constants.structure_evaluate(
        b.get("construct_id")),
    "/function/predict": lambda b: constants.function_predict(
        b.get("construct_id")),
    "/learning/calibrate": lambda b: constants.learning_calibrate(
        b.get("dataset_version"), b.get("model_family")),
}

# Every path the interface specification documents. D9 asserts each answers.
PATHS = (tuple(f"{PREFIX}{p}" for p in _POST)
         + tuple(f"{PREFIX}{p}" for p in _CONSTANTS)
         + (f"{PREFIX}/runs/{{run_id}}", f"{PREFIX}/binders/{{binder_id}}"))


def dispatch_post(path: str, body: dict):
    """Handle one POST, or return None when the path is not ours."""
    if not path.startswith(PREFIX):
        return None
    tail = path[len(PREFIX):]
    if tail in _CONSTANTS:
        return 200, _CONSTANTS[tail](body)
    handler = _POST.get(tail)
    if handler is None:
        return None
    created = tail in ("/projects", "/runs")
    return (201 if tail == "/projects" else 202 if tail == "/runs" else 200,
            handler(body)) if created else (200, handler(body))


def dispatch_get(path: str, query: dict):
    """Handle one GET, or return None when the path is not ours."""
    if not path.startswith(PREFIX):
        return None
    tail = path[len(PREFIX):]
    if tail.startswith("/runs/"):
        return 200, run_status(tail[len("/runs/"):])
    if tail.startswith("/binders/"):
        binder = tail[len("/binders/"):]
        project = query.get("project_id")
        if not project:
            raise _refuse(
                "project_id",
                "a binder id is scoped to the run that retrieved it; supply "
                "?project_id= so the run is named")
        return 200, binder_scorecard(project, binder)
    return None
