"""HTTP API over the pipeline. Standard library only.

**Job and poll, not request and response.** A screen takes minutes: it reads a
9 GB single-cell matrix, evaluates 19,900 pairs and makes a network call per pool
member. A synchronous endpoint would time out in every client, so a run is
submitted, a job identifier comes back, and the client polls.

**A zero-construct pipeline is a result.** Every collection endpoint returns a
`status` and a `reasons` list beside its rows. A run that assembles nothing
answers 200 with `status: NO_BUILDABLE_CONSTRUCT` and the reasons it did not —
never 404, never 500, and never a bare empty list. An empty list reads as "we
looked and there was nothing to say"; the truth is that something specific and
measured stopped each design, and the caller needs that rather than the absence
of it.

Deployment shape: one process, an in-memory job table, a thread per run.
Restarting loses jobs and loses nothing else, because every stage's real output
is on disk under its own manifest.
"""

from __future__ import annotations

import json
import re
import threading
import traceback
import uuid
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from car_pipeline.api import pipeline
from car_pipeline.data.source import CacheError  # noqa: F401
from car_pipeline.stages import stage6, stage10, stage11

_LOCK = threading.Lock()
PROJECTS: dict[str, dict] = {}
JOBS: dict[str, dict] = {}
RESULTS: dict[str, dict] = {}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --------------------------------------------------------------------------
# projects and jobs
# --------------------------------------------------------------------------


def create_project(cancer_type: str) -> dict:
    if not isinstance(cancer_type, str) or not cancer_type.strip():
        raise ValueError("cancer_type is required and must be a string")
    # Refused here rather than at run time, so a caller learns immediately that
    # only some indications are configured instead of receiving another
    # indication's results under their own name.
    pipeline.project_for(cancer_type)
    project_id = uuid.uuid4().hex[:12]
    project = {
        "project_id": project_id,
        "cancer_type": cancer_type.strip(),
        # The null is what selects discovery mode, and it stays null. Seeding it
        # would void the exercise.
        "target_antigen": None,
        "discovery_mode": "B",
        "created_at": _now(),
    }
    with _LOCK:
        PROJECTS[project_id] = project
    return project


def start_run(project_id: str) -> dict:
    with _LOCK:
        if project_id not in PROJECTS:
            raise KeyError(project_id)
        # One at a time. Two runs would race on RESULTS and on the shared binder
        # cache, whose writer unlinks the manifest before rewriting the payload —
        # a reader between those two steps sees an unblessed artifact.
        running = [j for j in JOBS.values()
                   if j["project_id"] == project_id
                   and j["status"] in ("queued", "running")]
        if running:
            raise RuntimeError(
                f"run {running[0]['job_id']} is already {running[0]['status']} "
                "for this project"
            )
        job_id = uuid.uuid4().hex[:12]
        job = {
            "job_id": job_id, "project_id": project_id, "status": "queued",
            "stage": None, "note": "", "started_at": _now(),
            "finished_at": None, "error": None,
            "stages": list(pipeline.STAGES),
        }
        JOBS[job_id] = job

    def worker():
        def progress(stage, note=""):
            with _LOCK:
                job["status"] = "running"
                job["stage"] = stage
                job["note"] = note
        try:
            result = pipeline.run(PROJECTS[project_id]["cancer_type"], progress)
            with _LOCK:
                RESULTS[project_id] = result
                job["status"] = "complete"
                job["stage"] = "ranking"
                job["note"] = result["status"]
                job["finished_at"] = _now()
        except Exception as exc:                       # noqa: BLE001
            with _LOCK:
                job["status"] = "failed"
                job["error"] = f"{type(exc).__name__}: {exc}"
                job["trace"] = traceback.format_exc()[-2000:]
                job["finished_at"] = _now()

    threading.Thread(target=worker, daemon=True).start()
    return job


# --------------------------------------------------------------------------
# views
# --------------------------------------------------------------------------


class RunNotComplete(LookupError):
    """No completed run for this project. A statement about the job."""


class NotFound(LookupError):
    """The run completed and does not contain what was asked for."""


def _result(project_id: str) -> dict:
    result = RESULTS.get(project_id)
    if result is None:
        raise RunNotComplete(project_id)
    return result


def targets_view(project_id: str, limit: int = 50) -> dict:
    r = _result(project_id)
    # Sorted here. stage3.rank() emits in universe order, not ranked order, and
    # an endpoint called /targets returning that would be unranked rows under a
    # ranked name. Same key the pool uses, so the two orders agree.
    scored = sorted(
        (x for x in r["ranked"] if x.composite is not None and x.gene),
        key=lambda x: (-x.composite, x.gene, x.accession),
    )
    rows = []
    for rank_index, entry in enumerate(scored[:limit], 1):
        rows.append({
            "rank": rank_index,
            "gene": entry.gene,
            "accession": entry.accession,
            "evidence_class": entry.evidence_class,
            "tier_rank": entry.tier_rank,
            "composite": entry.composite,
            "measured_weight": entry.measured_weight,
            "risk": entry.risk,
            "risk_organ": entry.risk_organ,
            "cleared": entry.cleared,
            "confidence": entry.confidence,
            "breakdown": {
                key: entry.component_value(key)
                for key in entry.components
            },
        })
    return {
        "status": "RANKED",
        "universe": len(r["ranked"]),
        "scored": len(scored),
        "ceiling": r["ceiling"],
        "returned": len(rows),
        "targets": rows,
        "reasons": [],
    }


def pairs_view(project_id: str, limit: int = 50) -> dict:
    r = _result(project_id)
    measured = [p for p in r["pairs"] if p.coverage.measured]
    rows = []
    # `or 9` would treat a combined risk of exactly 0.0 — the safest pair
    # reachable, since a per-organ score of 0.0 is a real measurement — as the
    # worst and sort it off the page.
    ordered = sorted(
        measured,
        key=lambda x: (x.risk.combined is None, x.risk.combined or 0.0),
    )
    for p in ordered[:limit]:
        rows.append({
            "gene_a": p.gene_a, "gene_b": p.gene_b,
            "combined_risk": p.risk.combined,
            "peak_organ": p.risk.organ,
            "cleared": p.cleared,
            "coverage_f_ab": p.coverage.f_ab,
            "coverage_span_kb": p.coverage.span_geomean_kb,
            "coverage_span_percentile": p.coverage.span_percentile,
            # Carried on the number, not in a footnote.
            "coverage_caveat": (
                "span-confounded; read the percentile beside the fraction"),
        })
    return {
        "status": "PAIRED",
        "evaluated": len(r["pairs"]),
        "measured": len(measured),
        "returned": len(rows),
        "pairs": rows,
        "reasons": [
            "Coverage does not gate: per-cell detection tracks genomic span "
            "(rho +0.68) more strongly than expression (+0.20).",
            "Stage 4 closed at 10 of 14 criteria with four documented limitations.",
        ],
    }


def constructs_view(project_id: str) -> dict:
    """Zero buildable constructs is a result, with reasons, at HTTP 200."""
    r = _result(project_id)
    constructs = r["constructs"]
    counts: dict[str, int] = {}
    for c in constructs:
        counts[c.verdict] = counts.get(c.verdict, 0) + 1
    buildable = [c for c in constructs if c.verdict == stage6.BUILDABLE]
    assembled = [c for c in constructs if c.amino_acid_sequence]

    rows = []
    for c in assembled:
        rows.append({
            "gene": c.gene, "partner": c.partner, "verdict": c.verdict,
            "architecture": c.architecture,
            "binder": c.binder_name, "partner_binder": c.partner_binder_name,
            "total_bp": c.total_bp, "budget_bp": stage6.BUDGET_BP,
            "headroom_bp": c.headroom_bp,
            "amino_acid_sequence": c.amino_acid_sequence,
            "dna": c.dna,
            "domains": [
                {"name": s.name, "provenance": s.provenance,
                 "accession": s.accession, "feature": s.feature,
                 "source_residues": (
                     f"{s.start_residue}-{s.end_residue}"
                     if s.start_residue else None),
                 "aa_start": s.aa_start, "aa_end": s.aa_end,
                 "bp_start": s.bp_start, "bp_end": s.bp_end}
                for s in c.segments
            ],
            "reason": c.reason,
        })

    # Reasons are computed from this run. Prose with a hardcoded fallback number
    # would present a remembered figure as a measurement.
    reasons = []
    if not buildable:
        switch_bp = sum(
            seg.bp_end - seg.bp_start
            for seg in (assembled[0].segments if assembled else [])
            if "fkbp" in seg.name.lower() or "caspase" in seg.name.lower()
            or seg.name == "SGGGS linker"
        )
        candidates = sum(len(b.structure) + len(b.sequence)
                         for b in r["binders"].values())
        single_domain = sum(
            1 for b in r["binders"].values()
            for c in (b.structure + b.sequence) if "single" in c.fmt.lower()
        )
        reasons.append(
            "Conservative safety tolerance mandates a safety switch"
            + (f" ({switch_bp} bp)." if switch_bp else ".")
        )
        if assembled:
            worst = max(assembled, key=lambda c: c.total_bp)
            reasons.append(
                f"The largest assembled design reaches {worst.total_bp} bp "
                f"against a {stage6.BUDGET_BP} bp payload budget, over by "
                f"{-worst.headroom_bp}."
            )
        else:
            reasons.append("No design assembled at all.")
        reasons.append(
            f"Single-domain binders would fit; {single_domain} of {candidates} "
            "retrieved candidates are single-domain."
        )
        reasons.append(
            "This is a constraint result, not a pipeline failure. The budget is "
            "Stage 1's and is doing what it exists for."
        )
    return {
        "status": "NO_BUILDABLE_CONSTRUCT" if not buildable else "BUILDABLE",
        "counts": counts,
        "buildable": len(buildable),
        "assembled": len(assembled),
        "constructs": rows,
        "reasons": reasons,
    }


def result_view(project_id: str) -> dict:
    r = _result(project_id)
    running = len(r["final"])
    chain = []
    for gate in stage11.GATES:
        n = r["attrition"][gate]
        running -= n
        chain.append({"gate": gate, "dropped": n, "remaining": running})
    return {
        "status": r["status"],
        "pool_size": len(r["final"]),
        "reached_the_end": sum(1 for x in r["final"] if x.survived),
        "attrition": chain,
        "reasons": [
            "Every drop is a measurement, not a failure of the stage that made it.",
            "An empty ranking would read as 'nothing ranked highly'; the true "
            "statement is that nothing arrived to be ranked.",
        ],
        "developability_status": r["developability_status"],
    }


def evidence_view(project_id: str, gene: str) -> dict:
    """Everything the pipeline holds about one gene, stage by stage."""
    r = _result(project_id)
    ranked = next((x for x in r["ranked"] if x.gene == gene), None)
    if ranked is None:
        raise NotFound(gene)
    decision = next((d for d in r["decisions"] if d["gene"] == gene), None)
    binder = r["binders"].get(gene)
    construct = next((c for c in r["constructs"] if c.gene == gene), None)
    gate = next((g for g in r["gated"] if g.gene == gene), None)
    dev = [d for d in r["developability"] if d.gene == gene]
    final = next((f for f in r["final"] if f.gene == gene), None)

    return {
        "status": "EVIDENCE",
        "gene": gene,
        "accession": ranked.accession,
        "stage3_screen": {
            "evidence_class": ranked.evidence_class,
            "composite": ranked.composite,
            "risk": ranked.risk, "risk_organ": ranked.risk_organ,
            "cleared": ranked.cleared, "confidence": ranked.confidence,
            "components": {k: ranked.component_value(k) for k in ranked.components},
        },
        "stage4_pairing": decision,
        "stage5_binders": None if binder is None else {
            "verdict": binder.verdict,
            "entries": len(binder.entries),
            "structure_candidates": len(binder.structure),
            "sequence_candidates": [
                {"name": c.name, "clinical_stage": c.clinical_stage,
                 "status": c.status, "format": c.fmt,
                 "affinity": c.affinity, "isoform": c.isoform}
                for c in binder.sequence
            ],
        },
        "stage6_construct": None if construct is None else {
            "verdict": construct.verdict, "reason": construct.reason,
            "total_bp": construct.total_bp,
        },
        "stage9_safety": None if gate is None else {
            "verdict": gate.verdict, "reasons": gate.reasons,
            "binder_origins": gate.binder_origins,
            "epitope_immunogenicity": gate.epitope_immunogenicity,
            "trials_total": gate.trials_total,
            "trials_stopped": gate.trials_stopped,
            "trials_truncated": gate.trials_truncated,
        },
        "stage10_developability": [
            {"binder": d.binder, "isoelectric_point": d.isoelectric_point,
             "net_charge": d.net_charge, "cysteine_parity": d.cysteine_parity,
             "flags": [{"kind": k, "detail": v} for k, v in d.flags]}
            for d in dev
        ],
        "stage11_outcome": None if final is None else {
            "survived": final.survived, "failed_at": final.failed_at,
        },
        "reasons": [],
    }


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    server_version = "car-platform/1"

    def _send(self, code: int, payload: dict) -> None:
        body = json.dumps(payload, indent=2, default=str).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):        # quieter default logging
        return

    def _body(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        if not length:
            return {}
        return json.loads(self.rfile.read(length) or b"{}")

    def do_POST(self):                        # noqa: N802
        path = self.path.split("?")[0]
        try:
            if path == "/projects":
                body = self._body()
                return self._send(201, create_project(body.get("cancer_type", "")))
            m = re.match(r"^/projects/([0-9a-f]{12})/runs$", path)
            if m:
                return self._send(202, start_run(m.group(1)))
        except ValueError as exc:
            return self._send(400, {"status": "BAD_REQUEST", "error": str(exc)})
        except RuntimeError as exc:
            return self._send(409, {"status": "RUN_IN_PROGRESS", "error": str(exc)})
        except KeyError as exc:
            return self._send(404, {"status": "NOT_FOUND", "error": str(exc)})
        except Exception as exc:                       # noqa: BLE001
            # Nothing leaves this handler without an HTTP response. A closed
            # connection is the one failure a client cannot interpret.
            return self._send(500, {
                "status": "INTERNAL_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            })
        return self._send(404, {"status": "NOT_FOUND", "path": path})

    def _limit(self, default: int = 50) -> int:
        """`?limit=` was accepted in the signature and never passed. Bounded."""
        query = self.path.split("?", 1)[1] if "?" in self.path else ""
        for part in query.split("&"):
            if part.startswith("limit="):
                try:
                    return max(1, min(500, int(part[6:])))
                except ValueError:
                    return default
        return default

    def do_GET(self):                         # noqa: N802
        path = self.path.split("?")[0]
        try:
            m = re.match(r"^/jobs/([0-9a-f]{12})$", path)
            if m:
                with _LOCK:
                    job = JOBS.get(m.group(1))
                    # Copied inside the lock: the worker adds a "trace" key on
                    # failure, and serialising the live dict can raise
                    # "dictionary changed size during iteration".
                    snapshot = dict(job) if job else None
                if snapshot is None:
                    return self._send(404, {"status": "NOT_FOUND"})
                return self._send(200, snapshot)

            for name, view, paged in (("targets", targets_view, True),
                                      ("pairs", pairs_view, True),
                                      ("constructs", constructs_view, False),
                                      ("result", result_view, False)):
                m = re.match(rf"^/projects/([0-9a-f]{{12}})/{name}$", path)
                if m:
                    if paged:
                        return self._send(200, view(m.group(1), self._limit()))
                    return self._send(200, view(m.group(1)))

            m = re.match(r"^/projects/([0-9a-f]{12})/evidence/([A-Za-z0-9_.-]+)$", path)
            if m:
                return self._send(200, evidence_view(m.group(1), m.group(2)))
        except RunNotComplete:
            # Not an error about the data. A statement about the job, naming the
            # call that would fix it.
            return self._send(409, {
                "status": "RUN_NOT_COMPLETE",
                "reasons": ["No completed run for this project. POST "
                            "/projects/{id}/runs, then poll /jobs/{job_id}."],
            })
        except NotFound as exc:
            # The run completed and does not contain this. Telling the caller to
            # re-run would send them round a loop that cannot help.
            return self._send(404, {
                "status": "NOT_FOUND", "error": str(exc),
                "reasons": ["The run completed; this identifier is not in it."],
            })
        except Exception as exc:                       # noqa: BLE001
            return self._send(500, {
                "status": "INTERNAL_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            })
        return self._send(404, {"status": "NOT_FOUND", "path": path})


def serve(host: str = "127.0.0.1", port: int = 8000) -> None:
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"  listening on http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    serve()
