"""HTTP API over the pipeline. Standard library only."""

from __future__ import annotations

import json
import re
import threading
import traceback
import uuid
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from car_pipeline.api import pipeline
from car_pipeline.data.source import CacheError
from car_pipeline.stages import stage4, stage6, stage10, stage11, validation

_LOCK = threading.Lock()
PROJECTS: dict[str, dict] = {}
JOBS: dict[str, dict] = {}
RESULTS: dict[str, dict] = {}


MAX_RESULTS = 8


def _evict_results() -> None:
    """Drop the oldest completed runs beyond MAX_RESULTS. Call under _LOCK."""
    while len(RESULTS) > MAX_RESULTS:
        RESULTS.pop(next(iter(RESULTS)))


def _now() -> str:
    """The current UTC time, to the second."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def create_project(cancer_type: str, target_antigen: str | None = None) -> dict:
    """Register a project and choose the discovery mode from the target."""
    if not isinstance(cancer_type, str) or not cancer_type.strip():
        raise ValueError("cancer_type is required and must be a string")
    if target_antigen is not None and not str(target_antigen).strip():
        raise ValueError(
            "target_antigen was supplied but is empty; omit it entirely to "
            "screen for targets")

    pipeline.project_for(cancer_type)
    project_id = uuid.uuid4().hex[:12]
    project = {
        "project_id": project_id,
        "cancer_type": cancer_type.strip(),
        "target_antigen": (target_antigen or "").strip().upper() or None,
        "discovery_mode": "A" if target_antigen else "B",
        "created_at": _now(),
    }
    with _LOCK:
        PROJECTS[project_id] = project
    return project


def start_run(project_id: str) -> dict:
    """Queue a screen for this project unless one is already running."""
    with _LOCK:
        if project_id not in PROJECTS:
            raise KeyError(project_id)

        running = [j for j in JOBS.values()
                   if j["status"] in ("queued", "running")]
        if running:
            other = running[0]
            same = other["project_id"] == project_id
            raise RuntimeError(
                f"run {other['job_id']} is already {other['status']}"
                + (" for this project" if same else
                   f" for project {other['project_id']}; the pipeline writes a "
                   "shared cache and runs one at a time")
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
        """Run the screen on a background thread and record the outcome."""
        def progress(stage, note=""):
            """Record which stage the run has reached."""
            with _LOCK:
                job["status"] = "running"
                job["stage"] = stage
                job["note"] = note
        try:
            project = PROJECTS[project_id]
            target = project.get("target_antigen")
            if target:
                validation = pipeline.validate(
                    project["cancer_type"], target, progress)
                result = validation["screen"]
                result["validation"] = validation
            else:
                result = pipeline.run(project["cancer_type"], progress)
            with _LOCK:
                RESULTS.pop(project_id, None)
                RESULTS[project_id] = result
                _evict_results()
                job["status"] = "complete"
                job["stage"] = "ranking"
                job["note"] = result["status"]
                job["finished_at"] = _now()
        except Exception as exc:
            with _LOCK:
                job["status"] = "failed"
                job["error"] = f"{type(exc).__name__}: {exc}"
                job["trace"] = traceback.format_exc()[-2000:]
                job["finished_at"] = _now()

    threading.Thread(target=worker, daemon=True).start()
    return job


class RunNotComplete(LookupError):
    """No completed run for this project. A statement about the job."""


class UnknownProject(LookupError):
    """No project with this id. A statement about the id, not the job."""


class NotFound(LookupError):
    """The run completed and does not contain what was asked for."""


def _result(project_id: str) -> dict:
    """The finished result for a project, or a refusal naming which it is."""
    if project_id not in PROJECTS:
        raise UnknownProject(project_id)
    result = RESULTS.get(project_id)
    if result is None:
        raise RunNotComplete(project_id)
    return result


def _evidence(r: dict) -> dict:
    """The run-level evidence judgement, attached to every view."""
    unavailable = r.get("unavailable") or []
    out = {
        "usability": r.get("usability", pipeline.USABLE),
        "unavailable": unavailable,
        "indication": (r.get("indication").cancer_type
                       if r.get("indication") else None),
    }
    if unavailable:
        out["degraded_note"] = (
            f"{len(unavailable)} component(s) could not be measured for this "
            "indication; each names the source it needed. The ranking is still "
            "supported, but it rests on less evidence than a complete run."
        )
    return out


def targets_view(project_id: str, limit: int = 50) -> dict:
    """The ranked targets, or the refusal when the evidence cannot support a ranking."""
    r = _result(project_id)
    if r.get("usability") == pipeline.NOT_USABLE:
        return {
            "status": pipeline.NOT_USABLE,
            **_evidence(r),
            "targets": [],
            "reasons": r.get("reasons", []),
        }

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
        **_evidence(r),
        "universe": len(r["ranked"]),
        "scored": len(scored),
        "ceiling": r["ceiling"],
        "returned": len(rows),
        "targets": rows,
        "reasons": [],
    }


def pairs_view(project_id: str, limit: int = 50) -> dict:
    """The admissible pairs, ordered by how far under the ceiling they sit."""
    r = _result(project_id)
    measured = [p for p in r["pairs"] if p.coverage.measured]
    rows = []

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
            "coverage_caveat": (
                "span-confounded; read the percentile beside the fraction"),
        })
    return {
        "status": "PAIRED",
        **_evidence(r),
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
    complete = [c for c in buildable if c.amino_acid_sequence]
    awaiting = [c for c in buildable if not c.amino_acid_sequence]
    over_budget = [c for c in constructs if c.verdict == stage6.BUDGET_EXCEEDED]

    def state(c):
        """How far a construct got: complete, awaiting a binder, or over budget."""
        if c.verdict != stage6.BUILDABLE:
            return "BUDGET_EXCEEDED"
        return "COMPLETE" if c.amino_acid_sequence else "AWAITING_BINDER"

    rows = []
    for c in buildable + over_budget:
        rows.append({
            "gene": c.gene, "partner": c.partner, "verdict": c.verdict,
            "state": state(c),
            "design_class": validation.design_class(c),
            "architecture": c.architecture,
            "binder_supplied": c.binder_supplied,
            "binder": c.binder_name, "partner_binder": c.partner_binder_name,
            "total_bp": c.total_bp, "budget_bp": stage6.BUDGET_BP,
            "headroom_bp": c.headroom_bp,
            "amino_acid_sequence": c.amino_acid_sequence or None,
            "dna": c.dna or None,
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
    if not buildable:
        status = "NO_BUILDABLE_CONSTRUCT"
    elif complete:
        status = "BUILDABLE"
    else:
        status = "BUILDABLE_AWAITING_BINDER"

    if buildable and not complete:
        reasons.append(
            f"{len(awaiting)} design(s) fit the {stage6.BUDGET_BP} bp budget "
            "but carry no binder sequence: the adaptor receptor binds a tag, "
            "and no anti-tag binder exists in the connected sources, so its "
            "size is declared and its sequence is not invented."
        )
    reasons.extend(adaptor_notices(r["constructs"]))
    classes = validation.design_class_summary(constructs)
    reasons.extend(classes["reasons"])
    return {
        "status": status,
        **_evidence(r),
        "design_classes": {k: classes[k]
                           for k in ("conservative_backup", "advanced")},
        "counts": counts,
        "buildable": len(buildable),
        "complete": len(complete),
        "awaiting_binder": len(awaiting),
        "over_budget": len(over_budget),
        "constructs": rows,
        "reasons": reasons,
    }


def validation_view(project_id: str) -> dict:
    """Mode A. What the platform concluded about a supplied target."""
    r = _result(project_id)
    v = r.get("validation")
    if v is None:
        return {
            "status": "NOT_APPLICABLE",
            "reasons": ["This project supplied no target, so it ran in "
                        "discovery mode. Create a project with target_antigen "
                        "to ask whether a specific target is suitable."],
        }
    return {
        "status": v["verdict"],
        "mode": "A",
        "cancer_type": v["cancer_type"],
        "target": v["target"],
        "accession": v.get("accession"),
        "rank": v.get("rank"),
        "of": v.get("of"),
        "composite": v.get("composite"),
        "measured_weight": v.get("measured_weight"),
        "risk": v.get("risk"),
        "risk_organ": v.get("risk_organ"),
        "evidence_class": v.get("evidence_class"),
        "architecture": v.get("architecture"),
        "reasons": v["reasons"],
    }


def adaptor_notices(constructs) -> list[str]:
    """What a reader must be told about a surviving adaptor design."""
    from car_pipeline.data.antitag import (
        DEPOSITION_ARTIFACTS, IDENTIFICATION, SPECIES_NOTICE)

    built = [c for c in constructs
             if c.verdict == stage6.BUILDABLE and c.outcome == stage4.ADAPTOR]
    if not built:
        return []
    notices = [TWO_PRODUCTS, IDENTIFICATION, SPECIES_NOTICE]

    sample = built[0]
    found = []
    for motif, what in DEPOSITION_ARTIFACTS:
        at = (sample.amino_acid_sequence or "").find(motif)
        if at >= 0:
            found.append(f"{motif} at residues {at + 1}-{at + len(motif)}, {what}")
    if found:
        notices.append(
            "The binder is emitted as deposited, including its crystallisation "
            "artifacts, because trimming them is a design decision this "
            "pipeline does not take silently. Each construct therefore carries "
            + "; and ".join(found)
            + ". As emitted these are not manufacturable: the first is a second "
            "leader sitting inside the mature protein, the second a His tag "
            "between the binder and the hinge. Removing them is a wet-lab step "
            "that has not been taken here."
        )
    return notices


TWO_PRODUCTS = (
    "Every surviving design routes to an adaptor architecture. That is two manufactured biologics, not one: the receptor and, separately, the tagged adaptor antibody that gives it its specificity. The second carries its own CMC package and its own regulatory path, and the payload budget the adaptor route saves is paid there instead."
)


def plan_view(project_id: str, gene: str) -> dict:
    """The validation plan for one candidate construct."""
    r = _result(project_id)
    symbol = gene.strip().upper()
    construct = next(
        (c for c in r["constructs"]
         if c.gene == symbol and c.verdict == stage6.BUILDABLE), None)
    if construct is None:
        assembled = sorted(c.gene for c in r["constructs"]
                           if c.verdict == stage6.BUILDABLE)
        raise NotFound(
            symbol if any(c.gene == symbol for c in r["constructs"])
            else symbol)
    ranked = next((x for x in r["ranked"] if x.gene == symbol), None)
    safety = next((s for s in r["gated"] if s.gene == symbol), None)
    target = {
        "risk": getattr(ranked, "risk", None),
        "risk_organ": getattr(ranked, "risk_organ", None),
        "evidence_class": getattr(ranked, "evidence_class", None),
    }
    return {**_evidence(r), **validation.plan(construct, target, safety)}


def result_view(project_id: str) -> dict:
    """The whole run: attrition chain, end states and recommendations."""
    r = _result(project_id)
    running = len(r["final"])
    chain = []
    for gate in stage11.GATES:
        n = r["attrition"][gate]
        running -= n
        chain.append({"gate": gate, "dropped": n, "remaining": running})
    survivors = [x for x in r["final"] if x.survived]
    complete = [x for x in survivors if x.binder_supplied]
    awaiting = [x for x in survivors if not x.binder_supplied]

    if not survivors:
        reasons = [
            "Every drop is a measurement, not a failure of the stage that made it.",
            "An empty ranking would read as 'nothing ranked highly'; the true "
            "statement is that nothing arrived to be ranked.",
        ]
    else:
        reasons = [
            f"{len(survivors)} design(s) reached the end: {len(complete)} "
            f"complete, {len(awaiting)} awaiting a binder sequence.",
            "Every drop is a measurement, not a failure of the stage that "
            "made it.",
        ]
        if awaiting:
            reasons.append(
                "A design awaiting a binder has a layout, a length and a "
                "domain map, and no residues for its anti-tag binder. It fits "
                "the budget and cannot be ordered yet, which is a different "
                "state both from a design that does not fit and from a "
                "finished one."
            )
    return {
        "status": r["status"],
        **_evidence(r),
        "pool_size": len(r["final"]),
        "reached_the_end": len(survivors),
        "complete": len(complete),
        "awaiting_binder": len(awaiting),
        "attrition": chain,
        "reasons": reasons + adaptor_notices(r["constructs"]),
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


class Handler(BaseHTTPRequestHandler):
    server_version = "car-platform/1"

    def _send(self, code: int, payload: dict) -> None:
        """Write one JSON response."""
        body = json.dumps(payload, indent=2, default=str).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        """Silence the default request logging."""
        return

    def _body(self) -> dict:
        """The decoded JSON request body, or an empty mapping."""
        length = int(self.headers.get("Content-Length") or 0)
        if not length:
            return {}
        return json.loads(self.rfile.read(length) or b"{}")

    def do_POST(self):
        """Route the write endpoints."""
        path = self.path.split("?")[0]
        try:
            if path == "/projects":
                body = self._body()
                return self._send(201, create_project(
                        body.get("cancer_type", ""),
                        body.get("target_antigen")))
            m = re.match(r"^/projects/([0-9a-f]{12})/runs$", path)
            if m:
                return self._send(202, start_run(m.group(1)))
        except ValueError as exc:
            return self._send(400, {"status": "BAD_REQUEST", "error": str(exc)})
        except RuntimeError as exc:
            return self._send(409, {"status": "RUN_IN_PROGRESS", "error": str(exc)})
        except KeyError as exc:
            return self._send(404, {"status": "NOT_FOUND", "error": str(exc)})
        except Exception as exc:
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

    def do_GET(self):
        """Route the read endpoints."""
        path = self.path.split("?")[0]
        try:
            if path == "/indications":
                from car_pipeline.configs.registry import registered
                return self._send(200, {
                    "status": "CONFIGURED",
                    "indications": registered(),
                    "reasons": [
                        "An indication needs a tumour cohort, a single-cell "
                        "atlas, a dependency lineage and a normal-tissue "
                        "denominator declared before it can be screened. None "
                        "of those is derivable from the cancer type, so an "
                        "unregistered one is refused rather than answered with "
                        "another indication's results.",
                    ],
                })

            m = re.match(r"^/jobs/([0-9a-f]{12})$", path)
            if m:
                with _LOCK:
                    job = JOBS.get(m.group(1))

                    snapshot = dict(job) if job else None
                if snapshot is None:
                    return self._send(404, {
                        "status": "NOT_FOUND",
                        "error": m.group(1),
                        "reasons": ["No job with this id. Jobs live in memory "
                                    "and do not survive a restart."],
                    })
                return self._send(200, snapshot)

            m = re.match(r"^/projects/([0-9a-f]{12})/validation$", path)
            if m:
                return self._send(200, validation_view(m.group(1)))
            for name, view, paged in (("targets", targets_view, True),
                                      ("pairs", pairs_view, True),
                                      ("constructs", constructs_view, False),
                                      ("result", result_view, False)):
                m = re.match(rf"^/projects/([0-9a-f]{{12}})/{name}$", path)
                if m:
                    if paged:
                        return self._send(200, view(m.group(1), self._limit()))
                    return self._send(200, view(m.group(1)))

            m = re.match(r"^/projects/([0-9a-f]{12})/plan/([A-Za-z0-9_.-]+)$", path)
            if m:
                return self._send(200, plan_view(m.group(1), m.group(2)))

            m = re.match(r"^/projects/([0-9a-f]{12})/evidence/([A-Za-z0-9_.-]+)$", path)
            if m:
                return self._send(200, evidence_view(m.group(1), m.group(2)))
        except UnknownProject as exc:
            return self._send(404, {
                "status": "NOT_FOUND",
                "error": str(exc),
                "reasons": ["No project with this id. It was never created, or "
                            "the service restarted: projects live in memory "
                            "and do not survive one. POST /projects to create "
                            "one. This is not a run that has yet to finish."],
            })
        except RunNotComplete:
            return self._send(409, {
                "status": "RUN_NOT_COMPLETE",
                "reasons": ["This project exists and has no completed run. "
                            "POST /projects/{id}/runs, then poll "
                            "/jobs/{job_id}."],
            })
        except NotFound as exc:
            return self._send(404, {
                "status": "NOT_FOUND", "error": str(exc),
                "reasons": ["The run completed; this identifier is not in it."],
            })
        except Exception as exc:
            return self._send(500, {
                "status": "INTERNAL_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            })
        return self._send(404, {"status": "NOT_FOUND", "path": path})


def serve(host: str = "127.0.0.1", port: int = 8000) -> None:
    """Start the HTTP server and block."""
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"  listening on http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Serve the design platform.")

    parser.add_argument("--host", default=os.environ.get("HOST") or "127.0.0.1")

    parser.add_argument("--port", type=int, default=None)
    options = parser.parse_args()

    port = options.port
    if port is None:
        raw = (os.environ.get("PORT") or "").strip()
        if not raw:
            port = 8000
        elif raw.isdigit() and 0 < int(raw) < 65536:
            port = int(raw)
        else:
            parser.error(f"PORT={raw!r} in the environment is not a port "
                         "number; pass --port instead")
    serve(options.host, port)
