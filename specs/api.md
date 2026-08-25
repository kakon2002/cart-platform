# The API

Standard library only. No framework, no new dependency.

```
.venv\Scripts\python.exe -m car_pipeline.api.server     # http://127.0.0.1:8000
.venv\Scripts\python.exe verify_api.py                  # end to end, 9 criteria
.venv\Scripts\python.exe run_all.py --fresh             # every stage, one report
```

**The interpreter must be the one in `.venv`.** A bare `py -3.13` resolves to
the system Python, which has none of the three dependencies and fails on the
first import. Nothing needs to be running first: the server loads the caches
under `data/` itself, and there is no database and no broker.

`--host` and `--port` override the defaults, falling back to `HOST` and `PORT`
from the environment. The default binds loopback only, so a demo is not also
serving the internet; a container has to be told `0.0.0.0` explicitly.

One process, an in-memory job table, a thread per run. Restarting loses jobs and
loses nothing else: every stage's real output is on disk under its own manifest.
One run per project at a time — two would race on the shared binder cache, whose
writer unlinks the manifest before rewriting the payload.

## Running it

**Nothing has to be running first.** No database, no broker, no queue, no
worker process. One command, and the caches under `data/` are read by the
server itself.

```
cd C:\cart-platform
.venv\Scripts\python.exe -m car_pipeline.api.server
```

`listening on http://127.0.0.1:8000`. A full screen takes about six minutes,
almost all of it Stage 5 making one network call per pool member, so the
machine needs to be online.

To demo it in one command, without driving the endpoints by hand:

```
.venv\Scripts\python.exe verify_api.py
```

That starts a server on 8137, creates a project, submits a run, polls the job
through each stage, checks nine criteria, and prints what the platform returns
for this indication. It is the whole flow end to end and it needs no arguments.

To drive it by hand against the server started above:

```
curl -s -X POST localhost:8000/projects -d "{\"cancer_type\":\"Pancreatic Ductal Adenocarcinoma\"}"
curl -s -X POST localhost:8000/projects/<id>/runs
curl -s localhost:8000/jobs/<job_id>            # poll until "complete"
curl -s localhost:8000/projects/<id>/constructs
```

## Job and poll, not request and response

A screen reads a 9 GB matrix, evaluates 19,900 pairs and makes a network call per
pool member. It takes minutes, so a synchronous endpoint would time out in every
client.

| | | |
| --- | --- | --- |
| `POST` | `/projects` | `{"cancer_type": "..."}` → `201` with a project id |
| `POST` | `/projects/{id}/runs` | → `202` with a job id |
| `GET` | `/jobs/{job_id}` | → `queued` / `running` + stage / `complete` / `failed` |
| `GET` | `/projects/{id}/targets?limit=` | ranked, with the six-component breakdown |
| `GET` | `/projects/{id}/pairs?limit=` | co-expression **and** its span percentile |
| `GET` | `/projects/{id}/constructs` | sequence, DNA map, domain provenance |
| `GET` | `/projects/{id}/result` | the end state and the attrition chain |
| `GET` | `/projects/{id}/evidence/{gene}` | every stage's view of one gene |

`target_antigen` is `null` at creation and stays `null`. The null is what selects
discovery mode.

## An empty pipeline is a result

**This is the part worth getting right.** For this indication the pipeline
produces no buildable construct, and that is a measurement, not a failure.

`GET /constructs` answers **`200`** with:

```json
{
  "status": "NO_BUILDABLE_CONSTRUCT",
  "counts": {"BUDGET_EXCEEDED": 2, "NO_CONSTRUCT": 198},
  "assembled": 2,
  "constructs": [ ... the two that assembled, with their overage ... ],
  "reasons": [
    "Conservative safety tolerance mandates a safety switch (1308 bp).",
    "The largest assembled design reaches 3894 bp against a 3500 bp payload budget, over by 394.",
    "Single-domain binders would fit; 1 of 735 retrieved candidates are single-domain.",
    "This is a constraint result, not a pipeline failure."
  ]
}
```

Never `404`, never `500`, never a bare `[]`. **An empty list reads as "we looked
and there was nothing to say."** Something specific and measured stopped each
design, and the caller needs that rather than the absence of it. The reasons are
computed from the run, not written down, so they cannot drift from it.

`GET /result` answers the same way with `NO_DESIGN_REACHES_THE_END` and the
attrition chain — every one of the 200 attributed to the first gate it failed,
summing to 200.

## Status codes, and what each means about whose problem it is

| code | status | meaning |
| --- | --- | --- |
| `200` | a named status | the run completed; the payload is the answer |
| `202` | — | the run was accepted; poll the job |
| `400` | `BAD_REQUEST` | the request is wrong — including an unconfigured indication, refused at creation rather than answered with another indication's results |
| `404` | `NOT_FOUND` | the run completed and does not contain this |
| `409` | `RUN_NOT_COMPLETE` | no completed run yet; the reason names the call that fixes it |
| `409` | `RUN_IN_PROGRESS` | a run is already going for this project |
| `500` | `INTERNAL_ERROR` | a real fault, always with a body — nothing leaves a handler without a response |

`404` and `409` are kept apart deliberately: telling a caller to re-run when the
run already completed and simply lacks their gene sends them round a loop that
cannot help.
