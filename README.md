# CAR-T design platform

Give it a cancer type. It screens the entire human surface proteome for
targets, pairs them, routes each to a receptor architecture, retrieves antibody
binders, assembles CAR constructs, gates them on safety and developability, and
ranks what survives — reporting at every step what it could not measure rather
than filling the gap in.

No target is ever seeded. You supply a diagnosis, not a hypothesis.

---

## What it produces

For pancreatic ductal adenocarcinoma, from 20,431 reviewed human proteins:

| | |
| --- | --- |
| **3,466** | pass the surface filter and are scored |
| **200** | carry forward as the screened pool |
| **19,900** | pairs evaluated for combinatorial designs |
| **191** | blocked on normal-tissue risk |
| **8** | designs fit the payload budget and reach the end |

The eight are `BUILDABLE_AWAITING_BINDER`. They route to an adaptor receptor,
which binds a tag rather than the antigen — and no anti-tag binder exists in the
connected sources, so the platform declares each construct's size and refuses to
invent its sequence.

**That is the result, served as HTTP 200 with reasons.** Never a 404, never an
empty list. An empty list would read as "we looked and had nothing to say";
something specific and measured stopped each design, and naming it is the answer.

A second indication, invasive breast carcinoma, runs the same way from the same
code — different cohort, atlas, dependency lineage and normal denominator, all
resolved from the cancer type.

---

## The rules it holds to

These are not style preferences. Each one exists because the obvious
alternative produced a confident wrong answer at some point in the build.

- **Missing is a third state.** An absent measurement is never zero, never the
  midpoint, never the mean. A protein nobody measured must not be scored as one
  measured and found absent.
- **A bound is not a measurement.** Where a number is the best the evidence
  supports rather than the thing itself, it carries a flag saying so.
- **The two scores stay independent.** Normal-tissue risk and evidence
  confidence are never combined into one number; a well-evidenced dangerous
  target and a poorly-evidenced safe one must not average into the same score.
- **Degradation refuses, it does not caveat.** An indication missing its atlas
  still scores 3,399 of 3,466 genes, but the ranking is unusable — one component
  is the only thing rejecting stromal and immune genes. The platform reports
  `NOT_USABLE` rather than returning numbers with a note attached.
- **Parameters are fixed before output exists**, and criteria are written before
  the run that tests them.

[specs/design-decisions.md](specs/design-decisions.md) records 485 of these
across 47 modules — what was chosen, and what breaks without it. The source
itself carries one-line docstrings and no commentary, so that file is the only
place the reasoning lives.

---

## Running it

Requires Python 3.13 and a provisioned `data/` cache (852 MB; not in git).

```bash
python bootstrap.py --from-release      # fetch and unpack the cache, ~2 min
python bootstrap.py                     # report what is on disk
```

**Use the interpreter in `.venv`.** A bare `python` or `py -3.13` is the system
Python, which lacks this project's dependencies and fails on the first import.

### Verify the whole pipeline

```bash
.venv\Scripts\python.exe run_all.py --fresh
```

Twelve stages end to end, **118 of 124 criteria clear**, about 35 minutes.
The six that trip are recorded limitations, each carrying the decision that
accepted it; a criterion tripping *without* one is reported as a regression
and exits non-zero.

`make_artifact.py` renders the whole run as a single page from the run's own
JSON.

### Serve it

```bash
.venv\Scripts\python.exe -m car_pipeline.api.server
```

`http://127.0.0.1:8000`, loopback only. Nothing needs to be running first — no
database, no broker, no worker. A screen takes minutes, so it is a job-and-poll
API rather than request-response.

---

## Documentation

| | |
| --- | --- |
| [API.md](API.md) | every endpoint, real captured responses, the states a frontend must handle |
| [FRONTEND.md](FRONTEND.md) | what a UI can display today, what is degraded or awaiting a part, and what nothing produces |
| [DEPLOY.md](DEPLOY.md) | fresh clone to a public URL, for someone who has never seen this project |
| [specs/design-decisions.md](specs/design-decisions.md) | why the code is shaped the way it is |
| [specs/](specs/) | one specification per stage, written before its implementation |

---

## Layout

```
car_pipeline/
  configs/     indications: cohort, atlas, lineage, normal denominator
  data/        one connector per source, each cache committed by a manifest
  schemas/     the project contract
  stages/      the pipeline, one module per stage
  api/         job-and-poll HTTP surface over the pipeline
verify_*.py    one script per stage; the criteria that decide whether it passed
run_all.py     runs all twelve as separate processes and folds the results
bootstrap.py   provisions the cache on a fresh clone
```

Stages 7 and 8 are absent from the pipeline, not from the numbering. The
reports show the gap rather than renumbering around it.
