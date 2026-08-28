# API

Every response below is captured verbatim from a live run against the
pancreatic ductal adenocarcinoma cache. Nothing here is illustrative.

- **No authentication.** The service has no notion of a user.
- **No database.** Projects, jobs and results live in the process. A restart
  loses them.
- **JSON in, JSON out.** `Content-Type: application/json` on every POST.
- **A screen takes minutes**, so starting one returns `202` immediately and you
  poll. A synchronous endpoint would time out in every client.

```bash
BASE=http://127.0.0.1:8000        # or your deployed URL
```

---

## Endpoints

| | | |
| --- | --- | --- |
| `GET` | `/indications` | what this deployment can answer for |
| `POST` | `/projects` | create a project from a cancer type |
| `POST` | `/projects/{id}/runs` | start the screen — returns `202` |
| `GET` | `/jobs/{job_id}` | poll progress |
| `GET` | `/projects/{id}/targets?limit=` | the ranked targets |
| `GET` | `/projects/{id}/pairs?limit=` | evaluated pairs |
| `GET` | `/projects/{id}/constructs` | assembled CAR constructs |
| `GET` | `/projects/{id}/result` | the attrition chain and end state |
| `GET` | `/projects/{id}/evidence/{gene}` | every stage's view of one target |
| `GET` | `/projects/{id}/validation` | Mode A verdict on a supplied target |

`{id}` is a 12-character hex project id; `{job_id}` likewise. Anything else on
those paths is a `404`.

---

## The flow

### 0 — Ask what it can answer for

```bash
curl -s "$BASE/indications"
```

```json
{
  "status": "CONFIGURED",
  "indications": [
    {
      "cancer_type": "Invasive Breast Carcinoma",
      "key": "brca",
      "cohort": "TCGA-BRCA",
      "atlas": "GSE176078",
      "dependency_lineage": "Breast",
      "normal_denominator": "Breast_Mammary_Tissue"
    },
    {
      "cancer_type": "Pancreatic Ductal Adenocarcinoma",
      "key": "pdac",
      "cohort": "TCGA-PAAD",
      "atlas": "GSE202051",
      "dependency_lineage": "Pancreas",
      "normal_denominator": "Pancreas"
    }
  ],
  "reasons": [
    "An indication needs a tumour cohort, a single-cell atlas, a dependency lineage and a normal-tissue denominator declared before it can be screened. None of those is derivable from the cancer type, so an unregistered one is refused rather than answered with another indication's results."
  ]
}
```

Call this first rather than guessing. Populate a dropdown from it — anything
not on the list is refused at creation.

### 1 — Create a project

```bash
curl -s -X POST "$BASE/projects" \
  -H 'Content-Type: application/json' \
  -d '{"cancer_type":"Pancreatic Ductal Adenocarcinoma"}'
```

```json
{
  "project_id": "944d0f1d854a",
  "cancer_type": "Pancreatic Ductal Adenocarcinoma",
  "target_antigen": null,
  "discovery_mode": "B",
  "created_at": "2026-08-28T19:50:59+00:00"
}
```

**`target_antigen` is null on purpose and stays null.** That null is what
selects discovery mode. Sending a target instead selects Mode A — see
[Mode A](#mode-a-validating-a-target-you-already-have) below.

Short aliases resolve: `PDAC`, `breast`. Matching is case-insensitive.

### 2 — Start the screen

```bash
PID=944d0f1d854a
curl -s -X POST "$BASE/projects/$PID/runs"
```

```json
{
  "job_id": "452515b67e75",
  "project_id": "944d0f1d854a",
  "status": "running",
  "stage": "sources",
  "note": "loading cached sources",
  "started_at": "2026-08-28T19:50:59+00:00",
  "finished_at": null,
  "error": null,
  "stages": ["sources", "screen", "pairing", "binders", "constructs",
             "safety", "developability", "ranking"]
}
```

HTTP `202`. The `stages` array is the full sequence, in order — use it to build
a progress indicator rather than hardcoding stage names.

### 3 — Poll

```bash
JID=452515b67e75
until curl -s "$BASE/jobs/$JID" | grep -qE '"status": "(complete|failed)"'; do
  sleep 5
done
```

> **Match both terminal states.** Waiting only for `complete` spins forever on a
> failed job while discarding the body that explains it.

Mid-run:

```json
{
  "job_id": "452515b67e75",
  "status": "running",
  "stage": "screen",
  "note": "ranking the surface proteome",
  "finished_at": null,
  "error": null
}
```

Finished — this run took **355 seconds**:

```json
{
  "job_id": "452515b67e75",
  "status": "complete",
  "stage": "ranking",
  "note": "RANKED_AWAITING_BINDER",
  "started_at": "2026-08-28T19:50:59+00:00",
  "finished_at": "2026-08-28T19:56:54+00:00",
  "error": null
}
```

On completion `note` carries the run's end state, so a client that only polls
already knows the outcome before fetching anything.

### 4 — Ranked targets

```bash
curl -s "$BASE/projects/$PID/targets?limit=3"
```

```json
{
  "status": "RANKED",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "universe": 3466,
  "scored": 3400,
  "ceiling": 0.15,
  "returned": 3,
  "targets": [
    {
      "rank": 1,
      "gene": "CEACAM5",
      "accession": "P06731",
      "evidence_class": "PROTEIN_CONFIRMED",
      "tier_rank": 1,
      "composite": 0.8769,
      "measured_weight": 0.55,
      "risk": 0.6,
      "risk_organ": "gi_tract",
      "cleared": false,
      "confidence": 0.685,
      "breakdown": {
        "malignant_expression": null,
        "malignant_vs_stroma": null,
        "tumour_vs_normal": 1.0,
        "patient_prevalence": 0.8820224719101124,
        "surface_accessibility": 1.0,
        "escape_resistance": 0.0
      }
    }
  ],
  "reasons": []
}
```

> ### `null` in `breakdown` means unmeasured. It never means zero.
>
> The top-ranked target has two null components — its transcripts were not
> captured in the atlas — and `measured_weight: 0.55` says the composite rests
> on 55% of the evidence. **Render nulls as "not measured", never as 0, and
> never omit the row.** A frontend that drops null components silently converts
> "we do not know" into "we looked and found nothing".
>
> `risk` and `confidence` are separate numbers and must stay separate. A
> well-evidenced dangerous target and a poorly-evidenced safe one must not
> collapse into the same badge.

`cleared: false` means the target exceeds the normal-tissue ceiling on its own.
That is not a rejection of the design — pairing and routing come next.

### 5 — Pairs

```bash
curl -s "$BASE/projects/$PID/pairs?limit=3"
```

```json
{
  "status": "PAIRED",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "evaluated": 19900,
  "measured": 14535,
  "returned": 3,
  "pairs": [
    {
      "gene_a": "MUC16",
      "gene_b": "MUC17",
      "combined_risk": 0.0048,
      "peak_organ": "endocrine",
      "cleared": true,
      "coverage_f_ab": 0.006321856890514115,
      "coverage_span_kb": 79.159,
      "coverage_span_percentile": 0.0785,
      "coverage_caveat": "span-confounded; read the percentile beside the fraction"
    }
  ],
  "reasons": [
    "Coverage does not gate: per-cell detection tracks genomic span (rho +0.68) more strongly than expression (+0.20).",
    "Stage 4 closed at 10 of 14 criteria with four documented limitations."
  ]
}
```

**Never show `coverage_f_ab` without `coverage_span_percentile`.** The raw
fraction correlates with how long the gene is more than with how much is
expressed. `coverage_caveat` ships on every pair for exactly this reason — it is
not decoration.

### 6 — Constructs

```bash
curl -s "$BASE/projects/$PID/constructs"
```

```json
{
  "status": "BUILDABLE_AWAITING_BINDER",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "counts": {"NO_CONSTRUCT": 190, "BUILDABLE": 8, "BUDGET_EXCEEDED": 2},
  "buildable": 8,
  "complete": 0,
  "awaiting_binder": 8,
  "over_budget": 2,
  "constructs": [
    {
      "gene": "CD207",
      "partner": null,
      "verdict": "BUILDABLE",
      "state": "AWAITING_BINDER",
      "architecture": "adaptor, anti-tag receptor, antigen on the adaptor",
      "binder_supplied": false,
      "binder": "anti-tag binder (sequence not supplied)",
      "total_bp": 2811,
      "budget_bp": 3500,
      "headroom_bp": 689,
      "amino_acid_sequence": null,
      "dna": null,
      "domains": [
        {
          "name": "CD8A leader",
          "provenance": "proteome",
          "accession": "P01732",
          "feature": "Signal",
          "source_residues": "1-21",
          "aa_start": 0, "aa_end": 21,
          "bp_start": 0, "bp_end": 63
        }
      ]
    }
  ]
}
```

Each construct is in one of three states, and they are **not** interchangeable:

| `state` | meaning |
| --- | --- |
| `COMPLETE` | has residues; can be ordered |
| `AWAITING_BINDER` | has a layout, a length and a domain map, and no residues |
| `BUDGET_EXCEEDED` | assembled, and does not fit the payload budget |

`AWAITING_BINDER` is the interesting one. The design routes to an adaptor
receptor, which binds a tag rather than the antigen, and no anti-tag binder
exists in the connected sources. The platform declares the construct's size and
**refuses to invent its sequence** — `amino_acid_sequence` and `dna` are null,
and that null is load-bearing. A frontend that renders it as an empty string,
or offers a download, is claiming something the platform explicitly declined to
claim.

### 7 — The end state

```bash
curl -s "$BASE/projects/$PID/result"
```

```json
{
  "status": "RANKED_AWAITING_BINDER",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "pool_size": 200,
  "reached_the_end": 8,
  "complete": 0,
  "awaiting_binder": 8,
  "attrition": [
    {"gate": "blocked on normal tissue risk", "dropped": 191, "remaining": 9},
    {"gate": "no design recommended",         "dropped": 0,   "remaining": 9},
    {"gate": "no binder retrieved",           "dropped": 1,   "remaining": 8},
    {"gate": "no construct assembled",        "dropped": 0,   "remaining": 8},
    {"gate": "construct over budget",         "dropped": 0,   "remaining": 8}
  ],
  "reasons": [
    "8 design(s) reached the end: 0 complete, 8 awaiting a binder sequence.",
    "Every drop is a measurement, not a failure of the stage that made it.",
    "A design awaiting a binder has a layout, a length and a domain map, and no residues for its anti-tag binder. It fits the budget and cannot be ordered yet, which is a different state both from a design that does not fit and from a finished one."
  ],
  "developability_status": "SCORED"
}
```

Every candidate is attributed to the **first** gate it failed, so `dropped` sums
to `pool_size` and never double-counts. `remaining` is the running total — it is
the series to chart.

### 8 — Evidence for one target

```bash
curl -s "$BASE/projects/$PID/evidence/MSLN"
```

Returns every stage's view of a single gene under `stage3_screen`,
`stage4_pairing`, `stage5_binders`, `stage6_construct`, `stage9_safety`,
`stage10_developability` and `stage11_ranking`, each with the provenance behind
it — including `route_reason`, the sentence explaining why the target routed
where it did:

```json
{
  "status": "EVIDENCE",
  "gene": "MSLN",
  "accession": "Q13421",
  "stage4_pairing": {
    "outcome": "NO_DESIGN",
    "architecture": "NO_ARCHITECTURE",
    "route_reason": "risk 0.6366 exceeds every declared ceiling; no measured pair",
    "failed_on": {"risk": 199, "coverage_below_floor": 0,
                  "unmeasured": 29, "partner_ineligible": 0}
  }
}
```

---

## Mode A: validating a target you already have

Supply `target_antigen` at creation. The pipeline runs the same screen, then
answers one question: is this target suitable for this indication?

```bash
curl -s -X POST "$BASE/projects" -H 'Content-Type: application/json' \
  -d '{"cancer_type":"Pancreatic Ductal Adenocarcinoma","target_antigen":"ERBB2"}'
# then POST /runs, poll, and:
curl -s "$BASE/projects/$PID_A/validation"
```

```json
{
  "status": "UNSUITABLE",
  "mode": "A",
  "cancer_type": "Pancreatic Ductal Adenocarcinoma",
  "target": "ERBB2",
  "accession": "P04626",
  "rank": 124,
  "of": 3400,
  "composite": 0.4984,
  "measured_weight": 1.0,
  "risk": 0.6565,
  "risk_organ": "kidney",
  "evidence_class": "PROTEIN_CONFIRMED",
  "architecture": "NO_ARCHITECTURE",
  "reasons": [
    "Ranks 124 of 3400 on tumour attractiveness (composite 0.4984, 1.00 of the evidence measured).",
    "Normal-tissue risk 0.6565, peak organ kidney.",
    "risk 0.6565 exceeds every declared ceiling; no measured pair",
    "57 binder candidate(s) retrieved; 20 carry a usable sequence."
  ]
}
```

| `status` | meaning |
| --- | --- |
| `SUITABLE` | routes to a conventional receptor |
| `CONDITIONAL` | routes only to a gated or adaptor architecture |
| `UNSUITABLE` | in the proteome, but no architecture admits it |
| `NOT_ASSESSED` | did not reach the screened pool, so no design was evaluated |
| `UNKNOWN_TARGET` | not a gene symbol in the reviewed proteome |
| `NOT_APPLICABLE` | this project supplied no target and ran in discovery mode |

Calling `/validation` on a discovery-mode project is not an error:

```json
{
  "status": "NOT_APPLICABLE",
  "reasons": ["This project supplied no target, so it ran in discovery mode. Create a project with target_antigen to ask whether a specific target is suitable."]
}
```

---

## States a frontend must handle

### On every collection response

`/targets`, `/pairs`, `/constructs` and `/result` all carry the same three
run-level fields, so a caller never has to know to ask somewhere else whether
the numbers beneath it are supported:

| field | meaning |
| --- | --- |
| `usability` | `USABLE` or `NOT_USABLE` |
| `unavailable` | components that could not be measured, each naming the source it needed |
| `indication` | the resolved cancer type |
| `degraded_note` | present **only** when `unavailable` is non-empty |

**`NOT_USABLE` is a refusal, not an empty result.** Where a missing source
undermines the ranking, `/targets` returns `"targets": []` with
`"status": "NOT_USABLE"` and reasons. The empty array is not "nothing scored
well" — most genes still score. It means the ranking that could be produced
would not be trustworthy, so it is withheld. Render the reasons; do not render
an empty table.

This happens when an indication has no single-cell atlas: one scoring component
is the only thing rejecting stromal and immune genes, so without it a ranking
fills with immunoglobulin and MHC-II genes and still looks like an answer.

### Terminal statuses

| endpoint | possible `status` |
| --- | --- |
| `/jobs/{id}` | `queued`, `running`, `complete`, `failed` |
| `/targets` | `RANKED`, `NOT_USABLE` |
| `/pairs` | `PAIRED` |
| `/constructs` | `BUILDABLE`, `BUILDABLE_AWAITING_BINDER`, `NO_BUILDABLE_CONSTRUCT` |
| `/result` | `RANKED`, `RANKED_AWAITING_BINDER`, `NO_DESIGN_REACHES_THE_END` |
| `/validation` | the six in the table above |
| `/evidence/{gene}` | `EVIDENCE` |
| `/indications` | `CONFIGURED` |

`NO_BUILDABLE_CONSTRUCT` and `NO_DESIGN_REACHES_THE_END` are **HTTP 200
results, not errors.** They arrive with counts and reasons attached. Never a
404, never a 500, never a bare `[]` — an empty list would read as "we looked
and had nothing to say", when something specific and measured stopped each
design.

### Errors

| code | `status` | when |
| --- | --- | --- |
| `400` | `BAD_REQUEST` | blank or missing `cancer_type`; an unconfigured indication |
| `404` | `NOT_FOUND` | no project with this id |
| `404` | `NOT_FOUND` | no job with this id |
| `404` | `NOT_FOUND` | a gene that is not in the completed run |
| `409` | `RUN_NOT_COMPLETE` | the project exists; its run has not finished |
| `409` | `RUN_IN_PROGRESS` | a second run started while one is active |

An unconfigured indication names what is available rather than failing blankly:

```json
{
  "status": "BAD_REQUEST",
  "error": "no configuration for 'Glioblastoma'. Configured indications are: Invasive Breast Carcinoma; Pancreatic Ductal Adenocarcinoma. An indication needs a tumour cohort, a single-cell atlas, a dependency lineage and a normal-tissue denominator declared before it can be screened; none of those is derivable from the name."
}
```

### `404` and `409` are different questions

A bad id and a run in progress are separate conditions and answer separately,
on every result endpoint — `/targets`, `/pairs`, `/constructs`, `/result`,
`/validation` and `/evidence/{gene}`.

**`404` — no such project. Stop polling.**

```json
{
  "status": "NOT_FOUND",
  "error": "ffffffffffff",
  "reasons": ["No project with this id. It was never created, or the service restarted: projects live in memory and do not survive one. POST /projects to create one. This is not a run that has yet to finish."]
}
```

**`409` — the project exists and its run has not finished. Keep polling.**

```json
{
  "status": "RUN_NOT_COMPLETE",
  "reasons": ["This project exists and has no completed run. POST /projects/{id}/runs, then poll /jobs/{job_id}."]
}
```

A UI must branch on these: `404` is a dead end and should surface "no such
project", while `409` means come back. Treating them alike either polls a typo
forever or tells someone their finished run does not exist.

`/jobs/{job_id}` answers `404` the same way, with the same distinction — a job
id that never existed, or one lost to a restart:

```json
{
  "status": "NOT_FOUND",
  "error": "ffffffffffff",
  "reasons": ["No job with this id. Jobs live in memory and do not survive a restart."]
}
```

> **Both are also what a restart looks like.** Projects and jobs are held in
> the process, so an id that worked a minute ago can become a `404`. That is
> the same response as a typo and is not distinguishable from one; a client
> holding an id across a possible restart should be prepared to create the
> project again.

### One run at a time

The run guard is global, not per project: a second `POST /runs` while any run is
active returns `409 RUN_IN_PROGRESS`. A deployment must therefore serve with
concurrency 1 — see [DEPLOY.md](DEPLOY.md). Jobs live in memory and do not
survive a restart.
