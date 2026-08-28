# What a frontend can display today

For whoever builds the UI. Three groups: what is there, what is there but
incomplete, and what nothing produces.

**Provenance.** Groups 1 and 2 are taken from [API.md](API.md) — every field
below appears in a response captured from a live run. Group 3 cannot come from
API.md, because a document describing what the API returns cannot describe what
it does not; it is taken from the specifications in [specs/](specs/) and
verified against the code that would have to produce it. Where a group 3 entry
names a file, that file was read to confirm the absence rather than inferred.

The reference document itself (`AI_Pipeline.pdf`) is **not in this repository**.
Group 3 is assembled from what the specs record of it — chiefly
[specs/gaps-architecture-and-scoring.md](specs/gaps-architecture-and-scoring.md),
which is titled "Two gaps against AI_Pipeline.pdf". If you need group 3 checked
against the document itself, that check has not been done here.

---

# 1. Data that exists

Complete, measured, and safe to render as a value.

## Indication picker — `GET /indications`

| field | type |
| --- | --- |
| `status` | `CONFIGURED` |
| `indications[].cancer_type` | string — **this is the value to POST** |
| `indications[].key` | string (`pdac`, `brca`) |
| `indications[].cohort` | string (`TCGA-PAAD`) |
| `indications[].atlas` | string (`GSE202051`) |
| `indications[].dependency_lineage` | string |
| `indications[].normal_denominator` | string |
| `reasons[]` | array of strings |

Two indications today. Populate a picker from this and do not offer free text —
anything else is a `400` at creation. Short aliases (`PDAC`, `breast`) resolve,
case-insensitively.

## Project — `POST /projects` → `201`

| field | type |
| --- | --- |
| `project_id` | string, 12 hex |
| `cancer_type` | string, the **resolved** name — echo this, not what was typed |
| `target_antigen` | null in discovery mode; the symbol in validation mode |
| `discovery_mode` | `A` or `B` |
| `created_at` | ISO-8601 with offset |

## Run progress — `POST /projects/{id}/runs` → `202`, then `GET /jobs/{job_id}`

| field | type |
| --- | --- |
| `job_id`, `project_id` | string |
| `status` | `queued` \| `running` \| `complete` \| `failed` |
| `stage` | current stage name |
| `note` | prose mid-run; **the run's end state on completion** |
| `stages[]` | the full ordered sequence — build the progress bar from this, do not hardcode |
| `started_at`, `finished_at` | ISO-8601; `finished_at` null until done |
| `error` | null unless the run failed |

Eight stages: `sources`, `screen`, `pairing`, `binders`, `constructs`,
`safety`, `developability`, `ranking`. A captured run took **355 seconds**.

## Ranked targets — `GET /projects/{id}/targets?limit=N`

| field | type |
| --- | --- |
| `status` | `RANKED` |
| `universe`, `scored`, `returned` | 3466, 3400, N |
| `ceiling` | the normal-tissue risk ceiling in force (0.15) |
| `targets[].rank`, `.gene`, `.accession` | integer, symbol, accession |
| `targets[].evidence_class` | e.g. `PROTEIN_CONFIRMED` |
| `targets[].tier_rank` | integer |
| `targets[].composite` | 0–1 — **read with `measured_weight`, see group 2** |
| `targets[].risk`, `.risk_organ` | 0–1 and the peak organ |
| `targets[].confidence` | 0–1 |
| `targets[].cleared` | boolean |
| `targets[].breakdown` | six named components — **nullable, see group 2** |

`limit` is a query parameter; `returned` reflects it.

## Pairs — `GET /projects/{id}/pairs?limit=N`

| field | type |
| --- | --- |
| `status` | `PAIRED` |
| `evaluated`, `measured`, `returned` | 19900, 14535, N |
| `pairs[].gene_a`, `.gene_b` | symbols |
| `pairs[].combined_risk`, `.peak_organ` | 0–1 and the organ |
| `pairs[].cleared` | boolean |
| `pairs[].coverage_span_kb` | genomic span |
| `reasons[]` | why coverage does not gate |

## Constructs — `GET /projects/{id}/constructs`

| field | type |
| --- | --- |
| `status` | `BUILDABLE` \| `BUILDABLE_AWAITING_BINDER` \| `NO_BUILDABLE_CONSTRUCT` |
| `counts` | `{NO_CONSTRUCT, BUILDABLE, BUDGET_EXCEEDED}` |
| `buildable`, `complete`, `awaiting_binder`, `over_budget` | integers |
| `constructs[].gene`, `.partner` | symbol; partner null on a single design |
| `constructs[].verdict`, `.state` | see the state table in group 2 |
| `constructs[].architecture` | **prose**, e.g. "adaptor, anti-tag receptor, antigen on the adaptor" |
| `constructs[].total_bp`, `.budget_bp`, `.headroom_bp` | 2811, 3500, 689 |
| `constructs[].domains[]` | the full domain map — **this is the richest thing to draw** |
| `constructs[].domains[].name`, `.provenance`, `.accession`, `.feature` | strings |
| `constructs[].domains[].source_residues` | e.g. `"1-21"` |
| `constructs[].domains[].aa_start`, `.aa_end`, `.bp_start`, `.bp_end` | integers |

The domain map is complete and contiguous even when the sequence is not
supplied — it is enough to render a scaled construct diagram with labelled,
sourced segments.

## End state — `GET /projects/{id}/result`

| field | type |
| --- | --- |
| `status` | `RANKED` \| `RANKED_AWAITING_BINDER` \| `NO_DESIGN_REACHES_THE_END` |
| `pool_size`, `reached_the_end`, `complete`, `awaiting_binder` | 200, 8, 0, 8 |
| `attrition[].gate`, `.dropped`, `.remaining` | five gates, in order |
| `developability_status` | `SCORED` |
| `reasons[]` | prose |

Every candidate is attributed to the **first** gate it failed, so `dropped` sums
to `pool_size`. `remaining` is the series to chart.

## Evidence — `GET /projects/{id}/evidence/{gene}`

`status: EVIDENCE`, plus `gene`, `accession`, and seven per-stage objects:
`stage3_screen`, `stage4_pairing`, `stage5_binders`, `stage6_construct`,
`stage9_safety`, `stage10_developability`, `stage11_ranking`.

`stage4_pairing` is the one API.md expands: `outcome`, `architecture`,
`route_reason` (a sentence explaining the routing decision), and `failed_on`
with `risk`, `coverage_below_floor`, `unmeasured`, `partner_ineligible`.

## Validation, Mode A — `GET /projects/{id}/validation`

`status`, `mode`, `cancer_type`, `target`, `accession`, `rank`, `of`,
`composite`, `measured_weight`, `risk`, `risk_organ`, `evidence_class`,
`architecture`, `reasons[]`.

> **On `NOT_APPLICABLE` only `status` and `reasons` are present.** Every other
> field is absent, not null. A UI that reads `target` unconditionally will
> throw on a discovery-mode project.

## Run-level, on every collection response

`usability`, `unavailable`, `indication` ride on `/targets`, `/pairs`,
`/constructs` and `/result`. **They are not on `/validation` or
`/evidence/{gene}`** — do not build a shared header component that expects them
everywhere.

---

# 2. Data that exists but is degraded, bounded, or awaiting a part

Present in the response, and **not** a finished measurement. Rendering any of
these as a plain value states something the platform declined to state.

## Unmeasured — null that must never become zero

| where | what |
| --- | --- |
| `targets[].breakdown.*` | any of the six components can be null. In the captured run the **top-ranked target has two nulls** |
| `targets[].measured_weight` | the share of evidence actually measured — `0.55` on that top target |
| `stage4_pairing.failed_on.unmeasured` | 29 pairs failed for want of a measurement, not on a value |
| `/pairs` `evaluated` vs `measured` | 19,900 evaluated, 14,535 measured |

**Render null as "not measured". Never as 0, never as a blank cell, never
omitted.** Dropping the row converts "we do not know" into "we looked and found
nothing", which is the one inversion this whole platform is built to prevent.

Show `measured_weight` beside `composite` wherever the composite appears. A
composite at 0.55 measured weight and one at 1.0 are not comparable numbers.

Keep `risk` and `confidence` visually separate. They are deliberately never
combined; a single badge merging them re-introduces the error.

## Bounded — a number that needs its qualifier beside it

| where | qualifier |
| --- | --- |
| `pairs[].coverage_f_ab` | `coverage_span_percentile` **and** `coverage_span_kb` |
| `pairs[].coverage_caveat` | ships on every pair: "span-confounded; read the percentile beside the fraction" |
| `targets[].composite` | `measured_weight` |
| `/targets` `scored` vs `universe` | 3,400 of 3,466 were scorable |

**Never show `coverage_f_ab` alone.** It tracks how long the gene is more than
how much is expressed. The caveat string is not decoration — it is shipped per
pair so it cannot be separated from the number.

## Awaiting a part — real objects missing something not yet supplied

| state | meaning |
| --- | --- |
| `constructs[].state = COMPLETE` | has residues; can be ordered |
| `constructs[].state = AWAITING_BINDER` | layout, length and domain map; **no residues** |
| `constructs[].state = BUDGET_EXCEEDED` | assembled and does not fit the budget |

All three are distinct and must render distinctly. In the captured run **all
eight surviving designs are `AWAITING_BINDER`**: `amino_acid_sequence` and
`dna` are null, `binder_supplied` is false, and `binder` reads "anti-tag binder
(sequence not supplied)".

> Do not render a null sequence as an empty string, a spinner, or a download
> button. The platform declared the construct's size and **refused to invent
> its sequence**; a download control claims the opposite.

`/result` mirrors this at run level: `status: RANKED_AWAITING_BINDER`,
`reached_the_end: 8`, `complete: 0`.

`/validation` `status: CONDITIONAL` is the same shape one level up — suitable,
but only under a restricted architecture set.

## Degraded — the run itself may not support a ranking

| field | meaning |
| --- | --- |
| `usability` | `USABLE` or `NOT_USABLE` |
| `unavailable[]` | components that could not be measured, each naming the source it needed |
| `degraded_note` | **present only when `unavailable` is non-empty** |

> **`NOT_USABLE` is a refusal, not an empty result.** `/targets` returns
> `"targets": []` with reasons. Most genes still scored — the empty array means
> the ranking that could be produced would not be trustworthy, so it is
> withheld. **Render the reasons; do not render an empty table.**

Both are `[]` / absent for the two configured indications today. You will only
see them on an indication lacking a source — but build for them, because
degradation is a headline capability and the next indication may trigger it.

## Not-a-result statuses that still arrive as HTTP 200

`NO_BUILDABLE_CONSTRUCT`, `NO_DESIGN_REACHES_THE_END`, and validation's
`UNSUITABLE`, `NOT_ASSESSED`, `UNKNOWN_TARGET`, `NOT_APPLICABLE`. All carry
counts and reasons. **None is an error.** Never a 404, never a bare `[]`.

## Error states

| code | status | what a UI does |
| --- | --- | --- |
| `400` | `BAD_REQUEST` | show `error` — it names every configured indication |
| `404` | `NOT_FOUND` | no such project / job / gene. **Stop polling.** |
| `409` | `RUN_NOT_COMPLETE` | project exists, run unfinished. **Keep polling.** |
| `409` | `RUN_IN_PROGRESS` | the run guard is **global**, not per project |

> **Both 404s are also what a restart looks like.** Projects and jobs are held
> in memory. An id that worked a minute ago can 404, indistinguishably from a
> typo, and it is unrecoverable — there is no database. A client holding an id
> across a possible restart must be prepared to create the project again.

---

# 3. What no stage produces

Nothing here can be displayed, because nothing generates it. Do not design
screens around these.

## Whole stages that were specified and never built

| | |
| --- | --- |
| **Stage 7 — Manufacturing** | schema only. No endpoint, no field, no job stage. |
| **Stage 8 — Trial design** | schema only. Same. |

Declared absent in `make_artifact.py` (`ABSENT`) and shown as a "not
implemented" row on the report page rather than renumbered away. The job
`stages[]` array goes `constructs → safety`, and `/evidence/{gene}` jumps
`stage6_construct → stage9_safety`.

## Fields that exist but are permanently a sentinel

These are served, so a UI *will* receive them — as a fixed string, never a
value. Treat them as absent, not as data with a placeholder.

| field | always |
| --- | --- |
| `stage5_binders.sequence_candidates[].affinity` | `"NOT_CONNECTED"` — no affinity source is wired |
| `stage5_binders.sequence_candidates[].isoform` | `"ISOFORM_UNRESOLVED"` |
| `stage9_safety.epitope_immunogenicity` | `"NOT_CONNECTED"` |

## Scoring terms the reference asks for that are not computed

The target `breakdown` has exactly **six** components. Two more are specified
and absent:

- **Shedding / soluble antigen (term R).** No field, no flag. A shed antigen
  meets its binder in plasma rather than on a cell; nothing scores this.
- **Antigen stability under treatment (term A).** No treated-versus-untreated
  ratio is computed anywhere.

Both are recorded as GAP 2 in
[specs/gaps-architecture-and-scoring.md](specs/gaps-architecture-and-scoring.md).

## Architectures that cannot appear

`architecture` is reported two different ways — an **enum** on
`stage4_pairing` and `/validation`, and **prose** on `constructs[]`. The enum
values a client can actually receive are only:

`CONVENTIONAL`, `AND_GATE`, `ADAPTOR`, `NO_ARCHITECTURE`, `NOT_CONFIGURED`.

Specified and unreachable:

| architecture | why it never arrives |
| --- | --- |
| `AND_NOT` / inhibitory (iCAR) | `route()` returns it only when passed an exclusion marker, and **no caller ever passes one** — there is no exclusion-antigen source |
| `SWITCHABLE` (rapamycin ON-switch) | defined with an unbuilt reason, and **no branch of `route()` returns it** |
| OR-gated, tandem, bicistronic | no representation at all — not offered, and not refused by name |

The first two are the ones to know about: they are *named* in the code with
reasons, which makes them look supported. They are not. A heterogeneity-driven
or exclusion-based design has no way to be expressed or reported.

## Inputs the API will not accept

`POST /projects` takes **`cancer_type` and `target_antigen` only.** The project
schema declares much more — `malignancy_type`, `cancer_subtype`,
`patient_subgroup`, `product_type`, `car_format`, `safety_tolerance`,
manufacturing limits, tissue overrides, and an `existing_binder` with a format
enum — and **none of it is reachable over HTTP**.

> Sending those fields returns `201` and **silently ignores them.** There is no
> error. Do not build a configuration form: every control on it would appear to
> work and change nothing.

In particular, a caller **cannot supply their own binder**, which is the most
likely thing a user would want given eight designs are awaiting one.

## Other absences

- **No project-spec or dataset-availability endpoint.** `ProjectSpec` declares
  `required_datasets` with per-dataset status (`available` / `unreachable` /
  `not_configured`); no endpoint serves it. `/indications` exposes what each
  indication is *backed by* (cohort, atlas, lineage, denominator) but not
  whether those datasets are readable right now.
- **No arbitrary cancer type.** The reference asks for a platform accepting any
  cancer type; two are configured. A third needs its cohort, atlas, lineage and
  denominator declared and its caches provisioned — it is not a UI-side gap,
  and the API refuses cleanly rather than guessing.
- **No hematological, allogeneic or armored design space**, and no genetic-edit
  count is reported. Both configured indications are solid and autologous.

---

## The five rules, condensed

1. Null in a `breakdown` means **not measured**. Never 0, never blank, never dropped.
2. `coverage_f_ab` never appears without `coverage_span_percentile`.
3. `risk` and `confidence` never merge into one badge.
4. A null sequence on an `AWAITING_BINDER` construct is a **refusal**, not a pending value.
5. `404` stop, `409` keep polling.
