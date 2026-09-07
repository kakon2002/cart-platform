# Orientation

For someone joining who has not seen this before. It ties the pieces together;
it does not restate the specifications, which are listed in §8.

Read §1 and §9 first. §1 is what the platform will and will not tell you. §9 is
the thing this project knows that is hardest to learn any other way.

---

## 1. What this is, and what it is not

A cancer-agnostic CAR-T design platform. Given an indication and no target, it
screens the human surface proteome, selects targets, routes an architecture,
retrieves binders, assembles constructs, gates them on normal-tissue safety,
and produces a ranked candidate package with the evidence behind every number.

**What it does not do, and will not pretend to:**

- It does not measure anything on cells. Normal-tissue risk is derived from
  expression data. Every validation figure reads `TO_BE_MEASURED` and every
  acceptance threshold reads `TO_BE_SET_BEFORE_THE_RUN`.
- It does not predict structure or function. Those stages do not exist (§10).
- It does not generate binders. Binder discovery retrieves and stops.
- It does not fill a gap with a plausible value. Missing is a third state
  everywhere, and a component with no measurement reads `UNKNOWN` with the
  reason, never `0.00`.

**Two indications are configured.** Pancreatic ductal adenocarcinoma returns
five candidate designs. Invasive breast carcinoma returns none and says why.
Both are correct outputs, and the second is the one that exercises the
platform's refusal paths.

---

## 2. The eleven stages

The client's reference document describes fourteen. Eleven exist. Stages 7, 8,
13 and 14 are covered in §10.

| | stage | takes | produces |
| --- | --- | --- | --- |
| 1 | Design spec | a project definition | the run's constraints: payload budget, risk ceilings, allowed architectures, required datasets |
| 2 | Surface proteome | the reviewed human proteome | ~3,466 surface records with membrane class, topology and mature chains |
| 3 | Target discovery | surface records, expression, atlas, cohort, dependency | a ranked pool of 200 with six weighted components, a risk figure and its peak organ |
| 4 | Target pairing | the pool | ~19,110 measured pairs with combined risk and coverage |
| 4a | Architecture routing | pairs and tolerances | one routed outcome per target: SINGLE, DUAL, ADAPTOR, NO_DESIGN |
| 5 | Binder discovery | routed decisions | retrieved binders by two routes, with provenance |
| 6 | Construct assembly | decisions and binders | a receptor with a domain map, an amino-acid sequence and a DNA map |
| 9 | Safety gate | constructs, risks, trials | a verdict against the ceiling the design is judged by |
| 10 | Developability | binder sequences | sequence liability flags, never summed |
| 11 | Final ranking | everything above | a Pareto front, a Level B score, and a decision per candidate |
| 12 | Candidate package | the run | one document per candidate, with a named-gaps section |

**Two things about the chain.** Each stage's configuration hash covers the stage
before it, so a change anywhere upstream invalidates everything downstream and
no cached artifact can be reused across a configuration change. And Stage 4a is
not a numbered stage in the reference document — it exists because routing an
architecture is a decision, and a decision needs its own criteria.

---

## 3. The nine connected sources

The brief said six; there are nine. Each is pinned to a release, and the pin is
part of the configuration hash.

| source | pin | why pinned there |
| --- | --- | --- |
| UniProt | `2026_02` | The proteome defines what "surface" means. A release change moves the eligible universe, so it is pinned and a change is a deliberate re-baseline |
| Human Protein Atlas | `v23` | Normal-tissue staining is the veto arm. Its levels are ordinal and re-calibrated between releases; mixing two would compare incomparable grades |
| GTEx | `v10` | The transcript baseline the staining arm is calibrated against. Pinned to the release the calibration was fitted on |
| GENCODE | `47` | Gene spans, used to detect span confounding in pair coverage. Pinned to match the cohort's annotation |
| SAbDab | `2.0.10` | Antibody structures and named therapeutics. New depositions change what is retrievable, so recall is only comparable within a release |
| ClinicalTrials.gov | `v2` | The registry API version. Query-keyed rather than immutable: rebuilt when the screened set changes |
| GDC TCGA | `TCGA-PAAD` | Per indication. The tumour cohort; there is no default, and falling back would screen one indication against another's data |
| DepMap | `24Q4 / Pancreas` | Per indication. Dependency lineage |
| Single-cell atlas | `GSE202051` | Per indication. Supplies malignant-versus-stroma, without which nothing distinguishes a tumour antigen from a gene expressed by the lymphocytes beside it |

**The cache is ~12 GB**, of which 11 GB is the single-cell atlas, used at build
time and not at serve time. `bootstrap.py` provisions it; `run_all.py` verifies
every payload against its manifest before any stage runs, and compares a
size-and-mtime fingerprint afterwards so "read-only" is evidence rather than an
assurance.

---

## 4. The standing rules, and what breaks without each

These are not style preferences. Each exists because its absence produced a
specific failure.

**Specification before implementation, reviewed.** Without it, the design is
whatever the code did, and the specification becomes a description of an
accident rather than a decision anyone took.

**Criteria written before the run, with positive pins.** A criterion written
after seeing the output encodes the output. Positive pins — a case that *must*
appear — are what stops a criterion passing on an empty set.

**Parameters fixed before the output exists.** Never tuned until a known target
ranks well. Without this, the platform reproduces the answer it was tuned
toward and calls it a result.

**A tripped criterion changes the specification and the run repeats.** It is
never explained away. Without this, the criteria become commentary.

**Missing is a third state, never imputed.** `MEASURED`, `UNKNOWN`,
`NOT_APPLICABLE`. Imputing a plausible value makes an unmeasured candidate
indistinguishable from a measured one, and the difference between them is the
whole product. This rule has been broken twice and caught twice, most recently
inside the Pareto front (§9).

**Normal-tissue risk and evidence confidence are never combined.** A
well-evidenced dangerous target and a poorly-evidenced safe one must not reach
the same number. Confidence is a multiplier outside the sum, never a summand.

**Report every count beside the count expected.** A number alone cannot be
wrong; a number beside its prediction can.

**Commit after each step, with the reasoning in the message.** The commit graph
is the audit trail — the binder benchmark's blind rule depends on it literally
(§6).

---

## 5. Running it

```
.venv\Scripts\python.exe run_all.py --fresh    # every stage, both indications
.venv\Scripts\python.exe serve.py --port 8080  # the HTTP surface
```

Setup, data provisioning and dashboard wiring are in `HANDOVER.md`, the
operational companion to this document. Use the
virtual environment's interpreter for everything; a bare `python` is the system
interpreter and has none of the three dependencies.

**Testing.** Sixteen stages, **216 criteria**, ~28 minutes from cleared derived
artifacts. Each verifier prints its criteria and stops on the first trip.

**Adding an indication.** Add a module under `car_pipeline/configs/` declaring
a tumour cohort, a single-cell atlas, a dependency lineage and a normal-tissue
denominator, register it in `registry.py`, and add an alias if the bare organ
name is what a reader would type. None of the four is derivable from the cancer
type, which is why an unregistered indication is refused rather than answered
with another one's data. `verify_indications.py` covers the multi-indication
properties.

---

## 6. The API surface

Two shapes on one port.

- **`/v1/cart/…`** — nineteen endpoints the client's dashboard expects. Served
  by an adapter that *translates and never decides*: it cannot reach a stage
  module and cannot spell a decision value, and both are enforced mechanically
  rather than asserted (`D1`, `D2`, each blinded).
- **`/projects/…`** — the platform's own surface, which the verifiers use.

**Fields are honoured or refused by name.** `objective` and `delivery_mode`
return 400 naming the field, because nothing reads either and accepting them
would tell a client an instruction was followed when it was discarded. Six
endpoints answer with a refusal and its reason rather than a result; those are
complete answers, not stubs.

**Client references are addressable aliases** resolved before any view is
called, so one id space has two names and every view sees one kind of id.

---

## 7. The binder benchmark and its blind rule

`run_benchmark.py` produces a blind output and freezes it; `compare_benchmark.py`
compares it against the literature panel and **refuses to run unless the frozen
output was committed before the answers were.**

Four layers, each blinded: the import graph is walked by parsing rather than
importing, held-out values are scanned for in every reachable module, the output
is fingerprinted, and the commit graph is the evidence of ordering.

---

## 8. Where the specifications live

Everything is in `specs/`. The ones to read first:

| | |
| --- | --- |
| `verification-sharing-assumptions.md` | **read this before writing a criterion** (§9) |
| `design-decisions.md` | every decision taken, with what it would cost to reverse |
| `stage3-target-discovery.md` | the six-component scoring frame and the risk model |
| `stage11-candidate-scoring.md` | the Level B frame: nine weights, three states, the floor |
| `stage4-target-pairing.md` | pairs, coverage, and the span confounding |
| `dashboard-surface-and-msln-benchmark.md` | the adapter contract and the benchmark |
| `decision-binder-count-objective.md` | **open, needs an answer** |
| `p0-`, `p1-`, `p2-`, `p13-` | the priority items, each closed with its evidence |

---

## 9. Failure shapes

The most useful thing this project knows. Two families, and a new reader will
reintroduce one within a week without this section.

### Family one: criteria that cannot fail

Thirteen instances are recorded in `specs/verification-sharing-assumptions.md`.
The check that would have caught every one is a single question:

> **What would this criterion report if the thing it tests were broken?**
> If the answer is "the same as it reports now", it is not a test.

The recurring shapes: a summary literal that hides the denominator; a criterion
clearing on an empty set; grepping source instead of executing it; a correct
check wired to the wrong input; a criterion recomputing a value by calling the
function that produced it.

**The strongest instance is the most recent, and it is worth knowing by name.**
`BM4` was written to prove that the benchmark refuses to freeze an output while
the answers are on disk. It wrote a probe to the real answers path and unlinked
it in a `finally` — so **every run of the benchmark verifier deleted the
literature panel it existed to protect**, and the next `git add -A` committed
the deletion, taking it out of the index. It survived only in history. A
criterion that destroys the evidence it was written to check is the failure this
family produces at its worst: it passed, every time, while removing the thing
that made passing mean anything.

### Family two: unexamined assertions

A different failure, and it deserves its own name because the question above
does not catch it. These are not computations that go wrong. **They are claims
printed on every run that nothing ever checked.** No criterion covers them
because no one thought of them as claims at all.

Four from a single session:

- **`run_all.py` ran thirteen verifiers and called it a full run.** Two existed
  and were never added to the stage list. Every total reported externally was
  over a subset, short by 33 criteria.
- **"raw caches read, not rebuilt" was printed on every run and was wrong about
  one directory in eleven.** The trials cache is keyed by the screened set and
  rewrote itself whenever that changed. Nothing compared before to after until
  a fingerprint was added, and then it showed up on the first run.
- **Two verifiers sat outside the inventory every audit read** while counting
  toward the suite total, because they used a different check idiom. The total
  meant two things at once.
- **`BM4` deleted its own evidence** (above). It belongs to both families: the
  criterion could not fail, *and* nothing ever asserted that running a verifier
  leaves the repository as it found it.

The general form: **a statement the system makes about itself, repeated often
enough to look verified, that no check has ever touched.** They are found by
asking a different question from the one above:

> **What does this line claim, and what would notice if it were false?**

Every one of the four was found by making the system prove a claim it had been
making for months — count the verifiers, fingerprint the caches, inventory the
criteria, hash the file before and after.

### What the audit found, and what it did not cover

A line tracer ran every verifier and counted how many times each criterion's
comparison body executed. **173 criteria traced; none cleared over ground it
never touched.** Eighty-two examined a non-empty collection, from 3 items to
317,712. Two reported unreached were alternate call sites — empty-set guards that
trip rather than clear, correctly dormant.

**The limitation matters as much as the result, and belongs here rather than in
a footnote.** The tracer counts comprehension bodies. A condition built by a
helper function reads as a "direct comparison" even when the helper examines a
collection inside itself. **Eighty-nine criteria fell into that bucket, and
"89 direct" does not mean "89 audited"** — it means the tracer could not speak
to them. The clean result covers the 82 it could count. The remaining 89 are
unexamined by this method and are exactly where a fourteenth instance would
hide.

---

## 10. What is not built, and why

Three different reasons, and conflating them is how a fundable gap gets treated
like an impossible one.

### Needs a decision

- **`binder_count` counts database hits, not target-matched binders.** 83 of 315
  structure-route binders are annotated against a different protein. One shipping
  design is affected, and counting matches rather than hits moves it from
  ADVANCE to BACKUP. Specified in `decision-binder-count-objective.md`, open.
- **The staining veto gates on presence rather than amount** (`R14`). A tolerance
  decision, priced, not taken.
- **`architecture_mode` values `OR` and `AND-NOT`** have no implemented routing
  row and are refused by name rather than silently substituted.

### Needs data nobody has connected

- **Stage 8, functional prediction.** Its training inputs are partner-generated
  experimental measurements. This is the one gap a decision alone cannot close.
- **Affinity.** `UNKNOWN` on all 422 retrieved candidates; no connected release
  carries the column. Fold error and rank correlation are therefore
  *structurally uncomputable*, not failed.
- **Cross-reactivity screening.** The proteome cache carries no sequence column
  and no alignment library is installed. Costed in the gaps section of every
  candidate package: one field on an existing query, plus a k-mer screen. What
  it would produce is a similarity screen, which is not a cross-reactivity
  measurement, so it stays unbuilt rather than shipping under that name.
- **Stage 13 and 14**, experimental recommendation and the learning loop. Both
  wait on the first ingested result.

### Needs hardware

- **Stage 7, structural evaluation.** No coordinate file is connected anywhere.
  Retrieval over deposited complexes is buildable and cheap; *prediction* is
  not — a single-chain folding model's weights alone exceed the available GPU
  memory, and the database-driven route exceeds available disk. Scoped in full
  in the Stage 7 report.

**One thing that reads like a gap and is not.** Developability emits no score.
Stage 10 counts liability flags and refuses to sum them, because a flag that
fires on every binder in the pool carries no information and a total would hide
that. That is a decision, recorded, and it is not waiting on anything.
