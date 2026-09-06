# The dashboard surface and the MSLN binder benchmark

Two documents arrived together and are specified together, because three of the
endpoints the dashboard is missing *are* the benchmark. Building the benchmark
delivers `binders/rank`, `binders/{binder_id}` and most of `binders/retrieve`.

Written before any code. Every criterion in §7 is fixed before the first run,
and every number this document predicts is stated so the run can be checked
against it rather than explained afterwards.

## 1. What this is, and what it is not

**It is** an adapter that renames and re-routes what the pipeline already
produces, three endpoints that return an honest constant, and a binder-ranking
benchmark over connected evidence.

**It is not** a second implementation of anything. The adapter decides nothing.
The benchmark introduces no source that is not already connected, and it
produces no affinity value, because none exists to produce.

---

# PART A — the three honest constants

`structure/evaluate`, `function/predict` and `learning/calibrate` do not exist,
and the dashboard's own example responses for all three are refusals:
`UNKNOWN`, `INSUFFICIENT_VALIDATED_DATA`, `BLOCKED_INSUFFICIENT_DATA`. An
endpoint that returns exactly that is honest and complete, not a stub.

| endpoint | response | why this is the true answer |
| --- | --- | --- |
| `POST /v1/cart/structure/evaluate` | `structure_status: "UNKNOWN"`, `structure_score: null` | No coordinate file exists anywhere in the cache. The structure source returns entry identifiers, chain labels and experimental methods, never atoms |
| `POST /v1/cart/function/predict` | `status: "INSUFFICIENT_VALIDATED_DATA"`, every prediction `null` | Stage 8's required training inputs are partner-generated experimental measurements, and none is connected |
| `POST /v1/cart/learning/calibrate` | `status: "BLOCKED_INSUFFICIENT_DATA"`, `required_action: "ADD_EXPERIMENTAL_RESULTS"` | No experimental result has ever been ingested; there is nothing to calibrate against |

**Each carries its reason in the response**, not only its status. A dashboard
rendering `UNKNOWN` with no explanation teaches a reader that the platform is
broken; one rendering `UNKNOWN` beside *no coordinate source is connected*
teaches them what to fund.

**These must never become computed.** Criterion D8 asserts each returns its
declared constant, so that the day one of them starts returning a number,
something deliberate has happened.

---

# PART B — the adapter

## 2. The constraint, and how it is enforced

**The adapter translates and never decides.** No gate logic, no ranking, no
vocabulary of its own. Path mapping, field renaming, and batch/single fan-out.
Anything it cannot do by translation is a gap in the pipeline and is built
there.

Stated as intention this is worthless — this repository has been burned by
criteria that assert good behaviour. It is enforced by two mechanical checks:

**D1 — the import graph.** The adapter module may import the view functions in
`api/server.py` and nothing else. It may not import `stage3`…`stage12`,
`scoring`, `routing` or `validation`. A walk of its import graph asserts this.
An adapter that cannot reach a decision-making module cannot make a decision.

**D2 — the vocabulary.** No value from the decision vocabulary
(`ADVANCE`, `BACKUP`, `VALIDATE`, `REQUIRES_EVIDENCE`, `EXCLUDED`), the gate
tokens, or the design classes may appear as a literal in adapter source. If the
adapter needs to name a decision, it is deciding; if it is only passing one
through, it never needs to spell it.

Both are blinded before being trusted: inserting an import of `stage11` must
trip D1, and inserting the literal `"ADVANCE"` must trip D2.

## 3. The path map

Nineteen endpoints under `/v1/cart/`, mapped to what exists.

### Translate directly — 8

| dashboard | pipeline | translation |
| --- | --- | --- |
| `POST /v1/cart/projects` | `POST /projects` | field rename; §4 on ids |
| `POST /v1/cart/runs` | `POST /projects/{id}/runs` | project id moves body → path |
| `GET /v1/cart/runs/{run_id}` | `GET /jobs/{id}` | stage name passed through; see §5 on `progress` |
| `POST /v1/cart/evidence/build` | `GET /projects/{id}/evidence/{gene}` | fan out over `target_ids[]`, collect |
| `POST /v1/cart/targets/discover` | `GET /projects/{id}/targets` | `max_targets` → `limit` |
| `POST /v1/cart/targets/pair` | `GET /projects/{id}/pairs` | `max_pairs` → `limit`; see §5 on `logic` |
| `POST /v1/cart/gates/evaluate` | rows from `GET /projects/{id}/contract` | select by `candidate_ids`; see §5 on `failed_gates` |
| `POST /v1/cart/candidates/rank` | `GET /projects/{id}/contract` | closest match of the nineteen; near-identical already |

### Need a shape change in the pipeline — 4

These are pipeline gaps, not adapter work, and are built in the pipeline:

| dashboard | what is missing | where it is built |
| --- | --- | --- |
| `POST /v1/cart/binders/retrieve` | no endpoint, no `binder_id`, no `target_match` | `binder_id` and `target_match` come from B2/B3 |
| `GET /v1/cart/binders/{binder_id}` | the whole scorecard | B3 + B4 + B5 |
| `POST /v1/cart/constructs/generate` | no `construct_id`; ours returns a run's constructs, not one built from (target, binder) | Stage 6 gains a stable id |
| `POST /v1/cart/manufacturability/evaluate` | ours is headroom/budget; theirs is differently defined | the definition is declared, not changed to match an illustrative number |

### Return an honest constant — 3

Part A.

### Remain unbuilt and named — 4

`POST /v1/cart/binders/optimize`, `POST /v1/cart/binders/rank`,
`POST /v1/cart/experiments/results`, and the ranking half of
`POST /v1/cart/experiments/recommend`. Of these, `binders/rank` is built here
as B10. The other three are named in §6.

**D9 asserts every one of the nineteen documented paths answers.** A documented
path that 404s is worse than one that refuses, because the client cannot tell
"not built" from "wrong URL".

## 4. Ids: one space, two ways to name it

The dashboard uses the client's own `WM-CART-PDAC-001` as the project id and
passes it to `/runs`. The platform generates a twelve-hex id and keeps the
client's string as `client_reference`.

**Client references become addressable aliases. There is one id space and two
ways to name a project in it.** No reverse-lookup layer, no second registry:
the route resolves an alias to the canonical id before any view is called, so
every view continues to see exactly one kind of id.

Two refusals keep the space unambiguous:

- A client reference matching `^[0-9a-f]{12}$` is **refused by name**. It would
  be indistinguishable from a canonical id, and resolution order would decide
  which project a caller reached.
- A client reference already in use is **refused by name**, with the existing
  project id reported. Silently returning the older project would hand a caller
  someone else's run.

**D6** requires an alias and its canonical id to return byte-identical bodies
from every view. **D7** requires both refusals.

## 5. Three places the dashboard's shape and ours genuinely differ

Named rather than papered over, because each is a real difference:

**`progress: 0.55`.** We report stage names, not a fraction. A fraction implies
a measured proportion of work done; what exists is an ordered stage list. The
adapter reports `stage` and `stages_complete`/`stages_total`, and sets
`progress` to `null` with a reason. Inventing 0.55 would be the adapter
deciding.

**`failed_gates: []`.** Stage 11 records the *first* gate a candidate failed and
stops, because attribution to one gate is what makes the attrition chain sum to
the pool. A list implies every gate was evaluated. The adapter emits a
one-element list and states that gating is first-failure attribution.

**`logic: "OR"`.** The dashboard's pair request offers `["AND","OR"]`. No OR
architecture is implemented; routing has no such row. `OR` is refused by name,
as it already is on `architecture_mode`.

## 6. The input contract stands

`objective` and `delivery_mode` keep returning **400 naming the field**.
Accepting and dropping them tells a client their instruction was followed.

`binder_mode` is **honoured for `RETRIEVAL_FIRST` only**, because that is
exactly what Stage 5 does: a structure route over deposited complexes and a
sequence route over named therapeutics, and nothing that generates a binder.
Any other value is refused by name, listing what is supported. Honouring a
value the platform genuinely implements is not the same as accepting one it
ignores. **D4 and D5** assert both halves.

---

# PART C — the MSLN binder benchmark

## 7. The blind rule, enforced in four layers

The benchmark's own critical rule: SS1, M5, M11, 15B6/h15B6, their affinities
and the expected CAR-suitability interpretation must not enter the algorithm.
Four layers, each mechanical:

**Layer 1 — the literal guard.** No held-out literal (`SS1`, `M5`, `M11`,
`15B6`, `h15B6`, `26.9`, `64.7`) may appear in any module reachable from the
benchmark entry point. A source scan over the walked import graph.

**Layer 2 — the import graph.** The answers module must be unreachable from the
ranking entry point. Asserted by walking imports, not by convention.

**Layer 3 — the freeze, which is the one that matters.** The blind run writes
its output with a configuration hash and **that artifact is committed before the
answers file is read**. The comparison refuses to run if re-executing the
algorithm does not reproduce the frozen hash. This makes the ordering a property
of the commit graph rather than of anyone's discipline, which is what the
pattern file already prescribes: *if the spec commit and the result commit are
the same commit, the exercise has already failed.*

**Layer 4 — blind the guard.** Inject a held-out literal into a reachable module
and confirm layer 1 trips; make the answers module importable and confirm layer
2 trips; perturb the frozen artifact and confirm layer 3 refuses. A guard that
has never failed is not known to work.

**One leak already avoided, recorded as precedent.** While scoping, source
coverage for MSLN was counted *without* looking up whether the four named
binders are present — 7 distinct entries, 14 antibody instances, 4 named
therapeutics. Knowing which entries correspond to the held-out names, and then
designing retrieval, would tune to the answer. That check belongs inside the
blind run.

## 8. B3 — antigen match without coordinates

§5.2 is built on this platform's own GPR35 finding, where a structural hit's
recorded antigen was a different protein. B3 is the check that catches it, and
it is decidable **because it asks a weaker question than the one previously
refused**.

Refused: what does this antibody *bind*? That needs coordinates.
Asked here: what did the depositors *record* as its antigen? That is annotation.

The mechanism: match the deposited entry's `antigen_name`, `antigen_type` and
`antigen_species` against a name set built from the target's own UniProt
record — recommended name, parenthesised synonyms, and Chain feature notes.

For MSLN that set is derived, not written down: *Mesothelin; CAK1 antigen; MPF;
Megakaryocyte-potentiating factor; Pre-pro-megakaryocyte-potentiating factor;
Mesothelin, cleaved form.*

**Predicted verdicts, fixed before the run:**

| entry | recorded antigen | type | expected |
| --- | --- | --- | --- |
| `4F3F` | Mesothelin | PROTEIN | `PASS` |
| `7U8C` | Mesothelin, cleaved form | PEPTIDE | `PASS` — matches Chain `PRO_0000253561` verbatim |
| `8H8J` | G-protein α-13 \| G(I)/G(S)/G(T) β-1 | PROTEIN | `FAIL`, `wrong_antigen_flag` |
| `1P4B` | GCN4(7P-14P) peptide | PEPTIDE | `FAIL`, tag |

**These pins are on antigen annotation, not on binder identity.** Which
deposited entry corresponds to which literature-named binder is a held-out
answer and must be produced by the blind run. Pinning that an entry's recorded
antigen is Mesothelin says nothing about which named binder it is.

Two design constraints, both from the data:

- **`antigen_species` is a soft signal, never a hard gate.** It is `NA` on two
  of the four entries above, and gating on it would fail 7U8C — a genuine MSLN
  entry. **BM7** asserts 7U8C survives.
- **Missing or ambiguous annotation yields `UNKNOWN`, never `PASS`.** This is
  where coordinates would remain the only resolution, and the output says so.

The catch on 8H8J works by **absence**: GPR35 appears nowhere in its own entry's
antigen annotation. Note also that a substring test on gene symbols would be the
wrong mechanism — it fails on MSLN, whose antigen is recorded as "Mesothelin".

## 9. B4 — sequence sanitation

Original and cleaned sequences are stored separately and **provenance is never
overwritten**. The worked case is already in the cache: the anti-tag binder's
chains carry `GGGGSGGGGSGGGGSGGGGS` and `ASGADHHHHHH`.

`sequence_original` is immutable. `sequence_clean` exists only where source
metadata justifies the removal, and every entry in `modifications[]` names the
motif, the span, and the metadata that justified it. Where a motif is detected
but nothing justifies removing it — a histidine run that may be germline-encoded
— it is **flagged and not cleaned**. **BM8** asserts both halves.

## 10. B9 — structural evidence, gated on B3

The benchmark's own rule: use deposited structural evidence only if target match
is verified, otherwise `UNKNOWN`/`null`. That is exactly the Stage 7 design
already scoped.

What B9 returns today: the verified entry, its method and its resolution — for
MSLN, 14 instances, all XRAY, **1.52 to 4.31 Å**. What it does not return:
interface geometry, which needs coordinates that are not connected. **BM9**
asserts no structural evidence is emitted for a binder whose B3 verdict is not
`PASS`.

**Resolution is carried per entry, never averaged, and it is not a pass/fail.**
The span is wide enough to matter: a 1.52 Å structure and a 4.31 Å structure
support different claims about the same thing, and a stage reporting only *a
verified complex exists* would flatten that. No resolution cutoff is imposed,
because choosing one after seeing this distribution would be a threshold fitted
to an observation.

## 11. B10 — two rankings, and the front as a third output

The document requires binding rank and CAR-suitability rank as separate outputs
and is explicit that highest affinity must not equal best CAR binder. The
existing machinery applies, with one correction.

`scoring.py` already implements what §6.2 describes almost verbatim: a weighted
sum over normalised components, **normalised over the measured subset**, a floor
below which no score is emitted, three states where missing is never a
favourable zero, a confidence adjustment kept **outside** the sum, and hard
gates no score can rescue.

**The correction:** the document's two rankings are two different *weighted*
questions. The platform's two existing outputs are a weight-free view (the
Pareto front) and one weighted score — not the same axis. So:

| output | mechanism | component set |
| --- | --- | --- |
| **binding rank** | `scoring.combine`, weight version `binder_binding_v1` | target match, binding evidence, structural verification, evidence confidence |
| **CAR-suitability rank** | `scoring.combine`, weight version `binder_suitability_v1` | epitope accessibility, membrane proximity, shedding suitability, specificity, developability, humanness, affinity |
| **the front** | `pareto_front`, unchanged | the union, weight-free |

**The front is retained as a third output because it is the strongest guarantee
of the document's own rule.** Front membership cannot be moved by either weight
set — criterion W7 already asserts that property for the candidate front — so
where binding rank and suitability rank disagree, the front shows which binders
are non-dominated across both. A reader who disagrees with either weight set
still has it.

### The two weight sets, declared before any run

Both sum to 1.00 exactly and are versioned. Neither is fitted: each follows the
question its ranking asks.

**`binder_binding_v1`** — *which binder appears to bind the target most
strongly?* Affinity is the direct answer, so it carries the most weight. The
other two are evidence that binding occurs at all, and evidence that the first
two can be interpreted.

| component | weight | why |
| --- | --- | --- |
| affinity | **0.55** | the direct measure of the thing this ranking asks about |
| structural verification | **0.30** | a deposited complex with the B3-verified antigen is experimental proof that binding occurs |
| evidence completeness | **0.15** | assay, format and resolution metadata, without which the other two cannot be read |

**`binder_suitability_v1`** — *which binder makes the best CAR?* Affinity is
deliberately **not** the largest weight. That is the document's own thesis, and
putting affinity on top here would contradict §6.2 in the act of implementing
it.

| component | weight | why |
| --- | --- | --- |
| epitope accessibility | **0.20** | a receptor that cannot reach its epitope does not work at any affinity |
| membrane proximity | **0.15** | synapse geometry, and the reason a juxtamembrane epitope is a distinct question |
| shedding suitability | **0.15** | soluble antigen competing for the receptor is a CAR-specific failure mode |
| specificity | **0.15** | on-target off-tumour risk carried at binder level |
| affinity | **0.15** | it matters, and it is one factor among several rather than the determinant |
| developability | **0.12** | manufacturable at all, given the hard gates already passed |
| humanness | **0.08** | lowest because what is connected is a naming convention, not a sequence measurement |

### Predicted arithmetic, fixed before the run

In the same form as the candidate frame's 0.96 / 0.60 / 0.625, which landed
exactly. Every component is `UNKNOWN` unless something connected measures it.

| ranking | binder route | applicable | measured | fraction | floor 0.50 |
| --- | --- | --- | --- | --- | --- |
| binding | structure route, B3 `PASS` | 1.00 | **0.45** | **0.45** | below → `null` |
| binding | sequence route | 1.00 | **0.15** | **0.15** | below → `null` |
| suitability | either route | 1.00 | **0.00** | **0.00** | below → `null` |

**Predicted outcome: 0 binders receive a score in either ranking, and the Pareto
front is empty**, because no component carries a value for the front to compare.

The reasons, component by component, so the prediction is checkable rather than
atmospheric: affinity is `UNKNOWN` on all 422 candidates because no connected
release carries the column; epitope accessibility and membrane proximity need an
epitope, which needs coordinates; shedding suitability is derivable at *target*
level for MSLN — the `296–598` cleaved chain is exactly that — but not per
binder without an epitope; specificity is `UNKNOWN` because B8 stays unbuilt;
developability is `UNKNOWN` by the standing decision that Stage 10 does not sum
its flags into a score, which this document does not overturn; humanness is
`UNKNOWN` because no germline reference is connected.

**What the run does produce**, and what the report leads with after the coverage
statement: retrieval counts with source coverage, B3 verdicts per entry,
original and cleaned sequences with every modification justified, and structural
evidence carrying its own resolution. Those are the results. The rankings being
`null` is a fourth result, not a failure to produce the first three.

**If any of this lands elsewhere it is the first thing reported**, before any
interpretation. A ranking that comes back populated means a component was
measured that this table says is `UNKNOWN`, and that is a finding about the
sources, not a success.

One consequence worth stating plainly for whoever reads the dashboard: **a
correct first run renders an empty ranking table.** That is the platform
declining to order binders on evidence it does not have, and it is the same
refusal that makes the candidate frame trustworthy.

## 12. B1, B5, B7 — build what is derivable, name what is not

| | derivable now | not derivable, named |
| --- | --- | --- |
| **B1** target normalisation | accession, synonyms, GPI anchor, TM count, topology, and Chain boundaries. For MSLN: `37–598 Mesothelin`, `37–286 Megakaryocyte-potentiating factor`, `296–598 Mesothelin, cleaved form` — **the shedding biology, encoded** | the residues themselves; the proteome cache has no sequence column. Epitope regions |
| **B5** antibody annotation | VH/VL assignment, already present as chain identifiers and heavy/light sequences | CDR boundaries beyond a stated approximate method — no numbering tool is installed. **Germline similarity and humanness: no reference is connected**, and the naming-convention signal already held is a convention, not a measurement |
| **B7** developability | pI, net charge, aggregation windows, glycosylation sequons, cysteine parity, plus oxidation, deamidation, isomerisation and complexity — all motif or entropy calculations over a sequence, no new source | nothing further. Note Stage 10's standing refusal to sum flags into a score is **not overturned**: the scorecard lists flags and carries no total |

## 13. What stays unbuilt, and what connecting it would cost

**B6 — affinity.** Unbuilt because there is nothing to build from. Verified:
`affinity` is the literal `NOT_CONNECTED` on **all 422** retrieved candidates,
and neither connected antibody release carries an affinity, KD or ΔG column.
`EXPERIMENTAL_RETRIEVED` is unavailable for every binder; `PREDICTED` needs a
model that does not fit the available hardware.

**B8 — specificity proxy.** Unbuilt, and the cost is worth stating because it is
small and specific:

1. **One field on an existing query.** The proteome cache holds 20,431 reviewed
   human entries with no sequence column. Adding the sequence field to the
   existing UniProt query is one line and roughly doubles that cache.
2. **An alignment method, which is the real cost.** No alignment library is
   installed. The CPU-friendly route is a k-mer identity screen in the numeric
   library already present — feasible, and it yields a *similarity screen*, not
   a cross-reactivity measurement. The output must be named as the former.

Until both land, B8 returns `specificity_risk: null` with the reason.

**`binders/optimize`.** Unbuilt. Stage 5 retrieves and stops; nothing generates
a binder, models a complex, optimises an interface or germlines a framework.
The `OPTIMIZED_VARIANT` and `CPU_GENERATED_VARIANT` origins the document names
have no producer.

## 14. What the output states up front

Not in a footnote, not in a note at the end — **in the response body and at the
top of the report**:

> Affinity is `UNKNOWN` on all 422 retrieved binder candidates. No connected
> evidence release carries an affinity, KD or ΔG column. §5.4's fold-error
> calculation and §7's Spearman rank correlation are therefore **structurally
> uncomputable, not computed and failed** — there is no value on either side of
> the comparison. §5.1 permits exactly this, and it should be read as a
> source-coverage fact rather than as a result.

**BM12** asserts this statement is present and that no affinity is ever reported
as anything but `UNKNOWN` while no source is connected.

---

# 15. Rejection criteria, fixed before any run

Each is written to fail. Each is blinded before it is trusted.

| | trips when |
| --- | --- |
| **D1** | the adapter's import graph reaches any stage module, `scoring`, `routing` or `validation` |
| **D2** | any decision, gate-status or design-class literal appears in adapter source |
| **D3** | any field in an adapter response has no corresponding field in the pipeline response it translates |
| **D4** | `objective` or `delivery_mode` is accepted, or refused without naming the field |
| **D5** | `binder_mode: RETRIEVAL_FIRST` is refused, or any other value is accepted |
| **D6** | an alias and its canonical id return different bodies from any view |
| **D7** | a client reference shaped like a canonical id is accepted, or a duplicate reference is accepted |
| **D8** | any of the three constant endpoints returns a computed value |
| **D9** | any of the nineteen documented paths returns 404 |
| **BM1** | a held-out literal appears in any module reachable from the benchmark entry point |
| **BM2** | the answers module is reachable from the ranking entry point |
| **BM3** | re-running the algorithm does not reproduce the frozen artifact's hash |
| **BM4** | the comparison runs against an artifact not committed before the answers were read |
| **BM5** | any of layers 1–3 fails to trip when deliberately blinded |
| **BM6** | B3 returns `PASS` for an entry whose recorded antigen is absent from the target's name set, or `FAIL` for one that matches, or `PASS` where annotation is missing |
| **BM7** | `antigen_species` is used as a hard gate — asserted by requiring 7U8C, whose species is `NA`, to survive |
| **BM8** | `sequence_original` is modified, or a sequence is cleaned without a recorded justification |
| **BM9** | structural evidence is emitted for a binder whose B3 verdict is not `PASS` |
| **BM10** | either ranking emits a score for a binder that failed a hard gate, or any component is imputed |
| **BM11** | the front's membership changes between the two weight versions, or the two rankings read the same component set |
| **BM12** | the affinity-coverage statement is absent from the output, or any affinity is reported as other than `UNKNOWN` while no source is connected |

## 16. Numbers predicted before the run

Checked first, before any interpretation:

| quantity | predicted |
| --- | --- |
| retrieved candidates carrying an affinity value | **0 of 422** |
| MSLN entries whose recorded antigen matches the UniProt name set | **7 distinct, 14 antibody instances** |
| MSLN-targeting named therapeutics carrying a heavy-chain sequence | **4** |
| B3 verdicts | `4F3F` PASS, `7U8C` PASS, `8H8J` FAIL, `1P4B` FAIL |
| binding rank, structure route | applicable 1.00, measured 0.45, fraction 0.45 → `null` |
| binding rank, sequence route | applicable 1.00, measured 0.15, fraction 0.15 → `null` |
| CAR-suitability rank, either route | applicable 1.00, measured 0.00, fraction 0.00 → `null` |
| binders receiving a score in either ranking | **0** |
| Pareto front over binders | **empty** |

If any lands elsewhere, that is the first thing reported, before the biology.

## 17. Order of work

1. This document, reviewed and committed before any code.
2. The three honest constants, with reasons in the response.
3. The adapter: the eight that translate, then the four shape changes in the
   pipeline; D1–D9 land with it.
4. The blind-freeze harness and its four layers, **before** any benchmark
   module, so nothing is written under an unenforced rule.
5. B3, B4, B9 — derivable, and B3 is the check §5.2 is built around.
6. B10: two weight versions, the front as third output.
7. B1, B5, B7 partial — building what is derivable, naming what is not.
8. The report, leading with the affinity-coverage statement.

## 18. What this does not license

Weights fitted to the benchmark outcome. A component imputed when absent. An
affinity value invented to make a correlation computable. A ranking emitted
below the floor. The adapter deciding anything. Reading the answers file before
the frozen artifact is committed. Each has a criterion above.
