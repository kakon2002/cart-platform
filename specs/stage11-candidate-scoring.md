# Stage 11 — the candidate scoring frame (§14.2, §14.3, §16)

Written before any implementation. The criteria in §11 are fixed here and
committed before any result.

Seven of the eleven ranking components do not exist. **That is the reason to
build the frame, not to wait.** A frame with four components filled and seven
named as absent is a true statement about the platform. Four components silently
carrying the whole score is not, and is what happens if the frame is deferred
until the components arrive.

---

## 1. What this is, and what it is not

This adds a **Level B score** beneath the hard gates that already exist. It does
not add a gate, does not move a threshold, and does not change which candidates
survive. Every candidate that reaches scoring has already passed Level A.

It is not a calibration. The reference document says weights "should be
calibrated prospectively as experimental data accumulate", and **there is no
outcome data to calibrate against**. A weight fitted today would be fitted to
nothing, which is the defect this repository has recorded twice — a bound chosen
after seeing the value it bounds is always satisfied by that value. Every weight
below is therefore **reasoned from a stated ordering rule and declared before any
score exists**, and the rule is written down so the weights can be argued with
rather than merely accepted.

## 2. Eleven components, nine of which are weights

The document's §14.1 table lists eleven ranking components, and its §14.2 formula
also carries a separate `confidence_adjustment`. Read literally, evidence
confidence and prediction uncertainty appear in both — once inside the sum and
once as the multiplier — and would be counted twice.

They are the multiplier, not summands. **Nine components carry weights; two
parameterise the adjustment.** This is not a reduction of the document's list:
every one of the eleven carries a declared, versioned coefficient, and none is
dropped. It is where each coefficient acts that differs.

The reason is the platform's standing rule, which this preserves: **normal-tissue
risk and evidence confidence are never combined into one number.** Putting
confidence inside the same sum as the safety term is exactly that combination — a
well-evidenced dangerous candidate and a poorly-evidenced safe one would reach
the same total. Keeping confidence outside as a multiplier scales the whole
score by how well founded it is, without ever letting confidence substitute for
safety.

## 3. The declared weights

**The ordering rule, stated before the numbers: measured outranks modelled, and
decisiveness outranks granularity.** Every component derived from a model sits
below every component derived from a measurement. Among measurements, a
component that changes whether the design can work at all sits above one that
changes how well it is made.

| # | component | weight | why it sits here |
| --- | --- | --- | --- |
| 1 | Tumour coverage | **0.18** | How much of the tumour the design reaches. The single largest determinant of whether a response is possible at all |
| 2 | Malignant-cell specificity | **0.16** | Whether the antigen is on tumour cells or on stroma. A design aimed at fibroblasts is a different product, not a worse one |
| 3 | Normal-tissue safety | **0.16** | Residual margin below the ceiling. Among gate-passers this is margin rather than pass or fail, and margin is what separates a design you would dose from one you would watch |
| 4 | Binder quality | **0.12** | No binder, no product. Among candidates that have one, its quality drives everything downstream of it |
| 5 | Manufacturability | **0.10** | A design that cannot be made consistently is not a candidate — but the hard packaging gate has already removed the impossible ones, so this is degree, not kind |
| 6 | Developability | **0.10** | The same argument one level more granular: sequence liabilities within a construct that already fits |
| 7 | Structural feasibility | **0.08** | Geometry matters, and it is modelled. The first component below the measured/modelled line |
| 8 | Functional prediction | **0.06** | The most model-dependent term in the set, and the document itself admits it is valid only inside a stated applicability domain |
| 9 | Pairing robustness | **0.04** | Applies only to multi-target designs. Lowest weight because for most candidates it does not apply at all, and an inapplicable component should not carry weight it can rarely spend |

**Σ = 1.00 exactly.** Criterion W1 asserts it.

The two adjustment coefficients:

| # | component | coefficient | effect |
| --- | --- | --- | --- |
| 10 | Evidence confidence | multiplier, exponent **1.0** | scales the score linearly by how well founded it is |
| 11 | Prediction uncertainty | penalty weight **0.5** | maximal uncertainty halves the score; it never zeroes it, because a score with wide error bars is still a score |

`confidence_adjustment = evidence_confidence × (1 − 0.5 × prediction_uncertainty)`

Where prediction uncertainty is UNKNOWN — it is, until Stage 8 exists — the
penalty term is **not applied and its absence is recorded**. It is not treated as
zero uncertainty, which would be the favourable imputation the document forbids.

## 4. Three states, not two

Every component on every candidate is exactly one of:

- **MEASURED** — a value in [0, 1], with the stage that produced it named.
- **UNKNOWN** — no measurement exists. Named per candidate, with the reason, so
  a reader sees *which* components are missing rather than inferring it from a
  low score. This is point 3 of the instruction and it is the difference between
  a frame and a black box.
- **NOT_APPLICABLE** — the component has no meaning for this design. Pairing
  robustness on a single-target or adaptor design is the worked case: it is not
  missing evidence, it is a question that does not arise.

**The three are never collapsed.** UNKNOWN and NOT_APPLICABLE both keep a
component out of the score, and they say different things: one is a gap in the
evidence, the other is a gap in the question. Reporting them as one would make a
platform that has never built a dual look identically uncertain to one whose
pairing data is missing.

## 5. Normalisation over the measured subset

```
applicable   = Σ wᵢ  over components that are not NOT_APPLICABLE
measured     = Σ wᵢ  over components that are MEASURED
fraction     = measured / applicable
overall      = ( Σ (wᵢ × cᵢ) over MEASURED / measured ) × confidence_adjustment
```

**The denominator is the measured weight, not 1.** Dividing by 1 with four
components present would hand those four the weight belonging to eleven, and a
candidate scored on coverage and safety alone would be directly comparable to one
scored on all nine. It would not be. `fraction` is reported beside every score
and is the number that says how much of the frame the score rests on — the same
rule and the same reasoning as `MINIMUM_MEASURED_WEIGHT` in Stage 3, which is the
precedent this follows.

## 6. The floor

**`MINIMUM_SCORED_FRACTION = 0.50`.** Below it no overall score is emitted:
the candidate carries `overall: null` with the reason and the list of UNKNOWN
components. It is not ranked by score, and it is not silently given a low one.

Reasoned, not fitted: a score resting on less than half the applicable frame is
not a summary of the frame, it is a summary of whichever half happened to be
measured. Stage 3 draws the same line at 0.40 for a six-component set; this set
has nine and a heavier tail of modelled components, so it draws it slightly
higher.

**What that yields today**, stated in advance so the first run can be checked
against it rather than explained afterwards: coverage, specificity, safety and
manufacturability are available (0.18 + 0.16 + 0.16 + 0.10 = **0.60**); binder
quality is UNKNOWN for an adaptor design, which carries no Stage 5 binder;
developability is UNKNOWN **by standing decision**, because Stage 10 refuses to
sum its flags into a score and this specification does not overturn that;
structural and functional are UNKNOWN for want of Stages 7 and 8; pairing
robustness is NOT_APPLICABLE to a single-target adaptor. Applicable = 0.96,
measured = 0.60, **fraction = 0.625**, which clears the floor. The five current
candidates should therefore receive a score, with four components named UNKNOWN
and one NOT_APPLICABLE.

## 7. The Pareto front stays

Both are emitted. They answer different questions: the front says which designs
are not beaten on every axis at once, the score says which is best under a
declared set of priorities. A reader who disagrees with the weights still has the
front.

**The front is computed from the component values and never from the weights**,
so no choice of weights can move it. Criterion W7 asserts this by re-running the
whole ranking under a deliberately different weight set and requiring the front's
membership to be identical while the score order is free to change. That is what
makes "the front cannot be gamed by weights" a claim rather than an assurance.

## 8. Amendment — N4 re-specified

**N4 as it stands:** *no weighted or summed score across objectives is emitted.*

**N4 as re-specified:** *no weighted sum rescues a hard-gate failure, and no
component is imputed.*

The original rule was not about arithmetic. It was written because a weighted
sum lets a strong score on one axis buy off a weak one, and the axis at risk was
safety: a design with excellent predicted efficacy could out-total a safer design
and become Rank 1. Forbidding the sum was the available way to prevent that when
ranking and gating were one step.

**§14.1's two-level structure is what makes the sum safe, and it is the whole
reason this amendment is defensible.** Safety is a Level A hard gate. A candidate
over the applied ceiling never reaches scoring, so no weight can rescue it and no
efficacy term can compensate for it. What the score now differentiates is
residual margin among candidates that have already passed — which is the thing
N4 was never trying to prevent.

Stage 11 was already two-level: it attributes every candidate to the first gate
it failed and computes objectives only for survivors. The change is confined to
Level B, and the guarantee N4 protected is now carried by three things rather
than one — the gate order, criterion W4, and the retained Pareto front.

**What the amendment does not license.** Weights fitted to an outcome; a
component imputed when absent; confidence folded into risk; a score emitted below
the floor; or a gate relaxed because a candidate scores well. Each has its own
criterion.

## 9. §14.3 — the candidate table, and a collision that must be removed first

The table needs four things that do not exist: `candidate_id`, `gate_status`,
`overall score`, and a `decision`.

**The decision vocabulary collides with the design-class vocabulary, by
substring, in both directions:**

| decision | design class | relation |
| --- | --- | --- |
| `ADVANCE` | `ADVANCED` | `"ADVANCE" in "ADVANCED"` is true |
| `BACKUP` | `CONSERVATIVE_BACKUP` | `"BACKUP" in "CONSERVATIVE_BACKUP"` is true |

This repository has found the same shape three times — `renal` matching
**adrenal gland**, `cortex` matching **Kidney_Cortex**, `data` matching
**car_pipeline/data/**. In every case a substring match produced a larger answer
that still typechecked, and nothing errored. Shipping two vocabularies in this
relation, into adjacent fields on the same object, is laying that trap
deliberately.

**Proposed resolution: rename the design classes to `CONSERVATIVE_DESIGN` and
`INNOVATIVE_DESIGN`.** "Innovative" is the reference document's own word — *"one
advanced or innovative design"*. That leaves the two sets disjoint under equality
and under substring in both directions, and it keeps the externally-specified
decision vocabulary exactly as §14.3 names it.

**The cost, stated:** `design_class` appears in `validation.py`, the constructs
view, the candidate package, DEPLOY.md and two specs. It is a rename with a
visible blast radius and it changes a value the API already returns. The
alternative — keeping both and relying on everyone using `==` forever — is
cheaper today and is the bet this repository has lost three times.

**Applied, step 2 of §12, before any decision-column code.** `CONSERVATIVE` is
now `CONSERVATIVE_DESIGN` and `ADVANCED` is now `INNOVATIVE_DESIGN`. The constant
`ADVANCED` was renamed to `INNOVATIVE` alongside its value, and the summary's dict
keys `conservative_backup` and `advanced` became `conservative_design` and
`innovative_design` — a key named `advanced` returning `INNOVATIVE_DESIGN` is the
same drift in a different field. The captured payload in DEPLOY.md still shows
the old value and is recaptured by the step-8 run rather than hand-edited, per
the standing rule on captured output.

**The two are computed independently and neither is derived from the other.**
Design class answers *what kind of design is this*, from the architecture table,
fixed before the run. Decision answers *what should happen to this candidate*,
from gate status and rank. Criterion W9 asserts the independence by requiring
candidates of the same design class to receive different decisions where rank
differs.

Decision assignment, declared here:

| decision | when |
| --- | --- |
| `ADVANCE` | passed every gate and is **not dominated** — on the Pareto front |
| `BACKUP` | passed every gate and is dominated by some other survivor |
| `VALIDATE` | passed every gate but the measured fraction is below the floor, so no score was emitted |
| `REQUIRES_EVIDENCE` | failed a gate that a measurement could later clear — no binder, no construct |
| `EXCLUDED` | failed a gate that evidence will not change — over the safety ceiling, over the payload budget |

`candidate_id` is `CAR-{indication key}-{position:03d}`, matching the document's
`CAR-PDAC-001`. **It is a within-run label, not a durable identity**; the durable
identity of a candidate is its gene plus the configuration-hash chain, and the
spec says so rather than letting a sequential id be mistaken for a key.

### Amendment — ADVANCE is front membership, not rank position

The first draft defined `ADVANCE` as *ranked first among scored candidates* and
`BACKUP` as *ranked below first*. Both rest on a total order over survivors, and
reading the code before implementing showed that no such order is computed
anywhere in Stage 11.

**What `position` actually is.** `stage4.py:278` sorts the pool by
`-composite_supported`. `stage11.rank` appends rows in that order and never
sorts. `stage12.build` enumerates the survivors it is handed. So the package's
`rank 1 of 5` is Stage 4's composite ordering, inherited three stages downstream
and never declared as a ranking rule — while `_ranking_payload`'s docstring says
the position is *"carried from the ranking with no re-ordering"*, which reads as
though Stage 11 ordered it. Stage 11 did not.

Defining `ADVANCE` on that position would have made the platform's headline
decision a function of an undeclared inherited sort. It would also have
contradicted the payload sitting beside it, which says *"No weighted total
across objectives is emitted"* — `composite_supported` is a weighted total, so
ranking first by it is the weighted sum re-entering through the display column.
That is N4 being defeated by the field next to the one N4 guards.

**The resolution.** `ADVANCE` is now non-domination, which Stage 11 genuinely
computes and which needs no total order to be well-defined. More than one
candidate can be told to advance, and in the current pancreatic run two are
(FER1L6 and GPR35, the two front members). `BACKUP` is a survivor some other
survivor beats on every objective at once — a defensible second choice, which is
what the word means.

`position` survives as a display index and nothing decisional hangs off it. Its
basis is now stated in the payload rather than implied, so a reader is told the
order came from Stage 4's composite and not from the front.

**`VALIDATE` cannot be reached until step 4**, because it is defined on a
measured fraction below the floor and no score is emitted yet. That is a branch
no run can enter, which is the shape this repository has now recorded twelve
times. It is handled by exercising `decision_for` directly with the state the
pipeline cannot yet produce, so the branch is tested before it is reachable
rather than after.

## 10. §16 — accept what the document names, reject the rest loudly

Today `POST /projects` reads two keys and **silently discards everything else**.
Verified: a POST carrying all eight §13.1 inputs returns `201` and drops eight of
them without a word.

| §16 field | after this change |
| --- | --- |
| `indication` | accepted, alias of `cancer_type`, both honoured |
| `target_mode` | accepted; `DISCOVER` or a gene symbol, replacing the derived-from-presence rule |
| `project_id` | accepted as `client_reference` and echoed; the server-generated id stays canonical, because ids must be unique and server-controlled |
| `max_final_candidates` | **honoured** — a slice on the ranked list |
| `architecture_mode` | **honoured where it maps to `CARFormat`**, which already accepts `AUTO`; values with no mapping are rejected by name |
| `objective` | **400**, naming the field and saying it is not yet honoured |
| `delivery_mode` | **400**, same |
| any unknown field | **400** |

A field the platform does not honour is **rejected by name** with
`UNSUPPORTED_INPUT`, saying what the platform does instead and what to remove.
Accepting a field and ignoring it tells a client their instruction was followed.
It was not. A bare refusal is only marginally better, so the body names the field.

### Amendment — two fields honoured rather than refused

The first draft returned `400` for all four unhonoured fields. Two are now
honoured, for reasons that hold independently of convenience.

**`max_final_candidates` is honoured.** It is a slice on an already-ordered list
and needs no science the platform lacks. Refusing a field that the document's own
example JSON sends would read as the platform being unable to count. It caps the
ranked list only; it never manufactures candidates to reach the number, which is
§17's explicit requirement, and the response reports both the cap and how many
were actually eligible.

**`architecture_mode` is honoured where it maps.** `CARFormat` already accepts
`auto`, `conventional`, `dual_target`, `logic_gated`, `switchable` and `armored`,
and the document's `AUTO`, `SINGLE`, `AND`, `OR`, `AND-NOT`, `ADAPTOR` overlap
that vocabulary in part. The mapping is declared here:

| `architecture_mode` | maps to | note |
| --- | --- | --- |
| `AUTO` | `CARFormat.AUTO` | routing chooses, as today |
| `SINGLE` | `CARFormat.CONVENTIONAL` | single-antigen receptor |
| `AND` | `CARFormat.LOGIC_GATED` | the AND-gate row of the architecture table |
| `ADAPTOR` | `CARFormat.SWITCHABLE` | the adaptor row |
| `OR` | **rejected** | no OR-gate row is implemented; routing has no such architecture |
| `AND-NOT` | **rejected** | the inhibitory row is `NOT_IMPLEMENTED` in routing, with a recorded reason |

A rejected value is named in the refusal along with the values that do map.
**Accepting `OR` or `AND-NOT` and then routing something else would be the
silent-substitution failure this platform exists to avoid**, which is why they
refuse rather than falling back to `AUTO`.

`objective` and `delivery_mode` still return `400`. Neither has anything behind
it: no optimisation objective is read anywhere, and delivery modality is fixed by
the vector payload budget in Stage 1. Honouring them would mean accepting a value
that changes nothing.

**This is a breaking change** for any caller sending extra keys. None is known —
the two verifiers and DEPLOY.md send only `cancer_type` — and it is stated here
rather than discovered.

## 11. Rejection criteria — fixed before the run

| id | trips when |
| --- | --- |
| **W1** | the nine weights do not sum to 1.0 within 1e-12, or any of the eleven components has no declared coefficient, or the weight set carries no version |
| **W2** | any component on any candidate is not exactly one of MEASURED, UNKNOWN, NOT_APPLICABLE; or an UNKNOWN or NOT_APPLICABLE component contributes to either the numerator or the denominator; or an UNKNOWN component is not named with its reason |
| **W3** | any candidate's overall score does not recompute from the component values, weights and adjustment recorded on that candidate, to within 1e-12 — the record must reconstruct its own score |
| **W4** | any candidate that failed a hard gate carries an overall score or appears in the ranked list. The §14.1 guarantee, asserted rather than assumed |
| **W5** | any candidate whose measured fraction is below the floor carries a numeric score, or any candidate above the floor carries `null` |
| **W6** | the safety component changes when evidence confidence changes, or the confidence adjustment changes when the safety component changes. The two must be independently derived, which is the standing rule made falsifiable |
| **W7** | re-running the ranking under a different declared weight set changes the Pareto front's membership. The score order may change; the front may not |
| **W8** | any decision value is absent from the declared vocabulary, or any decision value equals or is a substring of any design-class value in either direction |
| **W9** | any candidate's decision does not recompute from its gate status and front membership alone, or the survivors split on front membership yet all carry the same decision |
| **W10** | `POST /projects` accepts a field the platform does not honour, or rejects one it does, or fails to name the offending field in the refusal |
| **W11** | the Stage 11 configuration hash does not change when a weight changes — a run under different weights must not compare equal to this one |

### Amendment — W9 tested the rank rule, not the property

W9 was written against the first draft's `ADVANCE`, which was rank position.
With `ADVANCE` now defined on non-domination, the original wording — *two
candidates sharing a design class and differing in rank receive the same
decision* — trips on correct output: FER1L6 at position 1 and GPR35 at position
2 differ in rank, share `INNOVATIVE_DESIGN`, and both correctly read `ADVANCE`
because both are on the front. A criterion that fails when the code is right is
worse than no criterion, because the run gets amended to satisfy it.

This is the stale-pin shape for the fifth time, caught before the run rather
than by it: a criterion holding a *result* — here, a specific rank-to-decision
mapping — rather than the property it exists to test. The property is that
decision and design class are computed independently, so W9 now recomputes each
decision from gate status and front membership with no reference to design class
and requires the survivors' decisions to distinguish them whenever the front
does.

**W8 and W9 land with the decision column, not at step 6.** They are that
column's criteria, and shipping the column ahead of them would be the run
preceding its criteria. W1–W7, W10 and W11 stay with the steps whose behaviour
they describe.

**W7 is the criterion that matters**, and it is the same shape as S11 in the
construct-safety arm: it makes a claim falsifiable by changing the input the
claim is about. W4 and W6 carry the two guarantees the amendment in §8 rests on.

### Explicitly not grounds for rejection

- Seven components being UNKNOWN today. That is the finding, named per candidate.
- The five current candidates scoring closely. They are the same architecture
  with the same binder; a frame that separated them sharply would be suspect.
- The weights being reasoned rather than fitted. §1, and there is nothing to fit
  them to.

## 12. Order of work

Ordered so that the part a frontend renders lands before the part gated on
Stages 7 and 8. If the milestone runs short, the decision column and gate status
are worth more than the score, and the score is the piece whose components are
missing anyway.

1. This document, reviewed, committed before any code.
2. **The design-class rename. Blocking** — it must land before the decision
   vocabulary, not after, or the collision ships and has to be unpicked.
3. `stage11`: `candidate_id`, `gate_status`, and the decision column.
4. `stages/scoring.py` and Level B: components, states, weights, adjustment,
   floor, beside the retained front; the weight version into the Stage 11
   configuration hash.
5. `api/server.py`: the §16 input contract, honouring and refusing by name.
6. `verify_ranking_final.py`: W1–W11, and N4 re-specified.
7. The candidate package and the run report carry the scorecard.
8. Full suite, both indications, every count reported beside the count predicted.

**The first run reports its scoring arithmetic before any biology**: applicable
0.96, measured 0.60, fraction 0.625 are predicted in §6, and if the run lands
elsewhere that is the first thing said, not a footnote after the results.
