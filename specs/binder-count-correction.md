# Implementation spec: counting target-matched binders

**Status: criteria fixed before implementation. No code changed at the time of
writing.** The decision this implements was taken in
`decision-binder-count-objective.md` and approved, with one addition: both the
original and the corrected evidence must survive in the output, so that the
audit trail reads as a chain rather than as a corrected answer with its
history removed.

## The rule

One clause, as approved:

> **Count a binder unless its recorded antigen names something else.**

Formally: a retrieved binder contributes to the objective unless its
target-match verdict is `FAIL`. `PASS` counts. `UNKNOWN` counts.

`UNKNOWN` counting is deliberate and is the part most likely to be misread, so
it is stated plainly here. `UNKNOWN` arises three ways and the rule admits all
three:

| how UNKNOWN arises | why it still counts |
| --- | --- |
| sequence route: a named therapeutic recorded against the target, carrying no deposited antigen annotation | the check that catches a wrong-antigen structural hit does not apply. Its UNKNOWN is about an inapplicable check, not about doubt |
| structure route: the entry records no antigen at all | absent annotation is not evidence of a different antigen. Treating absence as disqualifying would be imputation in the unfavourable direction |
| no reference record for the target, so no name set could be derived | the failure is on our side, not the binder's |

The third row is the most generous reading and the one to watch. **The
implementation must report how many binders count through each path**, so a
reader can see whether a corrected count rests on annotation that matched or on
annotation that was absent. If most of a counted set is UNKNOWN-for-want-of-
annotation, the correction is weaker than it looks, and that has to be visible
rather than left to be inferred.

This is not a new rule. The benchmark's own binder ranking already partitions on
exactly `target_match != FAIL`. The change makes the pipeline consistent with a
rule the platform already applies elsewhere.

## What each candidate carries

The requirement is that the audit trail reads as a chain: original evidence,
validation finding, corrected evidence, ranking change. A record carrying only
the corrected number makes the first half of that chain unreproducible from the
platform.

Each candidate therefore carries all four links.

| field | link in the chain |
| --- | --- |
| `retrieved_binder_count` | original evidence — every row the search returned |
| `wrong_antigen_binder_count` | validation finding — rows annotated against another protein |
| `binder_count` | corrected evidence — the objective, and what the front is computed from |
| `on_front_under_retrieved_count`, `decision_under_retrieved_count` | ranking change — the front and the decision as they would have been under the old objective |

Both fronts are computed on every run. The old objective is not read back from a
stored answer; it is recomputed from the same inputs beside the corrected one,
so the comparison cannot drift away from the result it sits next to.

`binder_count` keeps its name. It is one of the four Pareto objectives, and
renaming it would break the output contract for a change of meaning that its two
new siblings already make explicit. What the name means is pinned by an emitted
constant, `BINDER_COUNT_BASIS`, which also enters the configuration hash.

## What moves

- **The Pareto front.** `binder_count` is one of four objectives and the front is
  computed from objective values.
- **The decision column**, for any candidate whose front membership moves.
- **The Stage 11 configuration hash** — but only because this spec makes it move.
  See W13. It does not move on its own, which is a correction to the decision
  document and is recorded below.
- **Prose in four code locations and one document** which state that the
  objective counts database hits.

## What does not move

Each of these is a claim, so each gets a criterion rather than a sentence.

- **The gates.** `no binder retrieved` still reads retrieval, not matching. A
  design holding one wrong-antigen binder still passes that gate and is still
  ranked; it simply stops scoring a point for it. Changing the gate as well
  would move GPR35 to `REQUIRES_EVIDENCE` rather than `BACKUP` — a larger claim
  than the one approved, and not made here.
- **Attrition**, which follows from the gates being untouched.
- **The Level B score.** `binder_quality` reads the binder record, not the count,
  and returns `UNKNOWN` whether or not a binder was retrieved, because no
  affinity measurement is connected either way. The component's *state* is
  identical under both objectives, so applicable weight, measured weight, scored
  fraction and `overall` are all unchanged.
- **Affinity.** Stays `UNKNOWN`. Fold error and correlation stay reported as not
  computable. No estimated value is introduced anywhere by this change.
- **Attractiveness, safety margin and cleanliness.**

## A correction to the approved decision document

The decision document states that the Stage 11 configuration hash changes, "so
no cached run compares equal to a run under the old objective."

**As the code stands, it would not change.** The hash payload contains the
Stage 9 hash, the gene list, the gate strings, the recommendation set, the
decision tokens, the gate-status tokens and the scoring hash. Nothing in it
describes how an objective is computed. Verified rather than inspected:

```
Any term describing how binder_count is COMPUTED? False
hash: 22831ea0509a5c59
```

Changing the objective would therefore have left the hash identical and a cached
run under the old objective would have compared equal to a new one — the exact
silent-equality failure the chain exists to prevent. The assertion was made in
the decision document and was not true of the code.

`BINDER_COUNT_BASIS` is added to the hash payload to make the assertion true.
W13 tests it by computing both hashes rather than by reading the payload.

## Criteria

Fixed before the run. Positive pins wherever a criterion could otherwise clear
by finding nothing.

| id | criterion |
| --- | --- |
| **W12** | Every survivor satisfies `retrieved == binder_count + wrong_antigen`. The three counts are arithmetically closed, so no binder is dropped or double-counted. **Positive pin:** at least one candidate must carry `wrong_antigen >= 1`, or the criterion has proved nothing and trips. |
| **W13** | The Stage 11 configuration hash under `BINDER_COUNT_BASIS` differs from the hash under any other basis string, asserted by computing both. |
| **W14** | Both fronts are recomputed on the same run, and at least one candidate's `decision` differs from its `decision_under_retrieved_count`. **Positive pin:** if no decision differs, the change did not reach the ranking and this trips rather than clears. |
| **W15** | The attrition counts and every candidate's `gate_status` are identical to the run under the old objective. The gate did not move. |
| **W16** | Every survivor's `overall`, `scored_fraction`, `applicable_weight` and `measured_weight` are identical to the run under the old objective, to within 1e-12. Level B did not move. |
| **W17** | No candidate anywhere in the output carries a numeric affinity, fold error or correlation. All three read as not computable, each with a reason. |
| **W18** | The counted set is broken down by path — matched, sequence-route, no-annotation, no-record — and the four sum to `binder_count`. **Positive pin:** the breakdown is emitted for every survivor with `binder_count >= 1`, so a design that counts a binder always says which path counted it. |

## Predicted outcome, fixed before the run

- The front becomes `['FER1L6']`.
- **GPR35 moves from `ADVANCE` to `BACKUP`.** Its `retrieved_binder_count` stays
  1, its `wrong_antigen_binder_count` is 1, its `binder_count` becomes 0.
- FER1L6, TMEM92, TNFSF9 and BTNL8 retrieved no binder, so all three counts are
  0 for each and no decision changes.
- Attrition is unchanged at every gate.
- Every Level B `overall` is unchanged.
- One design advances where two did.

**If the run lands anywhere else, that is reported before anything else, and
this spec is amended and the run repeated rather than the result explained.**
