# Stage 4 — target pairing and the single-versus-dual decision

Written before `stages/stage4.py` exists and before any pairing output exists.
Every pool definition, threshold and rejection criterion below is fixed by this
document. If a criterion trips, the correction is a change to this document
followed by a re-run — never a narrative explaining why the output was
acceptable after all.

Status: **awaiting review. Nothing in `stages/stage4.py` may be written until
this is approved.**

Section references of the form "Stage 3 §6" point at
`specs/stage3-target-discovery.md`.

---

## 0. A precondition that is not currently met

Stage 3 at revision `c21794c` **trips R13**. Measured on a re-run of
`verify_ranking.py` at that revision:

```
TRIPPED  R13: clearance rate PROTEIN_CONFIRMED 0.6% vs RNA_SUPPORTED 42.8%,
              ratio 75.32x against a limit of 5x
cleared at 0.15   646 of 3,480
12/13 criteria clear
```

R13 was added by the calibration commit itself, to police exactly the defect the
calibration was meant to remove. It still trips: the gate still selects for
absence of evidence, at a ratio of 75x against a limit of 5x. The calibration
moved the staining axis (High is now worth 0.4626, not 1.0) but not far enough to
make the two axes commensurable at a 0.15 ceiling.

The clearance table at HEAD, by evidence class:

| ceiling | DATA_INSUFFICIENT | PROTEIN_CONFIRMED | RNA_SUPPORTED | total |
| --- | --- | --- | --- | --- |
| 0.15 conservative | 0 / 62 (0.0%) | 11 / 1,935 (0.6%) | 635 / 1,483 (42.8%) | 646 |
| 0.35 moderate | 0 / 62 (0.0%) | 156 / 1,935 (8.1%) | 969 / 1,483 (65.3%) | 1,125 |
| 0.60 permissive | 0 / 62 (0.0%) | 1,353 / 1,935 (69.9%) | 1,325 / 1,483 (89.3%) | 2,678 |

And the fact that decides §3. **Not one of the five known targets for this
indication clears at the conservative ceiling** — nor do five further known CAR
targets checked alongside them:

| gene | tier rank | composite | risk | risk organ | cleared |
| --- | --- | --- | --- | --- | --- |
| CEACAM5 | 1 | 0.8769 | 0.600 | gi_tract | no |
| PSCA | 8 | 0.7233 | 0.600 | gi_tract | no |
| CEACAM6 | 11 | 0.7000 | 0.814 | lung | no |
| MSLN | 16 | 0.6652 | 0.637 | lung | no |
| CLDN18 | 22 | 0.6415 | 0.727 | lung | no |
| EPCAM | 37 | 0.6163 | 0.753 | kidney | no |
| TACSTD2 | 46 | 0.6045 | 0.841 | kidney | no |
| MET | 72 | 0.5470 | 0.463 | lung | no |
| MUC1 | 87 | 0.5177 | 0.922 | liver | no |
| ERBB2 | 106 | 0.4984 | 0.657 | kidney | no |

All ten rank inside the top decile of their tier — R1 clears — and all ten are
blocked by the gate. The 646 cleared targets contain none of them.

By this project's own rule — *a tripped criterion gets a spec change and a
re-run, never a narrative explanation* — **Stage 3 is not in a state to be built
on.** Stage 4 is specified here, but it must not be run until R13 clears, because
the input it consumes is a clearance flag that R13 says is not yet meaning what
it claims to.

**This also means the pool premise needs restating.** The figure of **154**
cleared targets does not appear anywhere in this repository. HEAD produces
**646** at the conservative ceiling, and there is no stash, no dangling object
and no unpushed branch holding a run that produces 154. Either 154 comes from
work that was lost, or from a Stage 3 whose R13 fix has not been written yet.

The part of the premise that **is** confirmed at HEAD is the part that matters:
MSLN, CLDN18 and CEACAM6 do not clear at conservative. So do CEACAM5, MUC1 and
every other known target checked. The pool argument in §3 rests on that, not on
the size of the cleared set.

§3 is therefore written so that **it does not depend on which number is right.**
That is not a workaround; it is the correct design, and the reasoning in §3 is
the same whether the cleared set has 154 members or 646.

---

## 1. Scope

Stage 3 produced, per surface protein, three numbers that are never combined: a
tumour-side composite, a normal tissue risk with a pass/fail against the project
ceiling, and an evidence confidence. Stage 4 takes that output and answers two
questions:

1. For a pair of antigens under **AND-gating** — the therapy engages a cell only
   when it carries both — what is the pair's combined normal tissue risk, and
   what fraction of malignant cells does the gate still reach?
2. For each individual target, is the right design **single** antigen or **dual**?

`target_antigen` is still `None` and is still not written. Stage 4 produces a
recommendation and the evidence behind it; it does not seed anything.

Stage 4 introduces exactly **two new measurements**. Everything else is inherited
from Stage 3 unchanged:

| new measurement | what it is |
| --- | --- |
| conjunction risk | normal tissue risk for a pair, with the conjunction taken per organ and per cell type *before* any maximum |
| conjunction coverage | fraction of malignant cells positive for both antigens, measured per cell |

**Stage 4 does not produce a second composite score.** Reusing the name
"composite" under different semantics would create two numbers that look alike
and mean different things, and the discipline of three separate numbers would not
survive it. Where a pair needs a tumour-side ordering, `min(composite_A,
composite_B)` is used and is labelled a **bound**, not a score: under AND-gating
the weaker antigen limits the pair, and the bound is reported as such.

## 2. Sources

No new external source. One new derived artefact.

| source | supplies | status |
| --- | --- | --- |
| Stage 3 output | composite, risk, per-organ scores, confidence, clearance, evidence class | inherited |
| Human Protein Atlas (v23) | per (organ, cell type) staining call | re-read one aggregation step earlier |
| GTEx (v10) | per-tissue bulk median transcript | inherited |
| GEO GSE202051 | **per-cell** expression on malignant cells | **new artefact** |

The criticality table, the tier weights, the pancreas override and its rationale,
the exact-match label tables for both tissue vocabularies, the cultured-cell-line
exclusions and the atlas calibration curve are all **inherited from Stage 3
unchanged**. Stage 4 writes no new tissue mapping.

That is deliberate. The three tissue-mapping bugs in Stage 3 §6 all came from
mapping code, and the surest way not to reintroduce them is not to write a second
mapping. Any change to those tables is a Stage 3 change, made there, and it moves
the Stage 3 configuration hash.

---

## 3. The pool — what Stage 4 draws from, and why

This is the first thing to settle, because getting it wrong makes the stage
incapable of its own purpose in a way that would not show up as an error.

### 3.1 The rule

**Pool `P` = every surface protein from Stage 3 satisfying both:**

- **(i) it has a composite score** — measured weight `W >= MINIMUM_MEASURED_WEIGHT
  = 0.40`, so it is above the evidence floor; and
- **(ii) its normal tissue risk is defined** — at least one organ resolves, from
  staining or from transcript.

**Neither condition mentions the risk gate.** `cleared` is carried as a label on
every member and reported on every row. It is never a filter.

Both sides of a pair draw from this same pool. There is no separate partner pool.

Measured at HEAD, so the size of what is being specified is not a surprise:

| quantity | value |
| --- | --- |
| universe | 3,480 |
| excluded by (i), no composite | 66 |
| excluded by (ii), risk undefined | 66 |
| **pool `\|P\|`** | **3,413** |
| pool by class | PROTEIN_CONFIRMED 1,935 · RNA_SUPPORTED 1,478 · DATA_INSUFFICIENT 0 |
| `\|P ∩ cleared\|` | 645 |
| pairs evaluated | 5,822,578 |

The two exclusion conditions overlap but are not identical: 66 each, 67 in union.
All 62 DATA_INSUFFICIENT proteins fall out, together with 5 others — the run must
report which condition removed each, since the two mean different things. That
class emptying itself is the evidence floor and the undefined-risk rule doing
their job at the pool boundary, not a third rule added here; DATA_INSUFFICIENT is
defined by neither source resolving, which is close to the negation of both
conditions.

5.8 million pairs is tractable and is not reduced. The per-organ conjunction is
chunked by anchor — one anchor against all of `P` is a 3,413-row elementwise
minimum over the organ x cell type matrix — and the co-expression conjunction is
a popcount over packed bits (§6.2).

### 3.2 Why the cleared set is not the pool

MSLN, CLDN18 and CEACAM6 do not clear at the conservative ceiling. They are three
of the five known targets for this indication, and a dual-antigen architecture is
the thing that exists to make targets like them usable.

At HEAD it is stronger than that. **None of the five clears, and none of five
further known CAR targets clears either** (§0). Every one of the ten ranks inside
the top decile of its tier on the tumour side and every one is blocked by the
gate. The cleared set is not merely missing three interesting proteins; it is
disjoint from the entire set of antigens anyone has built a CAR against.

The Stage 3 risk gate is a **single-antigen** gate. It asks one question: if a CAR
recognises this one antigen, does any critical organ carry it? A target that fails
only because one organ carries it is precisely the case an AND gate is built to
rescue — the gate's claim is that the *second* antigen is absent from that organ,
so the first antigen being present there is not, on its own, disqualifying.

A Stage 4 restricted to the cleared set could not reach MSLN, CLDN18 or CEACAM6
at all. It would then report that dual-antigen design is unnecessary for this
indication — from a pool constructed so that no other answer was reachable. That
is not a finding. It is the shape of the pool being read back out, and it would
be indistinguishable in the output from a real negative result.

The cleared set is not the pool. It is **one column of the output**, and the most
direct evidence that this stage did anything is a count that can only be produced
by carrying it as a label: *how many clearing pairs consist of two members,
neither of which cleared alone.*

This reasoning is why §0's discrepancy does not change anything here. Whether the
cleared set has 154 members or 646, it is the wrong object to draw from, for the
same reason.

### 3.3 Why the partner side is not filtered either

A partner's job is to be **absent where the anchor is present**. A partner that is
itself high-risk, in a *different* organ, is a good partner — that is the entire
mechanism. Restricting partners to low-risk proteins would discard exactly those
whose risk is concentrated somewhere the anchor is clean, which is the population
the architecture exists to exploit.

### 3.4 Why condition (ii) is not optional

This is the one place something gate-like does happen, and it follows Stage 3 §6
directly: undefined risk is not low risk. For a pair it is worse than for a single
target.

The pair's entire safety claim is *"these two antigens are not on the same normal
cell."* An unmeasured member cannot support that claim in any organ. Admitting
unmeasured proteins as partners would let the AND gate manufacture safety out of
ignorance — the same failure Stage 3's evidence floor was built to stop, one level
up, and more dangerous here because the architecture makes it look like a design
feature rather than a gap.

Condition (i) keeps the anchor side honest for the symmetric reason: a target with
no composite has no tumour-side case, and pairing does not create one.

### 3.5 What the run must report about the pool

- `|P|`, and `|P|` by evidence class
- counts excluded by (i) and by (ii) separately, and their union
- `|P ∩ cleared|`, reported as a label — and a large departure from the Stage 3
  clearance count is a sign the conditions were applied differently there. At
  HEAD the gap is exactly one protein: 646 cleared, 645 of them in the pool.
- pairs evaluated, `|P|(|P|-1)/2`
- pool membership for MSLN, CLDN18, CEACAM6, CEACAM5 and MUC1, individually

No pair is excluded from **computation**. Reporting is capped (§9), and the cap
and the number of omitted rows are printed.

---

## 4. What makes a pair better than either target alone

Under AND-gating a pair is better than either member in exactly **one** currency,
and worse or equal in every other. Stating this first, because it decides how the
rest of the stage is scored.

**Better in:**

- **Normal tissue risk.** A normal cell is engaged only if it carries both
  antigens. This is the entire point of the architecture and the only thing the
  pair buys.

**Worse or equal in:**

- **Coverage.** `f_AB <= min(f_A, f_B)`, always. The gate reaches fewer malignant
  cells than either antigen alone. Strictly worse, never better.
- **Escape resistance.** Escape requires losing **one** antigen, not both.
  AND-gating strictly *reduces* escape resistance relative to a single target.
  This is the opposite of OR-gating and tandem constructs, and conflating the two
  is the standard error in this area. The single-positive malignant cells are not
  projected future escape variants — they are pre-existing ones, present at
  diagnosis, and their fraction is measured directly (§6.4, `unaddressed`).
- **Manufacturing.** Two binders and gated signalling against the Stage 1 budget
  of 3.5 kb and two genetic edits. Strictly more expensive.
- **Tumour-side attractiveness.** Bounded by the weaker member.

### The consequence

**A pair must be justified by clearance, not by a smaller number.** A pair that
lowers risk from 1.00 to 0.72 against a ceiling of 0.15 has bought nothing and has
paid for it in coverage, escape resistance and construct budget. The ceiling is a
gate, not a score; crossing it is the only risk improvement that changes what can
be built.

The pairing is doing work in exactly one situation:

```
pair_cleared(A, B)  and  not cleared(A)      # the anchor is rescued
```

and its strongest form, where neither member cleared alone. Everything else
Stage 4 emits is a report, not a recommendation.

---

## 5. Combined normal tissue risk under AND-gating

### 5.1 The conjunction is taken before any maximum

Stage 3 §6 computes, for a single target `T`:

```
atlas_score_T(o)      = max over cell types c in o of calibrated_level(T, o, c)
expression_score_T(o) = max( atlas_score_T(o), baseline_score_T(o) )
risk(T)               = max over organs o of [ expression_score_T(o) x criticality(o) ]
```

For a pair the conjunction must be taken at **the finest resolution each source
offers, before any maximum is applied.**

Taking `min(risk_A, risk_B)` is the wrong answer, and it is wrong in a way that is
easy to miss because the number it returns looks reasonable:

> A is present only in lung. B is present only in liver. No organ carries both, so
> the true pair risk is zero and the design is safe. `min(risk_A, risk_B)` returns
> 1.0.

`min(risk_A, risk_B)` is not conservative-but-usable. It is **blind**: it returns
the same answer whether or not the two antigens are ever co-located, which is the
only question being asked. It is also exactly what a Stage 4 that does nothing
would produce, which is why it is the null hypothesis in P1.

### 5.2 The formulas

**Atlas arm, resolved to cell type:**

```
pair_atlas(o) = max over cell types c in o of
                    min( calibrated_level(A, o, c), calibrated_level(B, o, c) )
```

The `min` inside because a CAR meets a cell, and one cell type carrying both
antigens is a real target however the organ averages out. The `max` outside for
the reason Stage 3 §6 already measures and defends. The only change from Stage 3
is that the conjunction is taken one step earlier.

**Baseline arm, bulk, organ resolution only:**

```
pair_baseline(o) = min( baseline_score_A(o), baseline_score_B(o) )
```

The transcript baseline is a bulk median per tissue with no cell type axis, so no
conjunction is observable within it. `min` here is the Fréchet upper bound on the
co-expressing fraction: it cannot be exceeded, and it is attained only when one
antigen's positive cells nest entirely inside the other's. Recorded as a **bound**,
not as a measurement.

**Combining the arms and the organs:**

```
pair_expression(o) = max( pair_atlas(o), pair_baseline(o) )
pair_risk          = max over organs o of [ pair_expression(o) x criticality(o) ]
```

`max` between the arms, not `min`, for the same reason as Stage 3 §6: a
non-detection in one source must never cancel a positive reading in the other.
That understates risk, which is the dangerous direction.

### 5.3 The bulk arm can erase the benefit, and the run must say how often

The bulk arm is a looser bound than the cell-type-resolved arm. Wherever
`pair_baseline(o) >= pair_atlas(o)` at the organ that sets the pair's risk, the
cell type resolution bought nothing and the safety case rests on a bound that
cannot see co-expression at all.

**The run must report the number and share of pairs whose risk-setting organ is
decided by the baseline arm.** If that is most of them, the AND gate's safety case
is an upper bound wearing the clothes of a measurement, and the output must say so
rather than print the number unqualified.

The complementary count — pairs where cell type resolution strictly lowers risk
below the organ-level bound — is what says the atlas arm is earning its place.
Both are reported whether or not any criterion trips.

### 5.4 Unresolved organs — where this stage is most likely to go wrong

Stage 3 tolerates per-organ gaps. A target's risk is the maximum over the organs
that resolve, and it fails the gate only when **no** organ resolves anywhere.

That tolerance does not survive conjunction. If organ `o` is measured for A and
unmeasured for B, the `min` over `o` is undefined. Treating it as zero, or skipping
the organ, credits the pair with an absence nobody observed. Since the pair's whole
claim *is* an absence, this is the failure mode that matters most in this document.

**Rule.** For organ `o`, an arm is unresolved for the pair when either member has
no measurement on that arm. The organ is unresolved for the pair only when it is
unresolved on **both** arms — one resolved arm is enough to bound it.

**Two risks are computed for every pair:**

- `pair_risk_optimistic` — unresolved organs contribute nothing.
- `pair_risk_conservative` — each unresolved organ contributes
  `measured_member_score(o) x criticality(o)`; the missing member is assumed
  present wherever nobody looked. Where **neither** member is measured in `o`, the
  organ contributes `criticality(o)` in full.

**Clearance uses the conservative value:**

```
pair_cleared = pair_risk_conservative <= normal_tissue_risk_ceiling
```

Where `pair_risk_optimistic <= ceiling < pair_risk_conservative`, the pair is
flagged `RISK_UNRESOLVED` and **does not clear**. It is reported with the explicit
list of organs that would have to be measured to settle it.

This is deliberately strict, and the strictness is what makes it useful. A pair in
`RISK_UNRESOLVED` is not a rejection: it is a named, short, specific list of
experiments that would decide a specific design. That is the most actionable thing
this stage can emit, and it exists only because the two risks are computed
separately instead of one of them being quietly chosen.

### 5.5 The conservative pair risk can exceed both members' individual risks

This looks wrong and is not. `pair_risk_conservative` fills organs that Stage 3
never had to fill, because Stage 3 took a maximum over the organs it had and a gap
simply did not contribute. Conjunction cannot do that, so the gaps become visible
and can push the conservative value above `risk_A` and above `risk_B`.

Recorded here so it is not later "fixed". `pair_risk_optimistic` **is** bounded by
both members — invariant I5, §10.

### 5.6 Inherited, not re-derived

The calibration curve mapping staining levels onto the transcript axis is computed
in Stage 3 at run time and is part of the experiment, not of the code. Stage 4
**consumes the same curve** and does not recompute it. The Stage 4 configuration
hash includes the Stage 3 configuration hash verbatim; a Stage 4 result is not
interpretable without knowing which Stage 3 produced it.

The cell-type-level staining calls Stage 4 reads are the same rows Stage 3
aggregates. Stage 4 reads them one aggregation step earlier. It does not read a
different table.

---

## 6. Co-expression on malignant cells

### 6.1 Why marginals will not do

The cached cell atlas artefact holds group means and compartment means. Those give
the mean expression of A and of B over malignant cells. They cannot answer what
fraction of malignant cells carry **both**.

Estimating `f_AB` as `f_A x f_B` assumes independence, and the departure from
independence is precisely what the AND gate's coverage depends on. Two antigens
marking disjoint malignant subpopulations have `f_AB ≈ 0` and marginals identical
to two antigens marking the same cells.

**A pairing stage built on marginals ranks pairs by the product of each target's
individual score and therefore adds nothing to Stage 3.** That is rejection
criterion P4, and it is one of the three criteria written specifically to detect
this stage doing nothing.

### 6.2 The new artefact

Stage 4 requires a per-cell presence matrix over the pool. Built by streaming the
cell atlas matrix in row blocks — it is CSR, 224,988 x 22,164, and the existing
loader already streams at 8,192 rows — retaining only malignant cells and only
genes in `P`. Cached under the same manifest discipline as every other artefact,
with the pool identity in the fingerprint, so a changed pool invalidates it rather
than silently reusing the wrong columns.

Measured, not assumed: the malignant compartment holds **64,538 cells** across
**43 patients**.

Retained in two forms:

- the sparse submatrix in the atlas's own units, so the positivity threshold can be
  perturbed (P12) without re-streaming 8.3 GB;
- a packed bit matrix at the fixed threshold. 64,538 bits is 8.1 kB per gene, so a
  3,413-gene pool is about 28 MB and every pairwise conjunction is a popcount over
  an AND.

The dense float form is not retained: at 3,413 genes it would be 880 MB, and the
sparse form carries the same information because most entries are zero.

### 6.3 Per-cell positivity

```
positive(cell, gene)  <=>  expm1(X[cell, gene]) >= THETA,   THETA = 1.0
```

`X` is `log1p(CP10K)` and per-cell `expm1` sums to roughly 9,600, so `THETA = 1.0`
transcripts per 10k corresponds to approximately **one captured molecule**. That is
the correct floor for a dropout-limited nuclear assay: anything lower counts noise,
anything meaningfully higher discards real single-molecule detections and compounds
the dropout problem this atlas already has.

`THETA` is a free parameter and is perturbed as one (P12), for the reason Stage 3's
R12 gives: it looks like implementation detail and is actually policy.

### 6.4 The quantities

Computed on malignant cells, pooled and per patient:

| quantity | definition | reads as |
| --- | --- | --- |
| `f_A`, `f_B` | marginal positive fractions | reach of each antigen alone |
| `f_AB` | double-positive fraction | **coverage** — what the gate reaches |
| `retention` | `f_AB / max(f_A, f_B)` | what the gate gives up against the better single antigen |
| `dependence` | `f_AB / (f_A x f_B)` | lift over independence; `>1` co-occurring, `<1` mutually exclusive |
| `unaddressed` | `(f_A + f_B - 2 f_AB) / (f_A + f_B - f_AB)` | share of antigen-positive malignant cells the gate does not reach — pre-existing escape, measured rather than projected |

### 6.5 The threshold that matters is `retention`, not `coverage`

This is the substantive threshold question, and the answer is driven by a defect
Stage 3 already measured.

The assay deflates `f_A`, `f_B` and `f_AB` **together**, so a ratio between them
survives the deflation far better than any absolute fraction does. Stage 3 measured
the magnitude directly: CEACAM5 reads 0.000065 in this atlas while sitting at 299
transcripts and 409x normal in bulk. An absolute coverage floor would reject every
CEACAM5 pair on a capture artefact — the exact error Stage 3's dropout rule exists
to prevent, committed one stage later.

**Thresholds, fixed here:**

| threshold | value | role |
| --- | --- | --- |
| `retention` | `>= 0.50` | **primary.** The gate must reach at least half the malignant cells the better single antigen reached. |
| `f_AB` pooled | `>= 0.10` | secondary, deliberately low. Present only to exclude pairs of two rare antigens that score well on retention while addressing almost nothing. Explicitly deflated by dropout; not to be read as a coverage estimate. |
| per-patient | `f_AB >= 0.10` in `>= 60%` of evaluable patients | a pair double-positive in half the patients and absent in the rest is pooled-identical to one that is 50% double-positive in every patient. Those are different products. |

The per-patient floor is the cell-level analogue of C4 patient prevalence, and it
exists because pooling 43 patients hides exactly the failure a bulk cohort would
also hide.

### 6.6 Evaluable patients

All 43 patients carry at least one malignant cell, but the per-patient counts run
7,167 / 5,790 / 4,865 / ... / 9 / 6 / 4 / 3. A proportion cannot be estimated from
three cells.

**A patient is evaluable when it contributes at least 100 malignant cells.**
Measured: **29 of 43** qualify (31 at a floor of 50, 25 at 200). At 100 cells the
standard error on a proportion near 0.3 is about 0.046, tight enough for a 0.10
threshold; at 30 cells it is 0.084, which is not.

The 14 excluded patients are reported by identifier and cell count. Not silently
dropped.

### 6.7 Subsets

Computed on `all` and on `untreated` separately, as Stage 2 requires of this source.

- `all`: 43 patients, 64,538 malignant cells, 29 evaluable
- `untreated`: 18 patients, 52,999 malignant cells

**Scored on `all`.** A conjunction is the quantity most damaged by small `n`, and
the untreated subset holds 18 patients against 43 — it retains 82% of the cells but
loses more than half the patients, and the per-patient floor is a patient-count
test. The treated samples are also a real clinical population, not a contaminant.

`untreated` is reported beside it on every row, and `subsets_disagree` is raised
where `retention` differs by more than 0.20 between them. Same pattern as C3's two
denominators and for the same reason: two comparators that are not interchangeable,
neither trustworthy alone. The flag exists to be read.

### 6.8 The dropout rule, carried forward and made stricter

Stage 3 §4.1: a value at or below `DROPOUT_EPSILON = 0.001` never rejects a target
and never scores it as zero.

Stage 3 counted 358 proteins silent at that threshold **across every cell type
group**. Stage 4 needs a different and larger count, because it reads the
malignant compartment specifically. Measured at HEAD over the pool:

| | count |
| --- | --- |
| pool members with no cell atlas row at all | 119 |
| pool members at or below the epsilon on the malignant compartment mean | 675 |
| pool members with a usable malignant mean | 2,619 |
| **pairs measurable for co-expression** | **3,428,271 of 5,822,578 (58.9%)** |

**Two fifths of all pairs cannot be measured for co-expression at all.** That is
the single largest limitation on this stage and it is a property of the assay, not
of the pairs. It is reported in the header, not discovered in the tail of a table.

At pair level:

- If **either** member's malignant compartment mean is at or below the dropout
  epsilon, or either member has no row in the cell atlas at all, its per-cell
  positive fraction is not a measurement. The pair carries
  `CO_EXPRESSION_NOT_MEASURED` and receives no `retention`, `f_AB`, `dependence` or
  `unaddressed` value. **Not scored zero. Not dropped from the output.**
  The two causes are recorded separately — `dropout_silent` and `no_atlas_row` —
  because one is an assay failure on a protein that is present and the other is an
  absence of any observation, and only the first names an experiment worth doing.
- A pair carrying that flag **cannot be recommended** (P6), because the coverage
  claim underpinning a recommendation was never made.

It can still be reported as a pair whose safety case holds and whose coverage is
unknown. For CEACAM5 — the top-ranked target in Stage 3 — that is the honest
statement, and it is also the actionable one: it names a specific experiment, which
is to measure the two antigens on the same section by flow or immunohistochemistry,
in an assay that does not have this atlas's capture problem.

The run must report how many pool members are dropout-silent and how many pairs are
consequently unmeasurable.

---

## 7. The single-versus-dual decision

Taken **per target `T`**, over `T`'s admissible partners.

### 7.1 Admissibility

A partner `Q` is admissible for `T` when **all** hold:

1. `pair_cleared(T, Q)` — conservative pair risk at or below the project ceiling
2. `retention >= 0.50`, `f_AB >= 0.10` pooled, and `f_AB >= 0.10` in at least 60%
   of evaluable patients
3. not `CO_EXPRESSION_NOT_MEASURED`
4. `Q` is in `P`, so it has its own tumour-side case. The gate requires `Q` on the
   cell being killed; an antigen with no tumour-side evidence is not a partner, it
   is a liability.
5. `Q != T`

### 7.2 The four outcomes

Every target in `P` receives exactly one.

- **`SINGLE`** — `cleared(T)`. Recorded together with the best admissible partner
  if one exists, and the explicit note that dual was available and not taken. The
  reasoning is §4: the pair is strictly worse in coverage, escape resistance and
  construct budget, and the ceiling is a gate rather than a score. Nothing is
  bought by crossing it twice.
- **`DUAL`** — `not cleared(T)` and at least one admissible partner exists. The
  recommended partner is the admissible partner with the highest `retention`, ties
  broken by lower `pair_risk_conservative`, then by higher `composite(Q)`.
  **Ranked on retention, not on risk**: among admissible partners the risk question
  is already settled — all of them clear — and what separates them is how much of
  the tumour the gate still reaches.
- **`NO_DESIGN`** — `not cleared(T)` and no admissible partner. Reported with a
  breakdown of which admissibility condition the best candidates failed, counted. A
  stage that says "no" without saying which wall it hit cannot be acted on.
- **`UNRESOLVED`** — `not cleared(T)`, no admissible partner, but at least one
  partner would be admissible under `pair_risk_optimistic` and fails only on
  `RISK_UNRESOLVED`. Emitted with the specific organs that would have to be
  measured. This class converts a limitation into an experiment, and it is reported
  separately from `NO_DESIGN` precisely so it is not read as one.

The distribution across the four outcomes is reported. P11 trips if it is
degenerate.

### 7.3 What Stage 4 does not decide

- **Construct feasibility.** Stage 1 fixes the budget at 3.5 kb after the 1.2 kb
  backbone overhead, and `max_genetic_edits` at 2, and admits `dual_target` and
  `logic_gated` among the allowed formats. Stage 4 does **not** size a construct —
  that is not this stage's measurement — but a `DUAL` recommendation is meaningless
  if the format is disallowed, so the recommendation carries the allowed-format
  list from the resolved spec, and the stage fails loudly if both `dual_target` and
  `logic_gated` are absent. Sizing belongs to a later stage and is named here so it
  is not mistaken for having been done.
- **Binder availability.** That is Stage 5. A partner with no retrievable binder is
  **flagged, not filtered**. Excluding it here would import a Stage 5 result into a
  Stage 4 decision, and the pool would then be shaped by which proteins happen to
  have been crystallised — a property of the literature, not of the biology.
- **Which gating chemistry.** AND-gating is the architecture this document assumes.
  Whether it is realised by split signalling, a synthetic receptor cascade or
  otherwise is a later decision and changes no number here.

---

## 8. Evidence confidence for a pair

The third number, reported alongside the risk and coverage quantities and **never
combined with either**.

`pair_confidence` reflects: both members' Stage 3 confidence, the number of organs
resolved on both arms for both members, whether co-expression was measurable at
all, and the number of evaluable patients behind the per-patient test.

It is bounded by its weaker member — `pair_confidence <= min(confidence_A,
confidence_B)` — because a pair cannot be better evidenced than the least evidenced
antigen in it. Asserted as invariant I4, not assumed.

The reason it stays separate is the reason Stage 3 §7 gives, and pairing sharpens
it: a pair of two well-stained proteins measured and found **not** to co-occur, and
a pair of two proteins nobody has looked at, can produce the same
`pair_risk_optimistic`. One number cannot carry both, and the second of those is
exactly what §5.4 exists to catch.

---

## 9. Output and reproducibility

The header carries everything Stage 3's header carries, plus:

- **the Stage 3 configuration hash, verbatim**, and the Stage 3 R-criteria outcome
  — a Stage 4 run on top of a tripped Stage 3 must say so on its own first page
- the pool definition and `|P|`, by evidence class, with the cleared count as a
  label and the exclusion counts by condition
- `THETA`, the retention floor, the pooled coverage floor, the per-patient floor,
  the evaluable-patient minimum, and the resulting evaluable patient count
- the number of pairs evaluated
- pool members that are dropout-silent, and pairs consequently unmeasurable
- the share of pairs whose risk-setting organ is decided by the bulk arm (§5.3)
- the excluded patients by identifier and cell count
- Stage 4's own configuration hash, covering all of the above **including the
  Stage 3 hash**

The Stage 4 hash must be verified stable across processes, for the reason Stage 3
§9 gives.

**No silent caps.** The reported table is capped at the top 50 partners per anchor.
The cap and the number of omitted rows are printed, and the full table is written to
file. A bounded report that does not say it is bounded reads as complete coverage
when it is not.

---

## 10. Rejection criteria — fixed in advance

Prefixed `P` to keep them distinct from Stage 3's `R1`–`R13`. Stage 3's criteria
apply to the Stage 3 run that feeds this one and are not re-checked here.

A tripped criterion means this document changes and the run repeats. It never means
the result gets an explanation.

### Construction invariants — assert, do not report

A failure here is a bug, not a result. The run stops.

| id | invariant |
| --- | --- |
| I1 | `pair_risk(T, T) == risk(T)` for every `T` in `P` |
| I2 | `pair_risk(A, B) == pair_risk(B, A)` |
| I3 | `f_AB <= min(f_A, f_B)` for every measured pair |
| I4 | `pair_confidence <= min(confidence_A, confidence_B)` |
| I5 | `pair_risk_optimistic <= min(risk_A, risk_B)` |

I1 is the important one and is checked across the whole pool, not spot-checked.
Pairing a target with itself must reproduce Stage 3's risk exactly, since
`min(x, x) = x` at every level of the aggregation. If it does not, the conjunction
machinery and the single-antigen machinery disagree about what an organ score is,
and one of them is wrong.

I5 holds for the optimistic value only. See §5.5.

### Criteria

| id | criterion |
| --- | --- |
| P1 | `pair_risk_conservative` correlates above 0.95 (Spearman) with `min(risk_A, risk_B)` across all evaluated pairs |
| P2 | fewer than 1% of pairs achieve `pair_risk_optimistic < min(risk_A, risk_B) - 0.05` |
| P3 | no target that failed the single-antigen gate is rescued by any admissible partner |
| P4 | `f_AB` correlates above 0.98 (Spearman) with `f_A x f_B` across pairs where both members are measured |
| P5 | a pair clears whose clearance depends on an organ unresolved for either member |
| P6 | a pair carrying `CO_EXPRESSION_NOT_MEASURED` is recommended |
| P7 | a pair containing a ubiquitous immune protein (HLA-A/B, CD74, PTPRC) clears |
| P8 | more than 10% of clearing pairs clear because one member is unmeasured in the organ that sets the other member's risk |
| P9 | a recommended pair's coverage is concentrated in fewer than 60% of evaluable patients |
| P10 | `DUAL` is recommended for any target that already clears alone |
| P11 | the single-versus-dual decision returns the same outcome for more than 95% of the pool |
| P12 | a 2x change in `THETA`, the retention floor or the pooled coverage floor changes more than half of the `DUAL` recommendations |
| P13 | the same protein is the recommended partner for more than half of all `DUAL` targets |

### P1, P2 and P4 — the criteria that would tell us the pairing logic is doing nothing

Three criteria exist specifically to detect this stage adding nothing to Stage 3.
All three report their measured values whether or not they trip.

**P1 and P2 — the null result on the risk side.** `min(risk_A, risk_B)` is the
answer a stage produces if it never resolves organs: take each target's whole-organ
maximum and pick the smaller. It requires no conjunction, no cell type axis and no
second source. If `pair_risk` is rank-equivalent to it, then no pair is ever safe
*because the two antigens sit in different organs* — only because one of them was
safer to begin with — and Stage 4 has re-derived Stage 3's ordering at quadratic
cost. P1 tests the ordering; P2 tests whether **any** pair strictly beats its better
member by a margin large enough to matter. A stage can pass one and fail the other,
and failing either means the same thing.

**P4 — the null result on the coverage side.** `f_A x f_B` is what the marginals
alone give, and the marginals are already on disk in the cached compartment means.
If the measured double-positive fraction is rank-equivalent to the product, the
per-cell pass over the matrix produced no information that was not already
available, and the co-expression measurement is decorative.

The correlation is computed **only over pairs where both members are above the
dropout epsilon**. The rest have no measurement to correlate, and zero-filling them
inside the check that polices the co-expression measurement would be the same
imputation these documents forbid, committed by the auditor. Same rule and same
reason as Stage 3's R5 detail.

### P3 detail — deliberately universal

P3 is written as *"no target that failed the single-antigen gate is rescued"*, not
as *"MSLN is not rescued"*.

MSLN, CLDN18 and CEACAM6 are the three targets this architecture most obviously
exists for, and they are on a **reported watch list**: their outcome — pool
membership, best partner, pair risk both ways, retention, final decision — is
printed on every run whether or not they are rescued. But their individual rescue
is **not** a criterion. Requiring a named protein to survive is how a screen gets
tuned into a confirmation of what was already believed, and Stage 3 refused that
for the same reason.

What P3 asserts is weaker and sound: if the AND-gate machinery cannot rescue **any**
of the targets the single-antigen gate blocked, then either the conjunction is not
working or the architecture does not help this indication, and both of those need
finding out rather than reporting.

### P8 detail — the Stage 3 disease, one level up

Stage 3's open problem is a gate that selects for absence of evidence, and §0 says
it is still live. The pair version is worse, because it is disguised as a design
feature: a pair "clears" because member B was never measured in the organ where
member A is dangerous, and the output reads as a successful AND gate.

§5.4's conservative arm is what prevents it; P8 is what verifies §5.4 is actually
doing that. If P5 is correct by construction, P8 should be near zero. If P8 is
large while P5 passes, the conservative arm is being computed but not being used.

### Explicitly not grounds for rejection

- MSLN, CLDN18 or CEACAM6 specifically not being rescued
- every target resolving to `SINGLE`, if the risk gate genuinely clears them
- the best partner being a protein nobody has proposed as a CAR target
- a pair scoring lower on the tumour side than either member — that is arithmetic,
  not a defect
- a pair with high `dependence` but low absolute coverage, where both members sit
  near the dropout threshold

These are possible correct answers.

---

## 11. Expected results

**There is no reference run for this stage.** Nothing below is a number to hit.
This section lists what must be reported, so the absence of an expected value is
not mistaken for the absence of an obligation to report.

Required in the output regardless of what any of them turn out to be:

- `|P|`, the pool composition, and the pool's overlap with the cleared set
- pairs evaluated; pairs clearing; pairs clearing where **neither** member cleared
  alone
- the four-way outcome distribution over `P`
- P1, P2 and P4's measured correlations
- the share of pairs whose risk is set by the bulk arm (§5.3), and the share where
  cell type resolution strictly lowers risk
- dropout-silent pool members and unmeasurable pairs
- the watch list: MSLN, CLDN18, CEACAM6, CEACAM5, MUC1
- `RISK_UNRESOLVED` pairs with the organs that would settle them

**If a count comes out surprising, that is a result.** If a count comes out
impossible — a negative coverage, `f_AB > min(f_A, f_B)`, a self-pair risk that does
not match Stage 3 — that is a bug and the run stops.

---

## 12. Known open problems — carried forward, not fixed here

**1. Stage 3's R13 is still tripped.** See §0. Stage 4 must not be run for
interpretation until it clears. The specification is written now because the pool
question (§3) is independent of it, and because writing the spec first is the rule.

**2. The cell-type-count imbalance bites harder here than in Stage 3.**

Stage 3 §6 recorded that organs differ greatly in how many cell types the atlas
records: 41 for the gastrointestinal tract, 35 for brain, against 1 for heart and 2
for liver, and left it uncorrected.

Under conjunction that residual stops being a reporting caveat and becomes a
structural asymmetry. The atlas arm's benefit comes entirely from taking the `min`
*within* a cell type before the `max` across cell types. An organ with one recorded
cell type has no within-organ structure to exploit, so the pair can never be shown
to separate inside it. **Heart and liver therefore get almost no benefit from the
architecture, while the gastrointestinal tract gets the most — and heart and liver
are tier 1.** The safety benefit of AND-gating is unevenly distributed across
exactly the organs where it matters most.

Not fixable with this data. Printed with every result.

**3. Cell-type co-presence is still not cell co-presence.**

The atlas records a call per cell type, not a joint distribution over cells within
that cell type. Two antigens called present in the same cell type are not
necessarily on the same cell. Even the resolved arm is an upper bound.

This produces an asymmetry worth stating plainly: **Stage 4 measures co-expression
per cell where it strengthens the tumour case, and bounds it per cell type where it
would strengthen the safety case** — because a per-cell joint measurement exists for
the malignant compartment and does not exist for normal tissue. The asymmetry runs
in the conservative direction, which is why it is acceptable, and it is stated
rather than hidden because the two numbers look alike in the output and are not the
same kind of thing.

**4. Dropout bounds the entire coverage side.**

794 of the 3,413 pool members — 675 silent on the malignant compartment plus 119
with no atlas row — carry no usable per-cell measurement. That removes **41% of
all pairs** from consideration for co-expression, including every pair containing
CEACAM5, the top-ranked target in Stage 3. This is a property of a nuclear assay,
not of the biology, and the correct response is a different assay rather than a
lower threshold.

**5. The highest-value data addition for this architecture is a normal-tissue
single-cell atlas.**

With one, the safety arm could be measured per cell instead of bounded per cell
type, problems 2 and 3 would both dissolve, and `pair_risk` would become a
measurement rather than an upper bound. Named here so it is on the record as the
specific gap, rather than being rediscovered later as a surprise.

---

## Build note

Implementation order once approved, and **after R13 clears**:

1. the per-cell artefact in `data/singlecell.py` — malignant-cell submatrix over the
   pool, cached, with the pool in the fingerprint
2. `stages/stage4.py` — conjunction risk, then coverage, then the decision
3. `verify_pairing.py` — the five invariants first, then the thirteen criteria, and
   only then any reading of the biology

Same order and same rule as Stage 3: the invariants and criteria run before anyone
looks at which pairs came out on top.
