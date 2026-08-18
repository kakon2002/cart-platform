# Stage 4 — target pairing and the single-versus-dual decision

Written before `stages/stage4.py` exists and before any pairing output exists.
Every pool definition, threshold and rejection criterion below is fixed by this
document. If a criterion trips, the correction is a change to this document
followed by a re-run — never a narrative explaining why the output was acceptable
after all.

Status: **awaiting review. Nothing in `stages/stage4.py` may be written until this
is approved.**

Section references of the form "Stage 3 §6" point at
`specs/stage3-target-discovery.md`.

---

## 0. Preconditions and two numbers that do not reconcile

### 0.1 Stage 3 trips R13

Measured on a re-run of `verify_ranking.py` at `c21794c`:

```
TRIPPED  R13: clearance rate PROTEIN_CONFIRMED 0.6% vs RNA_SUPPORTED 42.8%,
              ratio 75.32x against a limit of 5x
cleared at 0.15   646 of 3,480
12/13 criteria clear
```

R13 was added by the calibration commit itself, to police the defect that
calibration was meant to remove. It still trips. The staining axis moved (High is
now worth 0.4626, not 1.0) but not far enough to be commensurable with the
transcript axis at a 0.15 ceiling.

By this project's rule — *a tripped criterion gets a spec change and a re-run,
never a narrative explanation* — **Stage 4 may be specified but must not be
interpreted until R13 clears.** The pool rule in §3 does not depend on the risk
gate at all, so the specification is not blocked by this; the reading of its
output is.

### 0.2 Two supplied figures do not match the repository

| supplied | measured at HEAD | status |
| --- | --- | --- |
| 154 cleared targets | **646** | no stash, no dangling object, no unpushed branch produces 154 |
| MSLN risk 0.487 | **0.6366**, peak organ lung | same provenance question |

Both are recorded rather than reconciled. Neither changes any decision in this
document: §3 discards the cleared set as a pool, and §9's rescue table reports
each member's risk as measured on the run that produced it rather than quoting a
constant. Where MSLN appears below as a worked example it carries **0.637**.

If 154 and 0.487 came from a Stage 3 with the R13 fix already applied, that work
is not in this repository and its absence should be treated as the more urgent
problem.

---

## 1. Scope

Stage 3 produced, per surface protein, three numbers that are never combined: a
tumour-side composite, a normal tissue risk with a pass/fail against the project
ceiling, and an evidence confidence. Stage 4 answers two questions:

1. For a pair of antigens under **AND-gating** — the therapy engages a cell only
   when it carries both — what is the pair's combined normal tissue risk, and what
   fraction of malignant cells does the gate actually kill?
2. For each target, is the right design **single** antigen or **dual**?

`target_antigen` stays `None`. Stage 4 produces a recommendation and the evidence
behind it; it does not seed anything.

Two new measurements. Everything else is inherited from Stage 3 unchanged:

| new measurement | what it is |
| --- | --- |
| combined risk | normal tissue risk for a pair, minimum per organ, maximum across organs |
| co-expression | fraction of malignant cells carrying both antigens, measured per cell |

**Stage 4 produces no second composite.** Where a pair needs a tumour-side
ordering, `min(composite_A, composite_B)` is used and labelled a **bound**: under
AND-gating the weaker antigen limits the pair.

## 2. Sources

No new external source. One new derived artefact.

| source | supplies | status |
| --- | --- | --- |
| Stage 3 output | composite, risk, per-organ scores, confidence, clearance | inherited |
| Human Protein Atlas (v23) | per-organ staining score | inherited |
| GTEx (v10) | per-tissue bulk median transcript | inherited |
| GEO GSE202051 | **per-cell counts** on malignant cells | **new artefact (§6.2)** |

The criticality table, tier weights, the pancreas override and its rationale, the
exact-match label tables for both tissue vocabularies, the cultured-cell-line
exclusions and the calibration curve are **inherited unchanged**. Stage 4 writes
no new tissue mapping: the three tissue-mapping bugs in Stage 3 §6 all came from
mapping code, and the surest way not to reintroduce them is not to write a second
mapping.

---

## 3. The pool

### 3.1 The rule

**Pool `P` = the top 200 surface proteins by Stage 3 composite, with risk ignored
entirely.** 200 x 199 / 2 = **19,900 pairs**.

Risk is not a filter, not a tiebreak, and not consulted at pool construction. It
is evaluated against the ceiling afterwards, on the pair.

### 3.2 Why risk is ignored and attractiveness is the filter

**Risk is the thing pairing exists to fix, so filtering on it first defeats the
stage.** The Stage 3 gate is a *single-antigen* gate: it asks whether any critical
organ carries this one antigen. A target blocked because one organ carries it is
exactly the case an AND gate rescues, since the gate's claim is that the *second*
antigen is absent from that organ. Drawing the pool from cleared targets would
return "no dual design is needed" from a pool built so that no other answer could
emerge — a result indistinguishable in the output from a real negative.

**Attractiveness is the right filter because a pair cannot rescue a target nobody
wants.** A pair is at best as attractive as its weaker member (§4), so a pool
ordered by composite is a pool ordered by the ceiling on what any pair containing
that member can achieve.

### 3.3 What the top 200 actually contains

Measured at HEAD, so the shape of the pool is known before the stage runs:

| quantity | value |
| --- | --- |
| scored universe | 3,414 |
| composite range in the pool | 0.4474 .. 0.8769 |
| by evidence class | PROTEIN_CONFIRMED 168 · RNA_SUPPORTED 32 · DATA_INSUFFICIENT 0 |
| **cleared at 0.15** | **1 of 200** |
| risk undefined | 0 of 200 |
| risk min / median / max | 0.0277 / 0.5930 / 1.0000 |
| pairs | 19,900 |

**One member of the pool clears the single-antigen gate. The other 199 are
blocked.** That is the strongest available confirmation that the pool is pointed
at the right population: it is almost entirely made of targets that only a dual
design could use, which is what Stage 4 exists to test.

**Risk is defined for all 200**, so the undefined-risk problem does not arise
inside this pool and needs no pool-level rule. §5.4 still handles per-organ gaps,
which do occur.

All ten known CAR targets checked in §0 fall inside the top 200. Nine further
known targets do not — FOLH1 at 1,221, ROR1 1,439, GPC3 893, TNFRSF17 1,375,
L1CAM 787, CD70 1,198, NCAM1 1,536, ALPP 1,156 — on tumour-side evidence in this
indication, which is the intended behaviour rather than a miss.

**The cut at 200 is arbitrary and IL13RA2 sits at 201.** One known target is
excluded by a single position. The run reports the 20 proteins immediately below
the cut so that a reader can see what the boundary discarded, and P15 perturbs the
pool size.

### 3.4 Reported per member

`cleared`, `risk`, `composite`, `evidence_class` and `tier_rank` travel as labels
on every pool member and every pair row. The most direct evidence that this stage
did anything is a count that requires them: *how many clearing pairs consist of
two members, neither of which cleared alone.*

---

## 4. What makes a pair better than either target alone

### 4.1 The claim, stated so it can fail

**Claim.** There exists a pair `(A, B)` such that

```
combined_risk(A,B) <= ceiling                          (1)  safety is won
risk_A > ceiling                                       (2)  and was not already there
combined_risk(A,B) < min(risk_A, risk_B) - 0.05        (3)  by conjunction, not by selection
f_AB >= COVERAGE_FLOOR                                 (4)  and the gate still kills enough
f_AB >= COVERAGE_FLOOR in >= 60% of evaluable patients (5)  in most patients, not on average
```

Every line is a number the stage emits. (3) is the one that separates a result
from a restatement: a pair that is safe only because one member was already safer
is the better of its two members wearing a different name.

**What refutes it.** Each is a possible output, and each falsifies a different
part:

| observation | what it refutes |
| --- | --- |
| no pair satisfies (1) and (2) together | AND-gating does not rescue anything here (P3) |
| pairs satisfy (1) but never (3) | the conjunction is inert; risk is selection (P1, P2) |
| pairs satisfy (3) but never (4) | the conjunction works and costs more tumour than it buys |
| `f_AB` tracks `f_A x f_B` | co-expression carries nothing beyond the marginals (P4) |
| the best pair is the top two singles | the ordering is the single ordering (P14) |

### 4.2 Better in one currency, worse in the rest

- **Better: normal tissue risk.** A normal cell is engaged only if it carries
  both. This is the only thing the pair buys.
- **Worse: coverage.** `f_AB <= min(f_A, f_B)` always. Strictly fewer malignant
  cells than either antigen alone.
- **Worse: escape resistance.** Escape requires losing **one** antigen, not both.
  AND-gating strictly *reduces* escape resistance — the opposite of OR-gating and
  tandem constructs, and conflating the two is the standard error here. The
  single-positive cells are pre-existing escape variants present at diagnosis, and
  `1 - f_AB` measures the reservoir directly.
- **Worse: manufacturing.** Two binders and gated signalling against 3.5 kb and
  two edits.
- **Bounded: tumour-side attractiveness.** By the weaker member.

**A pair is therefore justified by clearance, not by a smaller number.** Lowering
risk from 1.00 to 0.72 against a ceiling of 0.15 buys nothing and costs coverage,
escape resistance and construct budget. The ceiling is a gate, not a score.

### 4.3 The costs are columns, not concessions

| cost | quantity | where |
| --- | --- | --- |
| tumour not addressed | `1 - f_AB` | §6.4 |
| **coverage of A given up by adding B** | `sacrificed_A = 1 - P(B\|A) = 1 - f_AB / f_A` | §6.4 |
| **coverage of B given up by adding A** | `sacrificed_B = 1 - P(A\|B) = 1 - f_AB / f_B` | §6.4 |
| patients below the floor | share of evaluable patients | §6.7 |
| tumour-side ceiling | `min(composite_A, composite_B)`, a bound | §1 |

Note the pairing of subscripts, which is easy to invert: what A gives up is
governed by `P(B|A)`, the share of **A's** positive cells that also carry B. A
conditional named for A in front measures B's loss, not A's.

A pair that wins on (1)–(5) and leaves `1 - f_AB = 0.95` is a safe design that
kills a twentieth of the tumour. That has to be visible in the same row.

---

## 5. Combined normal tissue risk

### 5.1 The formula

```
combined_risk(organ) = min( score_A(organ), score_B(organ) ) x criticality(organ)
combined_risk        = max over organs of combined_risk(organ)
```

Same shape as the single-target gate — per organ, then maximum across organs — so
the two numbers are directly comparable and the self-pair identity holds (I1).

**`min`, because an AND-gate only fires where both antigens are present**, so the
organ's risk is bounded by whichever antigen is scarcer there. Taking `max` would
ignore the architecture entirely and reduce the pair to its more dangerous member.

`score_X(organ)` is Stage 3's `expression_score`, unchanged: the maximum of the
calibrated staining score and the transcript baseline score for that organ, with a
non-detection in one source never cancelling a positive reading in the other.

### 5.2 This is a bound, not a measurement — and the range is reported

**There is no per-cell normal tissue data.** The transcript baseline is bulk per
tissue. The atlas is per-cell-type staining, which is not a joint distribution
over cells — two antigens called present in the same cell type need not be on the
same cell. Neither source can say what fraction of cells in an organ carry both.

So `min` is not an estimate of co-expression. It is the **Fréchet upper bound**:
the largest co-expressing fraction consistent with the two marginals, attained
only when one antigen's positive cells nest entirely inside the other's.

**Naming the bounds precisely, because it is easy to get backwards.** `min`
assumes **perfect overlap** — maximal co-expression, not minimal. Saying it
assumes "worst-case co-expression" invites the opposite reading and should be
avoided. Perfect overlap is:

- **pessimistic for safety** — it puts the largest possible number of
  double-positive normal cells in the organ, so risk comes out as high as the
  marginals permit
- **optimistic for coverage** — the same assumption applied to a tumour would put
  the largest possible number of double-positive malignant cells there

**Those two coincide**, which is what makes `min` the right thing to gate on: the
single assumption that is least favourable on safety is also the one that flatters
coverage most, so a pair that clears under `min` clears under an assumption that
was working against it in the only direction that matters.

**The other bound, reported beside it:**

```
independence_risk(organ) = score_A(organ) x score_B(organ) x criticality(organ)
independence_risk        = max over organs
```

which is what the organ's risk would be if the two antigens were distributed
independently within it.

**Gate on `min`.** Report `independence_risk`, the gap between them, and the organ
that sets each.

| number | assumption | direction | used for |
| --- | --- | --- | --- |
| `combined_risk` (`min`) | perfect overlap — maximal co-expression | highest risk the marginals allow; pessimistic for safety | **the gate** |
| `independence_risk` | independence within the organ | lower risk; optimistic for safety | reported |
| `min(risk_A, risk_B)` | no conjunction at all | what a stage that never conjoined would produce | P1, P2 null |

The truth for any real organ sits between the first two unless the antigens are
exactly nested or exactly independent. Reporting both is what lets a reader see
how much of the safety case is architecture and how much is assumption.

### 5.3 A tighter bound is available and is reported, not gated on

The atlas records staining per (organ, cell type), and Stage 3 aggregates to organ
by maximum. Taking the minimum *within* a cell type before that maximum —

```
cellwise(organ) = max over cell types c of min( level_A(organ,c), level_B(organ,c) )
```

— is strictly tighter than organ-level `min` and still an upper bound, because a
cell type is still not a cell.

It is **reported as a third column and is not the gate**, since §5.2's objection
applies to it too. Its value is diagnostic: the number of pairs where it falls
below the organ-level bound says how much structure the cell type axis is
carrying. If that count is near zero the column should be dropped at the next
revision.

Say the word if you would rather this come out entirely and the stage carry two
risk numbers instead of three.

### 5.4 Unresolved organs

Stage 3 tolerates per-organ gaps: risk is the maximum over organs that resolve.
Conjunction cannot do that. If organ `o` is measured for A and not for B, the
`min` is undefined, and treating it as zero credits the pair with an absence
nobody observed — and the pair's whole claim *is* an absence.

**Rule.** Organ `o` is unresolved for the pair when either member has no
measurement there from either source. Two risks are computed:

- `combined_risk_optimistic` — unresolved organs contribute nothing
- `combined_risk` (**the gate**) — each unresolved organ contributes
  `measured_member_score(o) x criticality(o)`, the missing member assumed present
  wherever nobody looked

**Organs measured for neither member are outside the analysis entirely** and
contribute to neither number. This is not an oversight and it is forced: Stage 3
already tolerates organs nobody measured, taking its maximum over the organs that
resolve, and I1 requires `combined_risk(T,T) == risk(T)` exactly. Filling
never-measured organs with `criticality(o)` would break that identity for every
target and would introduce ignorance the single-antigen gate never charged for.
The pair adds no new tolerance; it only refuses to *benefit* from a gap where one
member is measured and the other is not.

Where `combined_risk_optimistic <= ceiling < combined_risk`, the pair is
`RISK_UNRESOLVED` and does not clear. It is reported with the specific organs that
would settle it — which is the most actionable output this stage produces, and it
exists only because the two are computed separately rather than one being quietly
chosen.

Note that `combined_risk` can exceed both members' individual risks, because
Stage 3 never had to fill those gaps. That is the gap becoming visible, not a bug.
`combined_risk_optimistic` **is** bounded by both (I5).

---

## 6. Co-expression on malignant cells

### 6.1 The primary metric is the absolute intersection

**Primary: `f_AB`, the fraction of malignant cells expressing both antigens.**
That is what an AND-gate kills. `1 - f_AB` is the escape population.

**Jaccard is not the primary and must not be used as one.** Two antigens each on
5% of malignant cells with perfect overlap score Jaccard 1.0, and an AND-gate on
them kills 5% of the tumour. That is a number that looks right for the wrong
reason: Jaccard measures agreement between two sets and is blind to how large they
are, and the gate's value is entirely a matter of how large they are.

Reported on every pair. Only the first is primary:

| quantity | definition | what it says |
| --- | --- | --- |
| **`f_AB`** | fraction of malignant cells positive for both | **what the gate kills** |
| `P(A\|B)` | `f_AB / f_B` | how much of B's coverage survives adding A |
| `P(B\|A)` | `f_AB / f_A` | how much of A's coverage survives adding B |
| `sacrificed_A` | `1 - P(B\|A)` | share of A's own coverage given up by adding B |
| `sacrificed_B` | `1 - P(A\|B)` | share of B's own coverage given up by adding A |
| `Jaccard` | `f_AB / f_(A or B)` | overlap shape, reported for comparability only |

`sacrificed_A` and `sacrificed_B` are the conditionals restated as losses, and
they are carried on the rescue row (§9.1) because that is where the trade has to
be legible. Watch the subscripts: what **A** gives up is set by `P(B|A)`.

`P(A|B)` and `P(B|A)` are what a designer reads to decide which antigen is paying:
an asymmetric pair, where one conditional is near 1 and the other near 0.1,
sacrifices almost none of one antigen's reach and almost all of the other's.

### 6.2 The cached atlas cannot answer this and must be re-derived

Settled here rather than discovered during implementation.

**It cannot.** `SingleCellSource.build_group_means()` accumulates per-group sums
and counts inside its streaming loop; the cell axis is consumed and discarded
before anything is written. `group_means.npz` is 5.7 MB and holds:

| array | shape | has a cell axis? |
| --- | --- | --- |
| `group_means` | 78 x 22,164 | no |
| `compartment_means` | 6 x 22,164 | no |
| `per_cell_total` | 224,988 | yes, but no gene axis |

No array carries both a cell axis and a gene axis, so `f_AB` is not recoverable
from it by any computation. Two consequences that are easy to miss: there is **no
patient axis either**, so §6.7's per-patient requirement is equally unanswerable;
and a compartment *mean* is not a positive *fraction*, so even `f_A` and `f_B` are
absent.

**What is required.** A second builder over the same `h5ad`, additive:

- one further pass in 8,192-row blocks over **`layers/counts`** — raw integer
  counts, which the file carries alongside the normalised `X` (§6.3)
- keeping rows where `Level 1 Annotation == "Epithelial (malignant)"` (64,538 of
  224,988) and columns for genes in `P` (200 of 22,164)
- carrying `pid` alongside, so the per-patient split needs no third pass
- written as a **new cache entry** with its own manifest, declared in
  `cache_entries()` beside the existing one

**`group_means.npz` is not modified and its fingerprint does not move**, so
Stage 3 is not invalidated — which matters because Stage 3 must be re-run for R13
and those two re-runs must not entangle.

At 200 genes the artefact is small: a 64,538 x 200 integer matrix, plus a packed
bit matrix per threshold at 8.1 kB per gene, so all 19,900 conjunctions are
popcounts over ANDs.

### 6.3 Two implementation traps, both found by probing the file

**Trap 1: `X` is not counts.** `X` is `log1p(CP10K)`. A fixed threshold on it is
not a fixed count threshold, because it is normalised per cell and the cells
differ in depth by two orders of magnitude — raw depth over malignant cells runs
**min 96, median 2,033, max 9,642**. One count in a median cell is about 4.9
CP10K; in the shallowest it is 104. Use **`layers/counts`**, which the file
carries and which is integer.

**Trap 2: the column indices are not sorted within a row.** This file stores CSR
`indices` in *descending* order per row — the first row begins `22151, 22127,
22107, 22104, ...`. Any lookup using `searchsorted`, or any code assuming the
ascending order that CSR conventionally has, **returns zero silently for every
gene**. This was not hypothetical: the probe written for this document hit it, and
reported every watch-list gene at exactly 0.0000 while row extraction and per-cell
depth were correct. Use a full-width lookup table indexed by column, never a
sorted search.

Trap 2 is the reason §6.6's sanity check is mandatory rather than advisory. It is
a silent-wrong-answer defect of exactly the class this project keeps finding, and
the only thing that catches it is a known-answer check.

### 6.4 Detection threshold

**A cell is positive for a gene when it carries at least 1 count.**

`DETECTION_COUNTS = 1`, with the run reporting sensitivity at **2** and **3**
counts and the ranking stability across the three.

Measured on the watch list, so the sensitivity is known to be real rather than
nominal:

| gene | >= 1 | >= 2 | >= 3 | total molecules |
| --- | --- | --- | --- | --- |
| KRT19 | 0.7294 | 0.5366 | 0.4063 | 230,536 |
| EPCAM | 0.1931 | 0.0456 | 0.0124 | 16,639 |
| CLDN18 | 0.1848 | 0.0781 | 0.0407 | 24,555 |
| CEACAM6 | 0.1439 | 0.0487 | 0.0235 | 18,146 |
| PSCA | 0.1393 | 0.0741 | 0.0489 | 27,273 |
| MSLN | 0.1197 | 0.0242 | 0.0070 | 10,001 |
| MUC1 | 0.0520 | 0.0047 | 0.0005 | 3,702 |
| CEACAM5 | 0.0000 | 0.0000 | 0.0000 | **2** |

**The ordering is not stable across thresholds and the run must say so.** At >= 1
the order is EPCAM > CLDN18 > CEACAM6 > PSCA > MSLN; at >= 3 it is
PSCA > CLDN18 > CEACAM6 > EPCAM > MSLN. EPCAM falls from first to fourth and PSCA
rises from fourth to first. MSLN loses 94% of its detected cells between 1 and 3.
Any conclusion that survives only at one threshold is a conclusion about the
threshold.

### 6.5 The dropout bias is stated, never corrected

**No dropout correction is applied.** No imputation, no scaling factor, no
inflation of `f_AB` toward an estimated truth.

The bias, stated so it travels with every number:

> This is a single-nucleus assay. A cell reading zero for a gene may still express
> the protein. Measured co-expression therefore **underestimates** true
> co-expression, and the underestimate is **worse for lowly-expressed genes**,
> because a gene at low abundance is more likely to be missed in any given cell.
> The bias is systematic and it favours abundant antigens: a pair of two abundant
> antigens will out-measure a pair of two scarce ones even if the scarce pair is
> biologically better co-located.

A correction factor we cannot validate would convert a stated, bounded limitation
into a hidden, unbounded one. It belongs in the output, not in the arithmetic.

The direction is worth being precise about: since the bias understates `f_AB`, it
understates coverage and **overstates** the escape population `1 - f_AB`. On the
coverage side this is the conservative direction. It is not conservative for the
single-versus-dual decision, because it makes every pair look worse than it is and
so biases the stage toward `SINGLE` and toward abundant antigens.

### 6.6 Mandatory sanity check on the derivation

Run before any pairing. **MSLN, CLDN18 and CEACAM6 must each be detected in a
substantial fraction of malignant cells.** Measured expectation at 1 count:

| gene | expected | tolerance |
| --- | --- | --- |
| MSLN | 0.1197 | must exceed 0.05 |
| CLDN18 | 0.1848 | must exceed 0.05 |
| CEACAM6 | 0.1439 | must exceed 0.05 |
| KRT19 (positive control) | 0.7294 | must exceed 0.50 |

KRT19 is included because it is the canonical malignant marker at a compartment
mean of 14.1, so it is the loudest possible signal — a derivation broken the way
§6.3's Trap 2 breaks one returns 0.0000 for it, which is unmistakable.

**KRT19 is requested as an extra column, not drawn from the pool.** It is a
cytokeratin, so it never passes the surface filter and can never be a pool
member. That is what makes it a good control rather than a problem: it is the
loudest signal in the atlas and it is independent of everything the ranking
stage decided. The per-cell artefact is therefore built over the pool **plus**
the control genes, and the controls take no part in any pairing.

**CEACAM5 is excluded from this check and must not be used for it.** CEACAM5 reads
**2 molecules across all 64,538 malignant cells**, a detection fraction of
0.0000. This is the dropout failure Stage 3 already documented and quantified
(0.000065 compartment mean, against 299 transcripts and 409x normal in bulk), not
a derivation error. Including CEACAM5 in a check whose failure means "the
derivation is wrong" would halt a correct implementation on a known assay artefact
— and would do so on the single highest-ranked target in the pool, which is the
worst possible place for a false alarm.

The distinction the check must preserve: **MSLN and CLDN18 failing means the
derivation is broken. CEACAM5 failing means the assay dropped it.** Only the
former is a stop.

### 6.7 Stratify by patient as well as compartment

**No downsampling is performed.** All 64,538 malignant cells are used, so there is
no sample to stratify and no sampling bias to introduce. If a sample is ever
introduced it must be stratified by **patient and compartment jointly**, and that
requirement is recorded here so it is not lost.

The concern it comes from is real and is present in the full data anyway, because
patients contribute wildly unequally — 7,167 malignant cells from the largest, 3
from the smallest. Pooling lets one tumour dominate, and co-expression measured
mostly within one tumour is a statement about that tumour.

**A patient is evaluable at >= 100 malignant cells: 29 of 43** (31 at 50, 25 at
200). At 100 cells the standard error on a proportion near 0.15 is about 0.036; at
30 cells it is 0.065.

Measured per-patient spread at 1 count, which shows why this is not a formality:

| gene | pooled | median patient | min | max | patients >= 0.10 |
| --- | --- | --- | --- | --- | --- |
| KRT19 | 0.7294 | 0.7768 | 0.0110 | 0.9405 | 28/29 |
| EPCAM | 0.1931 | 0.2200 | 0.0000 | 0.4763 | 23/29 |
| CLDN18 | 0.1848 | 0.0532 | 0.0000 | 0.5391 | 13/29 |
| CEACAM6 | 0.1439 | 0.1538 | 0.0000 | 0.5814 | 18/29 |
| PSCA | 0.1393 | 0.0590 | 0.0000 | 0.8360 | 8/29 |
| MSLN | 0.1197 | 0.0614 | 0.0000 | 0.5607 | 14/29 |

**CLDN18's pooled 0.1848 is a median patient value of 0.0532.** The pooled figure
is more than three times the typical patient, because a few patients carry most of
the signal. PSCA runs 0.0000 to 0.8360 across patients. A pooled number alone
would have presented both as uniform properties of PDAC.

Every pair therefore reports pooled `f_AB` **and** the per-patient distribution:
median, min, max, and the count of evaluable patients at or above the floor.

### 6.8 Coverage floor, and why it is set low

`COVERAGE_FLOOR = 0.02`, applied to pooled `f_AB` and, separately, required in at
least 60% of evaluable patients.

**This is set against the measured range and it is uncomfortably low on purpose.**
Pairwise intersections among the watch list at 1 count:

| pair | `f_AB` | `P(A\|B)` | `P(B\|A)` | Jaccard |
| --- | --- | --- | --- | --- |
| CEACAM6 + EPCAM | 0.0470 | 0.2432 | 0.3265 | 0.1620 |
| CLDN18 + EPCAM | 0.0466 | 0.2412 | 0.2521 | 0.1406 |
| CLDN18 + CEACAM6 | 0.0451 | 0.3137 | 0.2442 | 0.1592 |
| MSLN + EPCAM | 0.0332 | 0.1720 | 0.2774 | 0.1188 |
| MSLN + PSCA | 0.0295 | 0.2119 | 0.2465 | 0.1286 |
| MSLN + CEACAM6 | 0.0288 | 0.2003 | 0.2407 | 0.1227 |
| MSLN + CLDN18 | 0.0237 | 0.1284 | 0.1982 | 0.0845 |
| MSLN + CEACAM5 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

**The best pair of known targets reaches 4.7% of malignant cells.** A floor of
0.10 would admit none of them. A floor of 0.05 would admit none of them. The floor
is 0.02 so that the pairs this stage exists to evaluate are evaluated rather than
eliminated by a threshold set without looking.

What this means has to be said rather than buried: **on this assay, no pair of
known PDAC targets addresses more than about 5% of malignant cells, so `1 - f_AB`
exceeds 0.95 for all of them.** Under §6.5's bias that is a floor on the truth,
not an estimate of it — but it is the measurement we have, and a stage that
reported a comfortable-looking number here would be reporting the threshold rather
than the biology.

**A plausible correct outcome of this stage is that AND-gating cannot demonstrate
adequate coverage for any pair in this atlas.** That is a finding, not a failure,
and P16 exists so that it is reported as one rather than resolved by lowering the
floor again.

### 6.9 Subsets

Computed on `all` (43 patients, 64,538 malignant cells, 29 evaluable) and
`untreated` (18 patients, 52,999 malignant cells) separately. **Scored on `all`**:
a conjunction is the quantity most damaged by small `n`, the untreated subset
holds 82% of the cells but fewer than half the patients, and the per-patient floor
is a patient-count test. `subsets_disagree` is raised where pooled `f_AB` differs
by more than a factor of two.

---

## 7. The single-versus-dual decision

### 7.1 Admissibility

Partner `Q` is admissible for target `T` when all hold:

1. `combined_risk(T,Q) <= ceiling`
2. `f_AB >= 0.02` pooled, and `f_AB >= 0.02` in at least 60% of evaluable patients
3. the pair is not `CO_EXPRESSION_NOT_MEASURED` (§7.3)
4. `Q` is in `P`
5. `Q != T`

### 7.2 The four outcomes

- **`SINGLE`** — `cleared(T)`. The pair is strictly worse in coverage, escape
  resistance and construct budget, and the ceiling is a gate rather than a score,
  so nothing is bought by crossing it twice. The best admissible partner is still
  recorded, with an explicit note that dual was available and not taken.
- **`DUAL`** — `not cleared(T)` and an admissible partner exists. The recommended
  partner is the admissible partner with the **highest `f_AB`** — among partners
  that clear, the risk question is settled and what separates them is how much
  tumour the gate still kills. Ties broken by lower `combined_risk`, then higher
  `composite(Q)`.
- **`NO_DESIGN`** — `not cleared(T)`, no admissible partner. Reported with a count
  of which condition the best candidates failed on. A stage that says no without
  saying which wall it hit cannot be acted on.
- **`UNRESOLVED`** — no admissible partner, but at least one would be admissible
  under `combined_risk_optimistic` and fails only on `RISK_UNRESOLVED`. Emitted
  with the organs that would settle it.

Given §3.3 — 199 of 200 blocked — `SINGLE` should be rare by construction, and P11
trips if the outcome distribution is degenerate.

### 7.3 Pairs whose co-expression cannot be measured

Carrying forward Stage 3 §4.1: a single-cell zero never rejects a target. If
either member has fewer than 10 detected malignant cells in total, its per-cell
positive fraction is not a measurement. The pair carries
`CO_EXPRESSION_NOT_MEASURED`, receives no `f_AB`, and is **not scored zero and not
dropped**. It cannot be recommended (P6), because the coverage claim behind a
recommendation was never made.

Every pair containing CEACAM5 takes this path — 199 of the 19,900. For the
top-ranked target in the pool the honest output is a pair whose safety case holds
and whose coverage is unknown, which names a specific experiment: measure the two
antigens on the same section in an assay without this capture problem.

### 7.4 Not decided here

- **Construct feasibility.** Stage 1 fixes 3.5 kb and two edits and admits
  `dual_target` and `logic_gated`. Stage 4 does not size a construct, but fails
  loudly if both formats are disallowed.
- **Binder availability.** Stage 5. Flagged, never filtered — filtering here would
  shape the pool by which proteins happen to have been crystallised.
- **Gating chemistry.** AND-gating is assumed; how it is realised changes no
  number here.

---

## 8. Evidence confidence for a pair

Third number, never combined with the other two. Reflects both members' Stage 3
confidence, organs resolved for both, whether co-expression was measurable, and
the number of evaluable patients. Bounded by the weaker member (I4).

A pair of two well-stained proteins measured and found not to co-occur, and a pair
of two proteins nobody has looked at, can produce the same
`combined_risk_optimistic`. One number cannot carry both.

---

## 9. Output

### 9.1 The rescue table — the headline of this stage

**For every pair that clears, report which single-target risks it rescued and by
how much**, on the pair's own row rather than as something a reader reconstructs:

| column | meaning |
| --- | --- |
| `risk_A`, `risk_B` | each member's single-antigen risk, as measured on this run |
| `cleared_A`, `cleared_B` | whether each cleared alone |
| `combined_risk` | the pair's gated risk |
| `delta_A`, `delta_B` | `risk_X - combined_risk` — how far the risk moved |
| **`rescued`** | **which members go blocked -> cleared** |
| `organ_A`, `organ_B`, `organ_pair` | the risk-setting organ for each |
| `f_A`, `f_B`, `f_AB` | coverage of each alone, and of the pair |
| **`sacrificed_A`** | `1 - f_AB / f_A` — the share of A's coverage given up |
| **`sacrificed_B`** | `1 - f_AB / f_B` — the share of B's coverage given up |

**Rescue is a condition, not a statistic.** A member is rescued only when

```
risk_X > ceiling  and  combined_risk <= ceiling
```

A large `delta` that does not cross the ceiling is **not** a rescue and must not
rank as one. A pair moving MSLN from 0.637 to 0.20 has a `delta_A` of 0.437 and an
empty `rescued` field, and it sorts below a pair with a smaller delta that lands at
0.13. `delta` is reported because it says how far a pair got; it is never sorted
on and it never substitutes for the condition.

Sorting: `rescued` non-empty first, then `f_AB` descending. Within the rescued set
the safety question is already answered for every row, so what separates them is
coverage.

**The trade is visible in the same row.** A pair that rescues MSLN from 0.637 to
under 0.15 but reaches only 30% of the malignant cells MSLN reached alone has
traded a safety problem for an efficacy problem. `sacrificed_A = 0.70` says so on
the row where the rescue is claimed, rather than leaving a reader to reconstruct it
from `f_A` and `f_AB` elsewhere. A rescue reported without its cost is half a
result.

### 9.2 Header

Everything Stage 3's header carries, plus: the Stage 3 configuration hash verbatim
and its R-criteria outcome (a Stage 4 run on a tripped Stage 3 says so on its first
page); the pool rule, size, composite range, class composition and cleared count;
the 20 proteins immediately below the cut; `DETECTION_COUNTS` and the sensitivity
table at 1/2/3; `COVERAGE_FLOOR`, the per-patient floor and the evaluable-patient
count; the §6.5 dropout bias statement in full; the §6.6 sanity check results;
excluded patients by identifier and cell count; pairs unmeasurable for
co-expression; and Stage 4's own configuration hash covering all of it including
the Stage 3 hash. Verified stable across processes.

**No silent caps.** All 19,900 pairs are computed and written to file. Any capped
display states the cap and the number omitted.

---

## 10. Rejection criteria — fixed in advance

Prefixed `P`. Stage 3's `R1`–`R13` apply to the run feeding this one and are not
re-checked here.

### Construction invariants — assert, do not report

| id | invariant |
| --- | --- |
| I1 | `combined_risk(T,T) == risk(T)` for every `T` in `P` |
| I2 | `combined_risk(A,B) == combined_risk(B,A)` |
| I3 | `f_AB <= min(f_A, f_B)` for every measured pair |
| I4 | `pair_confidence <= min(confidence_A, confidence_B)` |
| I5 | `combined_risk_optimistic <= min(risk_A, risk_B)` |
| I6 | `f_AB` at 1 count `>=` `f_AB` at 2 counts `>=` at 3 counts |
| I7 | `independence_risk <= combined_risk` — the optimistic bound never exceeds the gate |

I7 follows from scores lying in `[0, 1]`, where `score_A x score_B <= min(score_A,
score_B)` organ by organ. It is asserted rather than assumed because it is the
arithmetic that makes §5.2's two bounds a range rather than two unrelated
numbers; if it ever failed, the pair of them would be reported as an interval
that does not contain what it claims to.

**Both risk numbers are carried at four decimals**, the precision the ranking
stage stores its own risk at. I1 and I5 compare across the two stages, and a
full-precision pair risk against a rounded single risk disagrees at the fifth
decimal for arithmetic reasons rather than substantive ones. The invariants
absorb that boundary and print the worst observed gap, so a real disagreement
cannot hide underneath the tolerance.

I1 is checked across the whole pool, not spot-checked: `min(x,x) = x` at every
level, so pairing a target with itself must reproduce Stage 3 exactly. If it does
not, the two machineries disagree about what an organ score is.

### Criteria

| id | criterion |
| --- | --- |
| P1 | `combined_risk` correlates above 0.95 (Spearman) with `min(risk_A, risk_B)` across all 19,900 pairs |
| P2 | fewer than 1% of pairs achieve `combined_risk < min(risk_A, risk_B) - 0.05` |
| P3 | **no blocked target is rescued by any pair** — no target moves from blocked to cleared |
| P4 | `f_AB` correlates above 0.98 (Spearman) with `f_A x f_B` across measured pairs |
| P5 | a pair is marked cleared on the optimistic arm — `optimistic <= ceiling < combined` |
| P6 | a pair carrying `CO_EXPRESSION_NOT_MEASURED` is recommended |
| P7 | a pair containing a ubiquitous immune protein (HLA-A/B, CD74, PTPRC) clears |
| P8 | more than 10% of clearing pairs stop clearing when unresolved organs are charged at full criticality |
| P9 | a recommended pair's coverage is concentrated in fewer than 60% of evaluable patients |
| P10 | `DUAL` is recommended for a target that already clears alone |
| P11 | the decision returns the same outcome for more than 95% of the pool |
| P12 | moving `DETECTION_COUNTS` from 1 to 2 changes more than half the `DUAL` recommendations |
| P13 | the same protein is the recommended partner for more than half of all `DUAL` targets |
| P14 | **the top-ranked pair is the top two singles put together** — the best pair is `(rank 1, rank 2)` by composite |
| P15 | halving or doubling the pool size changes more than half the `DUAL` recommendations |
| P16 | no pair anywhere reaches `f_AB >= 0.02` — the floor admits nothing |

### The four that say the stage is doing nothing

**P1 and P2 — risk is selection, not conjunction.** `min(risk_A, risk_B)` is what a
stage produces if it never conjoins: take each target's whole-organ maximum and
pick the safer. If `combined_risk` is rank-equivalent to it, no pair is ever safe
*because the two antigens sit in different organs* — only because one was safer to
start with, and the pair score is the better of its two members wearing a different
name. P1 tests the ordering, P2 whether any pair strictly beats its better member
by a margin that matters.

**P14 — the ordering is the single ordering.** If the best pair is simply the two
highest-composite singles, the pairing did no selecting: it inherited Stage 3's
order and joined the top of it. This is the sharpest single-shot version of the
same objection, and unlike P1 it is checkable by eye in the first row of output.

**P4 — coverage is the marginals.** `f_A x f_B` is the cheap half of the new
artefact: 200 column sums, no pairwise work. `f_AB` is the 19,900 conjunctions. If
they are rank-equivalent, the expensive half was implied by the cheap half.
Computed only over measured pairs; zero-filling unmeasurable ones inside the check
that polices the measurement would be the imputation these documents forbid,
committed by the auditor.

**P3 — nothing moved.** If AND-gating never takes a target from blocked to
cleared, the stage is doing nothing regardless of how its scores look. Written
universally, not as "MSLN must be rescued": requiring a named protein to survive is
how a screen becomes a confirmation of what was already believed. MSLN, CLDN18,
CEACAM6, CEACAM5 and MUC1 are a **reported watch list** whose outcome is printed
either way.

### P5 and P8 detail — two different questions about ignorance

Stage 3's open problem is a gate that selects for absence of evidence. The pair
version would be worse, because it is disguised as a design feature: a pair
"clears" because member B was never measured where member A is dangerous, and the
output reads as a successful AND gate.

**P5 is the wiring check.** Clearance must be decided on `combined_risk`, never on
`combined_risk_optimistic`. A pair that clears only on the optimistic arm is
`RISK_UNRESOLVED` by definition and must not be marked cleared.

**P8 is the substantive one, and it is not the same question.** Merely having an
unresolved organ is not a failure: the conservative arm already charges the
measured member's own score there, so such a pair cleared *despite* the gap. P8
asks whether clearance would survive the unmeasured antigen **saturating** the
organ nobody looked at. Written as a number that needs a third fill:

```
pessimistic(o) = criticality(o)                            o unresolved for the pair
               = min(score_A(o), score_B(o)) x criticality(o)   otherwise
```

computed and reported, and **never the gate**. Charging full criticality at every
gap would fail almost everything and would make the pool's measurement coverage,
rather than its biology, the thing under test.

An earlier draft phrased both loosely enough that they were first implemented as
"has any unresolved organ" — a property of the data rather than of the clearance.
The wording above is deliberately operational for that reason.

### P16 detail

§6.8 sets the floor at 0.02 against a measured best of 0.047 among known targets.
If nothing in 19,900 pairs reaches it, that is the stage's finding: on this assay
AND-gating cannot demonstrate adequate coverage. P16 exists so that outcome is
reported as a result and the floor is not quietly lowered a second time to produce
survivors.

### Explicitly not grounds for rejection

- MSLN, CLDN18 or CEACAM6 specifically not being rescued
- the best partner being a protein nobody has proposed as a CAR target
- a pair scoring lower on the tumour side than either member — arithmetic, not a
  defect
- low absolute `f_AB` across the board, given §6.5
- every target resolving to `SINGLE`, if the gate genuinely clears them

---

## 11. Expected results

**No reference run exists.** Nothing here is a number to hit. Required in the
output whatever the values:

- pool composition, composite range, cleared count, the 20 below the cut
- pairs clearing; pairs clearing where **neither** member cleared alone
- the rescue table (§9.1), sorted with rescues first then by `f_AB`, carrying
  `sacrificed_A` and `sacrificed_B` on every rescued row
- the four-way outcome distribution
- P1, P2, P4 and P14's measured values
- `f_AB` distribution over all pairs, pooled and per patient
- sensitivity at 1/2/3 counts and the ranking stability across them
- the §6.6 sanity check, printed
- `combined_risk` versus `independence_risk`, and how often they straddle the ceiling
- pairs unmeasurable for co-expression

**A surprising count is a result.** An impossible one — `f_AB > min(f_A, f_B)`, a
self-pair risk that does not match Stage 3, negative coverage — is a bug and the
run stops.

---

## 12. Open problems carried forward

**1. R13 is still tripped.** §0.1. Specify now, interpret later.

**2. Dropout bounds the whole coverage side, and it is not small.** §6.5. The best
known pair sits at 4.7% and CEACAM5 has 2 molecules in 64,538 cells. Every
coverage number in this stage is a floor on the truth. The correct response is a
different assay, not a correction factor.

**3. There is no per-cell normal tissue data, so the safety side is bounded where
the tumour side is measured.** §5.2. Stage 4 measures co-expression per cell where
it helps the tumour case and bounds it per organ where it would help the safety
case, because a per-cell joint measurement exists for the malignant compartment
and does not exist for normal tissue. The asymmetry runs conservative, which is
why it is acceptable, and it is stated because the two numbers look alike in the
output and are not the same kind of thing.

**4. Organs differ in how many cell types the atlas records** — 41 for the
gastrointestinal tract against 1 for heart and 2 for liver. This bites §5.3's
tighter bound specifically: an organ with one recorded cell type has no structure
to exploit, so heart and liver gain nothing from cell type resolution while the
gut gains most. Heart and liver are tier 1.

**5. The highest-value data addition is a normal-tissue single-cell atlas.** With
one, §5.2's bound becomes a measurement and problems 3 and 4 both dissolve.

---

## Build note

Once approved, and with §0.1 understood:

1. the per-cell artefact in `data/singlecell.py` — a **new** cache entry beside
   `group_means`, over `layers/counts`, malignant cells x pool genes, `pid`
   carried, fingerprinted on the pool. `group_means.npz` untouched.
2. **the §6.6 sanity check, run and passing, before anything else is written** —
   §6.3's Trap 2 returns zeros silently and only a known-answer check catches it
3. `stages/stage4.py` — combined risk, then co-expression, then the decision
4. `verify_pairing.py` — six invariants, then sixteen criteria, then the biology
