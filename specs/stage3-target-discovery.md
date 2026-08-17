# Stage 3 — target discovery, screen and ranking

Written before the ranking is implemented and before any output exists. Every
threshold, weight and rejection criterion below is fixed by this document. If a
criterion trips, the correction is a change to this document followed by a
re-run — never a narrative explaining why the output was acceptable after all.

Status: **awaiting review. Nothing in `stages/stage3.py` may be written until
this is approved.**

---

## 1. Scope

Input is a cancer type with no antigen supplied. `target_antigen` is `None` and
stays `None` through this stage; the null is what selects discovery mode, and
Stage 3 is the thing that would eventually fill it. Seeding it to unblock
anything downstream voids the exercise.

Stage 3 produces, for every protein in the surface universe:

- a **composite attractiveness score**, or an explicit refusal to score it
- a **normal tissue risk** value and a pass/fail against the project ceiling
- an **evidence confidence** value

These are three separate numbers. They are never combined into one. Combining
attractiveness with confidence lets a protein nobody measured score like a
protein measured and found good; combining risk with confidence lets a protein
nobody measured score like a protein measured and found safe. Both failures run
in the dangerous direction.

## 2. Sources and the names used below

| registry name | referred to below as | supplies |
| --- | --- | --- |
| UniProt | the proteome | topology, anchor, ectodomain annotation |
| Human Protein Atlas | the tissue atlas | normal tissue staining, localisation |
| GTEx | the normal baseline | per-tissue median transcript level |
| GDC TCGA | the tumour cohort | per-sample tumour transcript level |
| DepMap | the dependency screens | per-line gene effect |
| GEO GSE202051 | the cell atlas | per cell type expression |

## 3. Universe

The surface set as it stands after step 5: **3,480 proteins** (reference run:
3,496). Evidence classes as built in step 11:

| class | count | reference |
| --- | --- | --- |
| PROTEIN_CONFIRMED | 1,935 | 1,945 |
| RNA_SUPPORTED | 1,483 | 1,491 |
| DATA_INSUFFICIENT | 62 | 60 |

Ranking is reported within evidence class ("tier"). Cross-tier comparison of
composites is not meaningful and is not presented as though it were.

## 4. Composite score

Six components. **Weights are fixed here, before any output exists**, on stated
reasoning. They are not to be adjusted after seeing where anything ranks.

| # | component | weight | source |
| --- | --- | --- | --- |
| C1 | malignant cell expression | 0.25 | cell atlas |
| C2 | malignant versus stroma specificity | 0.20 | cell atlas |
| C3 | tumour versus normal margin | 0.25 | tumour cohort + normal baseline |
| C4 | patient prevalence | 0.15 | tumour cohort |
| C5 | surface accessibility | 0.10 | proteome + tissue atlas |
| C6 | escape resistance | 0.05 | dependency screens |

**This is the tumour side only. Safety is not a term in it.** Normal tissue risk
is a gate applied separately (§6), not a negative component folded into the sum.
A protein cannot buy its way past an unmanageable tissue by scoring well here.

### Why the dependency screens carry the lowest weight

Measured, not assumed: only **16 of 3,480** surface proteins (reference: 17 of
3,496) reach a median gene effect at or below −0.5 in this lineage. The
component cannot discriminate for 99.5% of the universe, so it cannot carry
real weight.

It is **not** inert, and is not dropped. The reference run recorded it moving
the ordering at 3,336 positions. The re-run must report the same measurement —
the count of positions whose rank changes when C6 is removed — and that number
belongs in the output.

### Component definitions

Each component maps a measurement to `[0, 1]`, or to `NOT_MEASURED`. Saturation
points are stated so the scale is absolute rather than relative to whatever
happens to be in the universe on a given run.

- **C1 malignant cell expression.** `mean expm1` expression in the malignant
  compartment, `m`. Score `clamp(log10(1 + m) / log10(1 + 100), 0, 1)`.
  Saturates at 100 transcripts per 10k. `NOT_MEASURED` when the protein has no
  row in the cell atlas.
- **C2 malignant versus stroma specificity.** `r = m / p`, where `p` is the peak
  of the fibroblast, immune and endothelial compartment means. Peak, not mean:
  one compartment expressing it is enough to matter. Score
  `clamp(log10(r) / log10(50), 0, 1)`. `NOT_MEASURED` when the protein has no
  row. Where `p` falls below the dropout epsilon the ratio is not computed from
  it — see §4.1.
- **C3 tumour versus normal margin.** **Both fold changes are computed and both
  are reported.** Neither denominator is trustworthy alone:

  - `fold_baseline` = median tumour level ÷ the normal baseline's pancreas
    median. **This is the primary**, because the cohort's own normal arm is n=4
    and every ratio built on it rests on four samples.
  - `fold_cohort` = median tumour level ÷ the cohort's solid-normal median.
    Same assay and same pipeline as the tumour side, so it carries no
    cross-cohort batch effect — but almost no samples.

  The two are not interchangeable: they come from different pipelines with
  different normalisation, so the primary carries a batch effect the secondary
  does not, and the secondary carries a sampling error the primary does not.
  Scoring uses `fold_baseline`; `fold_cohort` is reported beside it in every row.

  Score `clamp(log2(fold_baseline) / log2(64), 0, 1)`.

  **Disagreement is flagged, not resolved.** Where the two folds differ by more
  than 2× in either direction, the row carries `sources_disagree`. That flag is
  what R11 polices — it exists to be read, not to be silently overridden by
  whichever denominator was picked first.

  `NOT_MEASURED` when the tumour side or the baseline pancreas value is absent.
  Absence of the cohort normal alone does not block scoring; it suppresses the
  comparison and is recorded as such.
- **C4 patient prevalence.** Fraction of primary tumour samples at or above 10
  transcripts per million. Already `[0, 1]`. `NOT_MEASURED` when the protein has
  no column in the cohort.
- **C5 surface accessibility.** See §5. `NOT_MEASURED` where the ectodomain is
  unannotated.
- **C6 escape resistance.** Median gene effect `e` across screened lines in the
  lineage. Score `clamp(-e / 1.0, 0, 1)`. `NOT_MEASURED` where the gene has no
  row, **and separately where every line in the lineage was unscreened** — a
  column of not-a-number is not a zero.

### 4.1 The dropout rule

The cell atlas is a nuclear assay and drops transcripts that bulk measurement
finds abundantly present. In the reference run CEACAM5 reads 0.0001 there while
sitting at 299 transcripts and 409× normal in bulk; this run measures its peak
across all cell type groups at **0.000065**.

- The silence threshold is `0.001`, not zero. An exact-zero test misses this
  case entirely. Measured sensitivity in this run: 253 proteins silent at 0.0,
  **358** at 0.001, 530 at 0.01 (reference: 267 / 357 / 534).
- **A value at or below the threshold in the cell atlas never rejects a
  target and never scores it as zero.** C1 and C2 become `NOT_MEASURED` for that
  protein and their weight redistributes. This source separates compartments; it
  does not refute.
- The same applies to the C2 denominator: where the stromal peak is at or below
  the threshold, the ratio is not computed and C2 is `NOT_MEASURED` rather than
  infinite.

CEACAM5 retains 0.55 of measured weight under this rule, so it still scores.

### 4.2 Missing components and the evidence floor

`MINIMUM_MEASURED_WEIGHT = 0.40`.

Let `W` be the summed weight of components that are measured for a protein.

- If `W >= 0.40`: `composite = (Σ wᵢ · sᵢ) / W` over measured components only.
  Missing weight redistributes proportionally across the rest.
- If `W < 0.40`: **no composite is produced.** The protein receives a
  within-tier rank, its full component breakdown, and a
  `below_evidence_floor` flag. It is listed, not hidden.

The bound is the point of the rule. Unbounded redistribution lets a single
measured component stand in for the whole score: the reference run found **61
proteins** measured on one 0.10-weight component alone rescaling to composites
of 0.850 — above CEACAM6 at 0.709, CLDN18 at 0.678 and MSLN at 0.655. They were
overwhelmingly retroviral envelopes and pseudogene entries. That is absent
evidence manufacturing a high position, and it is the exact failure the floor
exists to stop.

**No component is ever imputed.** Not to zero, not to a midpoint, not to the
universe mean. Missing propagates to the evidence confidence number (§7) and to
the floor test, never into the composite or the risk value.

## 5. Surface accessibility — two traps

**Trap 1: unannotated is not small.** Annotated outward-facing residue count
comes from the proteome's topological features. Where a protein has no such
annotation the component is `NOT_MEASURED`. It must not be imputed to a midpoint:
doing so puts a protein nobody annotated (0.500) above one measured and found to
have a 10-residue ectodomain (0.420), which inverts the actual finding. The
reference run had **365 proteins** on this path (228 multi-pass, 137
single-pass); this run must report its own count.

**Trap 2: the anchored-without-a-segment class reads zero by construction.**
Proteins held by a lipid anchor have no transmembrane segment to annotate
topological domains around, so they report **zero annotated outward residues**
however large their real ectodomain. MSLN, CEACAM5 and PSCA are all in this
class. Anything that ranks on annotated ectodomain size scores the entire class
at zero unless it special-cases them.

Therefore C5 is computed as:

- for the lipid-anchored class: from confirmed plasma-membrane localisation
  alone, with the residue term excluded — not set to zero
- otherwise: from annotated outward residues, saturating at 200, combined with
  confirmed plasma-membrane localisation
- `NOT_MEASURED` where neither residue annotation nor a localisation call exists

Absence of a plasma-membrane call is **not refutation**. Only ~13,000 genes have
been imaged at all. It confirms a surface call; it never rejects one.

## 6. Normal tissue risk — a gate, not a term

`risk = max over organs of (expression_score(organ) × criticality_weight(organ))`

**Maximum, not mean.** One unmanageable tissue is disqualifying however clean
the other sixty-seven are. A mean lets breadth of safety dilute a single fatal
organ.

### Criticality tiers (platform defaults)

| tier | weight | organs |
| --- | --- | --- |
| 1 | 1.0 | brain, heart, lung, liver, kidney, pancreas |
| 2 | 0.6 | gi_tract, marrow_and_blood, bladder, endocrine, muscle, nerve |
| 3 | 0.3 | skin, adipose, breast, reproductive, salivary |

**Extensions.** The table above covers seventeen organs; the two tissue
vocabularies between them name several more. These assignments are additions
made here, not carried from any prior run, and are flagged for review:

| organ | tier | reasoning |
| --- | --- | --- |
| vascular (aorta, coronary, tibial artery) | 1 | on-target attack on large vessels is not survivable |
| eye (retina) | 2 | not life-threatening, but irreversible |
| mucosa (nasopharynx, oral) | 2 | aerodigestive lining, regenerates |
| connective (cartilage, soft tissue) | 3 | structural, tolerant of damage |

Where a label was genuinely ambiguous the higher tier was taken. Understating
risk is the failure that costs a patient; overstating it costs a candidate.

Cultured cell lines in the baseline vocabulary (cultured fibroblasts,
transformed lymphocytes) are **excluded from the risk computation entirely**.
They are not normal tissue and their expression says nothing about what a
therapy would encounter. The exclusion is reported rather than silent.

This indication overrides pancreas to tier 2 through the project config, with a
required rationale. **The rationale travels into the output header**, so a reader
sees which safety default was relaxed and why. An override without a rationale
fails schema validation — already enforced by the input model.

### Combining the two measurements per organ

`expression_score(organ) = max(atlas_score(organ), baseline_score(organ))`

A "not detected" staining call **must not zero out a positive transcript
measurement for the same organ.** Taking the minimum, or letting the atlas
override, understates risk — the dangerous direction. Where only one source
covers an organ, that one is used. Where neither does, see below.

### Undefined risk is not low risk

If risk cannot be computed — no staining and no baseline row anywhere — the
protein **fails the gate**. It is not scored as zero risk and it is not silently
cleared. This is the single mechanism that handles the unmeasured-protein
problem described in §8; it needs no blocklist.

### Three tissue-mapping bugs that will recur

These are not hypothetical. Each was found the hard way and each will come back
if the mapping is written with naive substring matching.

1. **`"renal"` also matches `"adrenal gland"`.** This placed the adrenal at
   kidney criticality for 2,228 targets in the reference run. Organ keywords
   must be matched on token boundaries, not as bare substrings.
2. **`"cortex"` also matches `Kidney_Cortex`.** This folded kidney into the brain
   bucket. Same fix; `cortex` alone is never a sufficient key.
3. **Six atlas labels map to nothing** and fall through to the tier-1 fail-safe:
   `cartilage`, `dorsal raphe`, `nasopharynx`, `soft tissue 1`, `soft tissue 2`,
   `sole of foot`. In the reference run 9,807 of 13,468 stained genes carried a
   Low-or-better call in at least one of them, so the fall-through was not rare.
   All six require explicit assignments.

**How all three are prevented rather than patched.** Both vocabularies are
finite and fully enumerable — 68 labels on the baseline side, 64 on the atlas
side. Every label is therefore assigned **explicitly, by exact match**, and no
keyword or substring test is used anywhere in the mapping. Substring matching is
what produces all three failures above; removing it removes the class. A label
absent from the table does not fall back to anything — it is counted as a
fall-through and reported.

**The output must report the fall-through count, and it must read 0.** A
fail-safe that is silently doing work is a mapping bug wearing a disguise.

### Clearance

`cleared = risk <= normal_tissue_risk_ceiling` from the project spec: 0.15
conservative, 0.35 moderate, 0.60 permissive. This indication runs conservative.

## 7. Evidence confidence — a third independent number

Reported alongside risk and composite, **never combined with either**.

Confidence reflects how much measurement stands behind a protein: which sources
resolved, how many components were measured, and the summed measured weight `W`.

The reason it stays separate is measurable within a single evidence class. In
the reference run, among PROTEIN_CONFIRMED alone, 1,436 proteins stain High
somewhere in normal tissue and 17 are never detected anywhere. Identical
evidence class, opposite risk. One number cannot carry both.

## 8. Immune-lineage proteins — no blocklist

Measured immune proteins sink themselves under a worst-organ risk score. In the
reference run HLA-B peaks at 8,754 transcripts across all 68 baseline tissues and
CD74 at 5,067 across 63. No special handling is needed and none is added.

The proteins that need handling are TRA, TRB, KIR2DL2 and HLA-DRB3, which have no
baseline row and no staining. They cannot be sunk by a measurement that does not
exist. **The real category is "unmeasured", not "immune"** — which is exactly why
framing it as an immune problem would miss the olfactory receptors, taste
receptors and retroviral envelopes sitting in the same position. Handled entirely
by §6's rule that undefined risk fails the gate.

A blocklist would also be a form of tuning: it encodes the answer rather than
measuring it.

## 9. Output and reproducibility

Every output opens with a header carrying:

- weight set and version identifier
- every threshold in this document
- the criticality table, plus any overrides with their rationale text
- release pins for all six sources
- universe size
- repository revision, including a dirty flag
- a configuration hash

**The configuration hash must cover the tissue keyword lists, not only the tier
assignments.** Adding one keyword moves thousands of risk scores. If the hash
does not move with it, two materially different experiments compare as the same
one, which is worse than having no hash.

The hash must be verified stable across processes — a hash seeded by anything
process-local silently defeats itself.

## 10. Rejection criteria — eleven, fixed in advance

A tripped criterion means the spec changes and the run repeats. It never means
the result gets an explanation.

| id | criterion |
| --- | --- |
| R1 | all five known targets for this indication fall outside the top decile of their tier |
| R2 | a ubiquitous immune protein (HLA-A/B, CD74, PTPRC) clears the ceiling |
| R3 | CEACAM5 fails to rank despite the dropout rule |
| R4 | a protein rejected by the surface filter appears in the ranking |
| R5 | the composite correlates above 0.95 (Spearman) with any single component |
| R6 | a ±20% weight perturbation reshuffles more than half of the top 50 |
| R7 | the risk gate blocks everything, or blocks nothing |
| R8 | the composite distribution spikes at one repeated value |
| R9 | the top of DATA_INSUFFICIENT scores above the top of PROTEIN_CONFIRMED |
| R10 | more than 10% of the top 100 were reached only through the identifier bridge |
| R11 | a `sources_disagree` flag goes unread in the top 25 |
| R12 | a 2× change in any single saturation point reshuffles more than half of the top 50 |

**R12 detail.** R6 perturbs the weights, which were fixed deliberately and are
visible as choices. The saturation points in §4 are equally free parameters and
equally capable of driving the ordering, but they look like implementation
detail rather than policy, so nothing would otherwise test them. Each is
perturbed independently, both doubled and halved, holding everything else fixed:

| parameter | value |
| --- | --- |
| C1 expression saturation | 100 transcripts per 10k |
| C2 specificity ratio saturation | 50× |
| C3 fold saturation | 64× |
| C4 prevalence threshold | 10 transcripts per million |
| C5 ectodomain saturation | 200 residues |
| C6 gene effect saturation | 1.0 |

The criterion trips if any single perturbation replaces more than half of the
top 50. Same shape and same threshold as R6, so the two are directly comparable.

**R5 detail.** Pair components **by name**, and compute the correlation **only
over targets that measured that component**. Zero-filling unmeasured components
inside the check that polices the composite would be the exact imputation this
document forbids, committed by the auditor.

**R10 detail.** "Bridge" means the symbol lookup failed and the row was reached
only through an identifier. The normal baseline is keyed by identifier — that is
its primary key, not a bridge, and must not be counted as one. This is why join
routes are recorded at join time (step 11) rather than reconstructed afterwards.

### Explicitly not grounds for rejection

- MSLN not ranking first
- a novel protein outranking a known one
- the weak tier containing nothing promising

These are possible correct answers, not failures. Treating them as failures is
how a screen gets tuned into a confirmation of what was already believed.

## 11. Expected results (reference run, for comparison only)

Not targets to hit. The universe differs by 16 proteins, so exact reproduction is
not expected and would be mildly suspicious.

| gene | rank (of 1,945) | composite | risk |
| --- | --- | --- | --- |
| CEACAM5 | 1 | 0.894 | 0.600 |
| CEACAM6 | 6 | 0.709 | 1.000 |
| CLDN18 | 8 | 0.678 | 0.667 |
| MSLN | 10 | 0.655 | 1.000 |
| MUC1 | 56 | 0.512 | 1.000 |

Clearance by ceiling, all / PROTEIN_CONFIRMED of 1,945:

| ceiling | all | protein confirmed |
| --- | --- | --- |
| 0.15 conservative | 547 | 1 |
| 0.35 moderate | 834 | 60 |
| 0.60 permissive | 1,332 | 269 |

## 12. Known open problem — carried forward, not fixed here

**The risk gate is miscalibrated and selects for absence of evidence.**

At a 0.15 ceiling the reference run cleared 546 RNA_SUPPORTED proteins and 1
PROTEIN_CONFIRMED. Only 0.2% of cleared targets had staining, against 67.4% of
blocked ones. Mean confidence was 0.49 among cleared and 0.86 among blocked.
Having been stained at all was effectively disqualifying.

**Mechanism.** The atlas's four-level ordinal is mapped onto the same 0–1 axis as
log-scaled transcript level. 1,436 of 1,945 stained proteins peak at "High", and
High × tier 1 = 1.000. 1,095 targets sit at exactly 1.0 and the median risk is
0.67, so a ceiling of 0.15 sits far below the bulk of the distribution and cannot
discriminate within it. It excludes everything anyone has looked at.

**The scores are directionally right, which is why this is calibration and not
collapse.** CEACAM5's peak risk organ is the gastrointestinal tract, and CEACAM5
trials caused severe colitis. MSLN's is lung, and mesothelin sits on mesothelium.
CEACAM6's is liver, MUC1's is kidney. It reproduces real, clinically observed
on-target toxicity. It is not inventing risk; it is failing to grade it.

**Two specific defects.**

1. An immunohistochemical intensity call is being treated as ratio-scale
   abundance. "High" is not maximal danger, and it is the modal value.
2. `peak_level()` takes the maximum over every tissue **and every cell type**, so
   staining in any one cell type of an organ scores as staining in that organ.

**Intended fix, not applied in this stage:** recalibrate the atlas level mapping
against the baseline transcript distribution so the two axes are commensurable.

**Do not move the ceiling.** That is fitting a threshold to a broken scale, and
it would hide the defect rather than correct it.

---

## Build note

Implementation order once approved: `stages/stage3.py`, then
`verify_ranking.py` checking all eleven criteria, and only then any reading of
the biology.
