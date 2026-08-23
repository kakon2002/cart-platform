# MSLN and CLDN18.2 — where the measurements actually stand

Compiled 24 August 2026 from the existing pipeline run. No figure below is new;
each comes from a stage that has already run and been checked against its criteria.

---

## Before the table: one thing that changes how to read it

**The pipeline did not recommend this pair.** Stage 4's own output pairs CLDN18
with LAMP5, and returns `NO_DESIGN` for MSLN — meaning MSLN cannot be made to
clear the safety ceiling by pairing with anything in the pool. The highest-ranked
pair the pipeline produced is NPSR1 + PTPRN2.

MSLN + CLDN18.2 is therefore assessed here as a **conditional candidate proposed
from outside the pipeline**, measured against the same yardstick as the
alternatives. That is the honest framing and it is not a threshold question: the
pair is not close to the ceiling, so relaxing a threshold is not what stands
between it and a recommendation.

---

## The table

| | **MSLN alone** | **CLDN18.2 alone** | **MSLN + CLDN18.2** |
|---|---|---|---|
| **Tumour coverage** | 11.97% of malignant cells | 18.48% of malignant cells ⚠️ | **2.37%** raw ⚠️ |
| *span-matched percentile* | — | — | **75th** for genes this short (17 kb) |
| **Normal tissue risk** (ceiling 0.15) | **0.6366** — peak organ lung | gene-level **0.7271** lung → isoform-resolved **0.5225** gi_tract | gene-level **0.6366** lung → isoform-resolved **0.2277** gi_tract |
| **Verdict against ceiling** | blocked, 4.2× over | blocked, 3.5× over | **blocked, 1.5× over** |
| **Evidence class** | protein-confirmed | protein-confirmed | both members protein-confirmed |
| *numeric confidence* | not persisted by any stage | not persisted | not persisted |
| **Isoform status — expression/risk** | **isoform-robust.** Single promoter; one transcription start site | **resolvable, and resolved above.** Transcript medians exist on the pinned release | **resolvable, and resolved above** |
| **Isoform status — co-expression** | isoform-robust | ⚠️ **isoform-summed, uncorrectable** | ⚠️ **isoform-summed, uncorrectable** |

⚠️ **Every co-expression figure involving CLDN18 is a sum over both isoforms and
cannot be corrected.** The single-cell atlas is 3′-capture single-nucleus
sequencing: it reads a short window at the 3′ end of each transcript, and
CLDN18.1 and CLDN18.2 differ by mutually exclusive *first* exons. They are
identical in the region sequenced. No reference annotation, threshold or added
depth changes this — it is a property of the chemistry. The direction is **not
conservative**: it overstates what an isoform-specific construct would reach,
because cells carrying only CLDN18.1 are counted as hits.

**What isoform resolution did change.** CLDN18.2 carries **1.4%** of the gene's
lung signal (2.05 of 150.96 TPM) and **95.6%** of its stomach signal. Gene-level
scoring therefore put CLDN18's peak risk in the lung — an organ its therapeutic
isoform is barely in. Resolved, the peak moves to the stomach, which is where the
approved CLDN18.2 antibody's dose-limiting toxicity actually is. The pair's risk
falls from 0.6366 to 0.2277, a 2.8× improvement, and **it still does not clear.**

---

## The alternatives, measured the same way

| pair | co-expression | span | span %ile | member risks | combined risk |
|---|---|---|---|---|---|
| CEACAM6 + CLDN18 ⚠️ | **4.51%** | 27 kb | 75th | 0.8135 / 0.7271 | not computed |
| CEACAM6 + MSLN | 2.88% | 13 kb | **94th** | 0.8135 / 0.6366 | not computed |
| **MSLN + CLDN18.2** ⚠️ | **2.37%** | 17 kb | 75th | 0.6366 / 0.5225 | **0.2277** |
| CLDN18 + MUC1 ⚠️ | 1.45% | 16 kb | 63rd | 0.7271 / 0.9215 | not computed |
| CEACAM6 + MUC1 | 0.93% | 12 kb | 77th | 0.8135 / 0.9215 | not computed |
| MSLN + MUC1 | 0.87% | 8 kb | 76th | 0.6366 / 0.9215 | not computed |

*Combined risk was computed only for the pair under discussion. Reporting the
others would require a new run and none was made for this brief.*

**MSLN + CLDN18.2 is not the best-covering pair** — CEACAM6 + CLDN18 reaches
roughly twice as much of the tumour. It is the pair with the **lowest measured
combined risk** of those examined, and the only one whose risk has been
isoform-resolved. Both of its members are among the five targets with clinical
precedent in this indication.

---

## Three caveats that belong with the numbers, not under them

**1. Coverage no longer gates anything, and 2.4% does not mean what it looks
like.** The per-cell detection rate tracks *genomic span* more strongly than it
tracks expression — rank correlation **+0.68** against **+0.20** for bulk tumour
expression, and the effect holds inside every quartile of expression. The cause is
known: the atlas was quantified against a pre-mRNA reference, so intronic reads are
counted and intronic content scales with gene length. Coverage was therefore
removed from partner selection and is now reported only.

So the raw figure and the percentile must be read together. **2.4% reads as
hopeless in absolute terms and sits at the 75th percentile among measured pairs of
comparable gene length.** The pairs that look best on the raw number — 65% and
above — are all pairs of genes an order of magnitude longer, sitting at the 100th
percentile. That is the artefact, not a result.

**2. Stage 4 closed complete-with-limitations: 10 of 14 criteria, four
documented.** Co-expression is barely separable from statistical independence
(the same capture artefact dominates both sides); a ubiquitous immune protein
reaches a cleared pair because the risk gate is organ-level and that protein's
problem is cell-type-level; **48.6% of all clearances depend on an organ nobody
measured for one member**; and the specific partner named is unstable under pool
changes even though the *set* of admissible partners is not.

The conclusion that stage supports is *"a dual design is worth pursuing for this
target"* — **not** *"pair it with this gene"*.

**3. The rescue count is 99, not 33.** Ninety-nine targets that cannot clear the
safety ceiling alone can clear it paired. The earlier figure of 33 was reported
from a run the current code does not reproduce; it should not be quoted again, and
neither should the earlier claims that MSLN + CLDN18.2 ranked first or that the
pair reached 6.1% coverage.

---

## What this means for architecture design

**Construct assembly produces zero buildable designs for this pool**, and MSLN +
CLDN18.2 is not among the two that assemble at all.

The reason is arithmetic and upstream: conservative safety tolerance mandates a
safety switch; the switch costs 1,311 bp; two single-chain binders plus that
switch reach 3,837 bp against a 3,500 bp payload budget. Single-domain binders
would fit — and across 720 retrieved binder candidates, **one** is single-domain,
on a target the pipeline did not recommend, with a placeholder where its second
chain should be.

The exits, priced: a smaller safety switch would have to free **337 bp**, leaving
it at **74%** of its current size. A vector large enough for the current design is
**5.04 kb** against the 4.7 kb assumed, a **7%** increase.

**Binders are not the obstacle.** MSLN has four clinically-staged binders with
sequences (amatuximab, anetumab and two others) and CLDN18 has eleven, including
the approved zolbetuximab. The obstacle is that this pair does not clear the
safety ceiling, and that the design that would carry it does not fit the vector.

---

## One-line summary

MSLN + CLDN18.2 is the lowest-risk pair we have measured and the only one whose
risk is isoform-resolved, but it sits **1.5× over the safety ceiling**, its
coverage figure cannot be corrected on this assay, and no construct for it fits
the payload budget. It is a defensible candidate to keep open; it is not a result.
