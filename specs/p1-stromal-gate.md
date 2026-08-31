# P1 — the stromal component becomes a gate

Written before implementation. The rule and its rejection criteria are fixed
here; the implementation and its run follow in a separate commit.

## The defect

`malignant_vs_stroma` carries 0.20 of the weight in a renormalised mean. A mean
cannot disqualify. LRRC15 scores 0.0 on it — the worst value the component can
take — and survives as the only conventional design in the pool, because five
other components outvote it. The component is working correctly and has no
authority.

This is the same shape as risk before it was a gate: a quantity that should
decide admissibility expressed as a term that can be averaged away.

## Which rule, and why not the other

Two candidates were considered.

**A. `malignant_vs_stroma` becomes a gate.** The component is
`log10(malignant / stromal_peak) / log10(50)`, clamped below at zero. The clamp
means a score of exactly 0.0 is not a low score, it is the statement *malignant
does not exceed the strongest stromal compartment*. The natural boundary is the
ratio at which the two are equal.

**B. A floor on malignant expression.** A target at 0.0043 transcripts per 10k
in malignant cells is not a tumour target whatever else it scores, which is
true. But it requires choosing a magnitude, and no magnitude in this pipeline is
already justified for that purpose. `DROPOUT_EPSILON` is 0.001 and marks where
measurement stops, not where biology starts; LRRC15 sits 4.3x above it and is
still a stromal marker. Any floor would be a number chosen to sit between
LRRC15 and the real targets, which is fitting.

**Rule A is adopted.** It needs no new constant. Measured afterwards and
recorded here for completeness: every pool member that passes rule A already
carries a malignant expression of at least 0.0212, against LRRC15's 0.0043, so
rule B would reject nothing that rule A does not already reject. The floor is
not merely harder to justify, it is redundant.

## The rule

> A target is **stroma-dominant** when its malignant-to-stromal ratio is
> measured and does not exceed 1. A stroma-dominant target is not a tumour-side
> candidate and does not enter the screened pool.
>
> The gate fires only on a measurement. Where the ratio is unmeasured, for any
> reason, the gate does not fire and the target is unaffected.

**The threshold is 1.0 because that is where malignant expression equals
stromal expression.** It is not a tuned parameter. It is the point the component
already clamps at, and it was chosen before any pool count was taken.

**Robustness, measured after the rule was fixed.** The nearest legitimate target
to the boundary across both indications is SLC39A6 in breast at 9.0, and ERBB2
at 10.8. LRRC15 sits at 0.105. The boundary could be placed anywhere between
0.11 and 9.0 — a factor of eighty — without changing a single verdict that
matters. A rule whose answer is stable across two orders of magnitude is not one
fitted to exclude a particular gene.

## Interaction with the dropout rule

This is the part that decides whether the change is safe, and it is where the
two kinds of zero must not be confused.

`_score_c2` returns unmeasured in three circumstances: no single-cell row exists
for the gene, the malignant mean is at or below `DROPOUT_EPSILON`, or the
stromal peak is. A measured ratio of 0.1 means *we looked and the gene is
ten times higher in stroma*. An unmeasured ratio means *we did not see it*, and
the capture chemistry of a 3-prime single-nucleus assay guarantees that some
genuinely expressed genes will not be seen.

**The gate must therefore fire only on the measured case.** An absent
measurement must leave a target exactly where it was, never reject it.

The proof cases are already in the pool and are pinned by name in the criteria
below:

| target | state | note recorded | must |
| --- | --- | --- | --- |
| CEACAM5 | unmeasured | `below capture threshold` | survive |
| MUCL3 | unmeasured | `no row` | survive |

CEACAM5 is the case that proves it. Its malignant mean is 5.31e-05, below the
capture floor, so both single-cell components are unmeasured. Under a gate that
treated absence as failure it would be rejected — the highest-ranked target in
the reference run, removed for a dropout. Under this rule it is untouched.

## What it does to the pool

Measured on the current state, after the rule was fixed:

| | pancreatic | breast |
| --- | --- | --- |
| pool | 200 | 200 |
| ratio measured | 198 | 199 |
| **rejected, ratio <= 1** | **47** | **47** |
| unmeasured, exempt | 2 | 1 |

The rejected set is what the rule was written for. In the pancreatic pool it
contains HLA-DQB2, HLA-DQA1, HLA-DQA2, HLA-DPB1, HLA-DRB1, HLA-DRB5 — six MHC
class II genes — together with IGHG2, CD52, CDH11, ADAM12, ANTXR1, PLXDC2 and
LRRC15. That is the immune and fibroblast contamination the earlier degradation
finding described, and it is being removed by the component built to remove it.

Known targets survive in both indications, all far from the boundary:

| | pancreatic | breast |
| --- | --- | --- |
| CEACAM6 | 193.4 | 138.7 |
| MUC1 | 34.8 | 49.3 |
| CLDN18 | 43.3 | not in pool |
| MSLN | 19.7 | 26.5 |
| ERBB2 | — | 10.8 |
| TACSTD2 | — | 12.5 |
| CEACAM5 | unmeasured, exempt | unmeasured, exempt |

**If the implementation removes CEACAM5 or MSLN, the rule is wrong and must not
ship.** That is criterion G2 below.

## Rejection criteria

Written before the implementation runs.

**G1 — the gate never fires on an absent measurement.** No target whose
malignant-to-stromal ratio is unmeasured may be rejected. Trips if any is.

**G2 — positive pin on the known targets.** CEACAM6, MUC1, CLDN18 and MSLN must
survive the gate in the pancreatic pool, and CEACAM6, MUC1, ERBB2, TACSTD2 and
MSLN in the breast pool. CEACAM5 must survive in both, by the exemption rather
than by its ratio. Trips if any is rejected.

**G3 — negative pin, the gate must actually subtract.** LRRC15 must be
rejected, and so must at least one MHC class II gene and one immunoglobulin
gene present in the pool. A gate that rejects nothing is not a gate, and this
half is what makes G2 more than a tautology.

**G4 — the gate must not empty the pool.** The rejected count must be neither
zero nor the whole measured population.

**G5 — every rejection is attributable to a measurement.** The number of
rejected targets must equal the number whose ratio is measured and at most 1,
exactly. No rejected target may carry the note `no row` or
`below capture threshold`. This is the criterion that would catch the gate
firing on missing data rather than on a measurement, and it is stated as an
equality rather than an inequality so that a rejection arriving by any other
path is visible.

**G6 — the gate changes tumour-side admissibility and nothing else.** No
target's `normal_tissue_risk`, `cleared` flag or component values may differ
before and after. The gate decides what is a tumour antigen; it does not decide
what is safe, and the two must not be allowed to leak into one another.

## What this does not do

It does not touch pairing, the risk gate, or the ceiling. It does not change
any component's value. It removes targets from the tumour-side candidate pool
on the grounds that they are not tumour-side candidates.

The consequences for pairing are expected and are P1's second half: LRRC15 is
currently the partner in 72.0% of dual recommendations and P13 trips on it. If
LRRC15 leaves the pool, that hub cannot form. Whether P12 and P13 clear as a
result is a measurement to report, not a claim made here.
