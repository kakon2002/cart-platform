# The staining gate blocks on presence, not on amount

A decision for whoever owns the safety tolerance. Every figure below was printed
by the pipeline's own checks on the reference state. Nothing here is estimated,
and nothing here is a recommendation.

---

## What the pipeline does today

Normal-tissue risk for a target is scored from two sources: **antibody staining**
of normal tissue, graded Low / Medium / High, and **bulk transcript** levels in
the same tissues. Each organ takes the higher of the two, and a target is
blocked if its worst organ exceeds a risk ceiling of 0.15.

Organs are grouped into three criticality tiers, weighted 1.0, 0.6 and 0.3, so
that the same amount of antigen matters more in the brain than in skin.

## The finding

**A Low staining call blocks a target exactly as hard as a High one.**

The three staining grades are calibrated to transcript-equivalent values, and
those land at scores of 0.288, 0.379 and 0.460. Against a 0.15 ceiling:

| | Low | Medium | High | effect |
| --- | --- | --- | --- | --- |
| **tier 1** (brain, heart, lung, liver, kidney, pancreas, vascular, eye) | 0.288 | 0.379 | 0.460 | all block |
| **tier 2** (gut, marrow, bladder, endocrine, muscle, nerve, mucosa) | 0.173 | 0.227 | 0.276 | all block |
| **tier 3** (skin, adipose, breast, reproductive, salivary, connective) | 0.086 | 0.114 | 0.138 | **none can block, at any grade** |

Three consequences follow, and all three are measured rather than argued.

**The grades do no work.** Within a tier, Low, Medium and High are
interchangeable. The outcome depends only on which tier the organ sits in and on
whether staining was detected at all. A three-by-three grid collapses to a
single tier-level rule.

**Nearly half of all blocks come from the weakest evidence.** Of the 72 targets
that pass on transcript alone and fail once staining is added, **34 — 47.2% —
are blocked by a Low call.** A Low call is the least reliable grade an
antibody-based assay produces.

**Tier 3 cannot contribute.** No staining grade reaches the ceiling in a tier-3
organ. The tier structure is currently carrying all of the discrimination and
the grades none of it.

**The calibration does not change who clears.** Because the grades sit on the
same side of each tier threshold, any assignment of transcript-equivalent values
to Low, Medium and High produces the same set of cleared targets, provided the
assignment keeps them on the same side. The calibration exercise, and the check
that defended it, are downstream of nothing that gates.

---

## The price list

The behaviour is controlled by a single scale parameter. Below is what each
setting buys and costs. **No value is recommended here.** The current setting is
marked.

| scale | grade scores | grades that block, by tier | targets clearing | of which protein-confirmed | of the 72 blocks, still blocked | grades separable? |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 0.830 / 1.000 / 1.000 | all / all / all | 468 | 0 | 19 | no |
| 50 | 0.506 / 0.665 / 0.808 | all / all / all | 519 | 0 | 39 | no |
| 60 | 0.484 / 0.636 / 0.773 | all / all / Med+High | 526 | 0 | 40 | yes |
| 100 | 0.431 / 0.567 / 0.689 | all / all / Med+High | 548 | 0 | 46 | yes |
| 250 | 0.360 / 0.473 / 0.575 | all / all / High | 584 | 2 | 55 | yes |
| 500 | 0.320 / 0.421 / 0.511 | all / all / High | 608 | 3 | 62 | yes |
| **1000** | **0.288 / 0.379 / 0.460** | **all / all / none** | **646** | **11** | **72** | **no — current** |
| 2000 | 0.262 / 0.344 / 0.418 | all / all / none | 678 | 15 | 72 | no |
| 5000 | 0.234 / 0.307 / 0.373 | all / Med+High / none | 716 | 23 | 69 | yes |
| 10000 | 0.216 / 0.284 / 0.345 | all / Med+High / none | 748 | 26 | 69 | yes |
| 50000 | 0.184 / 0.242 / 0.294 | all / High / none | 827 | 49 | 65 | yes |
| 100000 | 0.173 / 0.227 / 0.276 | all / High / none | 845 | 55 | 65 | yes |
| 577000 | 0.150 / 0.197 / 0.240 | all / none / none | 946 | 85 | 57 | no |

Read the last column as: does any tier put two grades on opposite sides of the
gate, so that the grade matters at all. The current setting sits in a gap
between two ranges where it would.

Two features of the table are worth noticing. Making the grades matter is
possible in **two directions**, and they are not equivalent: tightening the scale
towards 60–500 makes the gate harsher overall — as few as 0 protein-confirmed
targets clear — while loosening it towards 5,000–100,000 makes the gate more
permissive and raises protein-confirmed clearance from 11 to between 23 and 55.
And the number of blocks removed is modest in either direction: even at 100,000,
65 of the 72 are still blocked.

---

## Does the staining arm earn its place?

Four cases where a directed therapy caused observed toxicity in a known organ.
**This is four data points. It cannot support a threshold and is not used as
one.** It is reported because it bears on the question and nothing else
available does.

Since the gate turns on presence, the question is coverage: is the antigen
detected in the organ where toxicity occurred?

| case | toxic organ | detected by staining | detected by transcript | which arm ranks it worst |
| --- | --- | --- | --- | --- |
| Zolbetuximab, CLDN18.2 | stomach | yes, High | yes, 551 TPM | neither — both rank lung first |
| Mesothelin trials | serosal, sampled with lung | yes, High | yes, 80 TPM | both rank lung first |
| ERBB2, the 2010 case | lung | yes, Medium | yes, 49 TPM | **staining only** |
| CD19 | marrow and blood | yes, High | yes, 183 TPM | neither — both rank gut first |

**In all four, both arms detect the organ.** There is no case here where
staining finds a toxic organ that transcript misses entirely. In one case,
ERBB2, the staining arm ranks the toxic organ as the worst while the transcript
arm ranks kidney first — so staining points at the clinically correct organ and
transcript does not.

On this evidence the staining arm changes the ordering in the right direction
once, and adds no coverage the transcript arm lacks. That is a weak positive,
on four cases, and should be read as such.

One limitation is structural rather than statistical: the organ vocabulary has
no separate category for serosal tissue, which is sampled together with lung.
The mesothelin row therefore cannot distinguish the two, and that is a known
open question rather than a result.

---

## What this is, and what it is not

**It is not a defect.** Under a conservative safety tolerance, vetoing on any
detectable antigen in a critical organ — regardless of how much — may well be
the correct design. A Low staining call is weak evidence of *how much* antigen
is present, but it is still evidence that the antigen *is* present, and for a
therapy that kills what it binds, presence may be the right thing to gate on.

**It is a decision about tolerance, and it belongs to whoever owns that
tolerance.** The pipeline now reports it as an open criterion rather than
silently accepting it, so it cannot be passed over. The internal check that
flags it is failing deliberately and will continue to fail until the tolerance
question is answered. That failure is not a broken build.

The three things a decision needs are all above: what the current setting does,
what any other setting would do, and the only external evidence available on
whether the staining arm is pulling its weight.
