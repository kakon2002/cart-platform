# P13 — partner concentration, and what actually causes it

Written before any change. The measurements below were taken to answer the
question rather than to justify a fix, and they refute the explanation both of
us had been working from.

## The question

Three genes have now held the hub: NRG3, then NPSR1, then LRRC15, now PRSS21.
The obvious reading is that the selection rule has a systematic preference and
the identity of the gene is incidental. The rule sorts admissible pairs by
combined risk and takes the minimum, and combined risk is the minimum per organ
before the maximum across organs, so a partner that is near absent everywhere
minimises it for every target it is offered to. That is a rule-level cause and
it would need a rule-level fix.

**That explanation is wrong.** It is a coherent story about this scoring
function and it is not what the data does.

## What the measurement shows

**The selection rule is not exercised in the cases that produce the
concentration.**

| | |
| --- | --- |
| pairs evaluated | 19,900 |
| **admissible** (combined risk at or under 0.15, coverage measured) | **290, 1.46%** |
| dual recommendations | 30 |
| **dual targets with exactly one admissible eligible partner** | **22 of 30, 73.3%** |
| targets that chose PRSS21 | 21 |
| **of those, how many had exactly one option** | **21 of 21** |

Every target that took the hub had no alternative. For those 21 recommendations
the sort key never ran: there was one candidate and it was selected because it
was the only one. The distribution of options is 22 targets with one, four with
two, one with three, and three outliers with 15, 19 and 29.

**Consequence, and it decides the shape of any fix.** Even if every target that
*did* have a choice were made to choose differently, PRSS21 would still hold 21
of 30 recommendations — **70.0%**, unchanged, against P13's 50% limit. No
change to the selection objective can move this number.

## PRSS21 is not the profile the old explanation predicts

| | |
| --- | --- |
| normal-tissue risk | 0.3786, peak organ lung |
| clears the ceiling alone | **no** |
| bulk tumour median | 6.39 TPM |
| malignant-to-stromal ratio | 15.78, passes the P1 gate |
| evidence class | PROTEIN_CONFIRMED |
| rank among eligible partners by own risk | **11th** |

It is not the most absent protein. Ten eligible partners carry a lower
normal-tissue risk than it does, headed by TMEM92 at 0.2881. If the rule
favoured absence, one of those would hold the hub.

What PRSS21 has is **supply**: it appears in 26 admissible pairs against seven
for the next best, COL17A1 and CACNG4 and VSIG1. Its organ profile happens to be
complementary to many targets' profiles, so the pair clears where others do not.
That is a property of one gene's biology meeting a scarce admissible set, not a
property of the rule.

## Removing it does not promote an equal successor

Removing the hub and re-deciding, repeatedly:

| removed | duals remaining | top partner | share |
| --- | --- | --- | --- |
| nothing | 30 | PRSS21 | 70.0% |
| PRSS21 | **8** | CACNG4 | 50.0% |
| + CACNG4 | **5** | TMPRSS3 | 40.0% |
| + TMPRSS3 | **2** | CLDN18 | 50.0% |
| + CLDN18 | 2 | PSCA | 50.0% |
| + PSCA | 2 | AMN | 50.0% |

The share stays near half while the denominator collapses from 30 to 2. The
successor is not equally dominant; there is simply almost nothing left to
recommend. **The pool is not the problem and removing genes from it is not a
fix** — it destroys the recommendations rather than redistributing them.

## So what is P13 measuring

P13 asks what share of dual recommendations name the same partner. It cannot
distinguish a rule that prefers one gene from a supply in which most targets
have exactly one admissible partner. Today it is entirely the second. The
criterion is measuring the scarcity of the admissible set through a statistic
about names.

That scarcity is real and is worth a criterion. It is not the one P13 states.

## Options, and what each assumes

**Option 1 — change the selection objective.** Require the partner to
contribute rather than merely to minimise risk. *Measured to be ineffective*:
the objective is not exercised for 21 of the 21 recommendations that create the
concentration, so any objective leaves 70.0%. It also has nothing left to
require. Coverage is span-confounded and was removed from selection for that
reason; tumour expression is already a gate at 5 TPM; the malignant-to-stromal
ratio became a gate under P1. Adding any of these as an objective term changes
at most the nine recommendations that had a choice.

**Option 2 — tighten partner eligibility.** *Assumes the wrong direction.*
Raising the 5 TPM floor shrinks the admissible set further, which increases the
share of targets with exactly one option and makes the concentration worse.
This should be verified rather than asserted if it is ever attempted.

**Option 3 — report the forced choice, and gate on it.** Add a criterion on the
fraction of dual recommendations that had exactly one admissible eligible
partner. Today 73.3%. **A recommendation that names a partner when no
alternative existed is reporting a choice the stage did not make**, and that is
a defect worth failing on, distinct from concentration. This measures the thing
that is actually true and does not require re-fitting P13 to an observation.

**Option 4 — condition P13 on a choice existing.** Restrict the concentration
statistic to targets with two or more admissible eligible partners. That is
eight targets today, which is too small a population to support a 50% limit;
the criterion would become noise. Not recommended.

**Option 5 — accept P13 with a documented rationale**, under section 9. The
rationale would be this document: the concentration is forced by an admissible
set covering 1.46% of pairs, the rule is not exercised in 21 of the 21 cases
that produce it, and no rule change moves the number.

**Recommendation: option 3 together with option 5.** Add the forced-choice
criterion because it states the real property and can fail informatively.
Accept P13 with this document as the rationale, because it is measuring
something true through a statistic that cannot isolate it, and because the
alternative is fitting either the rule or the criterion to a number.

Nothing here is implemented. This document is for review.

## Recorded separately: the atlas composition finding

The P1 stromal gate rejects **1,665 of the 2,565 targets that carry a measured
malignant-to-stromal ratio** — roughly two thirds. Read as a statement about
targets that number is alarming and wrong. It is a statement about the atlas.

A pancreatic ductal adenocarcinoma sample is mostly not tumour cells. The
malignant compartment sits alongside fibroblasts, immune infiltrate and normal
ductal tissue, and a gene expressed anywhere in that mixture at a level
exceeding its malignant expression is stroma-dominant by this measure. Two
thirds of the surface proteome being higher somewhere in the microenvironment
than in the malignant compartment is the expected shape of a desmoplastic
tumour, not evidence that the gate is over-reaching.

The 900 targets whose ratio is unmeasured are a separate population again and
are exempt by construction, not by judgement.

Anyone reading the rejection count later should read it as: the gate removes
what the atlas says is not tumour-dominant, in a tumour whose cellularity is
mostly not tumour.
