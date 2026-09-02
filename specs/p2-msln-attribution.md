# P2 — risk attribution

Written before the implementation and committed before any result. The criteria
in §6 are identities rather than thresholds, which is deliberate: there is no
number here to fit, and a bound chosen after seeing MSLN's decomposition would
be a bound chosen to make MSLN's rejection look defensible.

MSLN is the worked example, not the subject. The subject is every target. This
is stated first because the platform is cancer-agnostic and seeds no target, and
an attribution facility written around one gene would be the seeding rule broken
by the back door.

## 1. What the client is actually asking

The platform declines MSLN. MSLN is the most clinically-exercised solid-tumour
CAR-T antigen in the registry the pipeline reads — S1 returns 169 trials for it
— and this pipeline routes it `NO_ARCHITECTURE` at Stage 4a on a measured risk
of 0.6366 in lung, against a declared terminable ceiling of 0.35. A6 already
records that as an accepted, deliberate outcome.

What A6 does not do, and what nobody reading the API can currently do, is say
**which measurement caused it**. That is the request. A reader who disagrees
with the verdict has, today, nothing to disagree with: `risk: 0.6366`,
`risk_organ: "lung"`, and no way to reach the observation underneath.

The honest framing is not "show that the platform is right about MSLN". It is
"make the verdict falsifiable". If the number rests on one staining call in one
tissue, that must be visible, including in the case where seeing it makes the
verdict look weak.

## 2. The defect

Risk is computed and then flattened. From `stage3.py:583`:

```
risk = max over organs of ( per_organ[organ] x tier_weight(organ) )
per_organ[organ] = max over measurements mapped to that organ of
                     staining: calibration.score(level)
                     baseline: log10(1 + tpm) / log10(1 + 1000), clamped
```

Three reductions happen in sequence — a max across measurements within an organ,
a multiplication by a tier weight, a max across organs — and only the final
scalar and the winning organ's name survive into `Ranked`. The evidence endpoint
at `server.py:493` faithfully reports both and cannot report more, because by
then the inputs are gone.

Consequences, all of which apply to every target and are merely most visible on
MSLN:

- A reader cannot tell whether 0.6366 came from the protein arm or the
  transcript arm. These are different kinds of claim with different failure
  modes, and R14 is an open decision about exactly one of them.
- A reader cannot tell whether the winning organ won by a wide margin or by
  0.001 over the runner-up. A verdict resting on a tie is not the same verdict.
- Where the protein arm was never measured, its absence contributes nothing to
  the inner max, which is arithmetically identical to contributing zero. P0
  fixed the labelling of this at the target level. It is not fixed per organ,
  and per organ is where the max actually happens.

## 3. What attribution means here

Attribution is a **reconstruction**, not a commentary. For each target it
records, for every organ that scored, the inputs to the three reductions above,
such that the reported risk can be recomputed from the record alone.

Per organ:

| field | meaning |
| --- | --- |
| `organ` | the organ key |
| `score` | `per_organ[organ]`, before criticality |
| `tier`, `weight` | the consequence tier and the weight it carries |
| `weighted` | `score x weight` — the quantity the outer max ranks |
| `arm` | which arm produced `score`: `STAINING`, `BASELINE`, or `TIED` |
| `staining` | the winning atlas tissue, its level name, and its calibrated score — or `NOT_MEASURED` |
| `baseline` | the winning normal-tissue label and its TPM — or `NOT_MEASURED` |

Per target: the ordered organ list, the winning organ, and `margin`, the
difference between the best and second-best `weighted` values where at least two
organs scored.

`NOT_MEASURED` is a value, not a zero. An organ with no staining call and a
baseline reading carries `staining: NOT_MEASURED` and is not recorded as having
been stained and found clean. This is the same third state P0 establishes, moved
to the level the max is taken at.

## 4. Where it is exposed

The attribution block hangs off the existing evidence trail rather than becoming
a new endpoint: `GET /projects/{id}/evidence/{gene}` gains
`stage3_screen.risk_attribution`. The endpoint already answers per gene, already
carries `risk` and `risk_organ`, and A9 already pins its shape. A separate
endpoint would let the two drift, and the number they disagree about would be
the one under dispute.

The report gains a worked example. Which target is worked is a review-time
choice and MSLN is the obvious one, but the selection lives in the report
generator, never in the pipeline.

## 5. Independence

Risk attribution carries no confidence, no `measured_weight`, no evidence class
and no composite term. Those belong to the other score and the two do not
combine. The temptation here is specific and worth naming: a decomposition that
shows a verdict resting on one measurement invites a confidence-weighted
softening of that verdict, and that is the one move this platform does not make.
Attribution reports what the measurement was. Whether the measurement is
trustworthy is the other score's question, answered separately, and the reader
holds both.

## 6. Rejection criteria — fixed before the run

Written before implementation. None is a fitted threshold; all but T5 and T8 are
identities, which is what makes them worth running.

| id | trips when |
| --- | --- |
| **T1** | for any target carrying a risk figure, `max(weighted)` over the attributed organs does not equal the reported `risk` to within 1e-12 — the attribution does not reconstruct the number it explains |
| **T2** | the named `risk_organ` is not among the organs attaining that maximum, or a tie is reported as a single winner |
| **T3** | for any attributed organ, recomputing `score` from the named raw input — the level through the calibration curve, or the TPM through the baseline curve — does not reproduce the recorded `score`, or `arm` names an arm whose value is not the maximum |
| **T4** | any organ reports a numeric staining score for a target whose protein arm was never measured, or reports `NOT_MEASURED` on an arm that did supply the winning value |
| **T5** | across the run, no target's winning organ is attributed to `STAINING`, or none to `BASELINE` — a positive pin, derived from the run rather than named in advance, because an attribution exercised on only one arm has not been shown to work on the other |
| **T6** | any target with two or more scored organs omits `margin`, or reports a `margin` that disagrees with its own organ list |
| **T7** | the attribution payload contains any confidence, evidence-class or measured-weight field |
| **T8** | any gene symbol appears as a literal anywhere in the attribution code path |

**T5 is the one that would catch a dead facility.** T1 through T4 are satisfied
by an attribution that is correct about every target it describes, including one
that describes none. T5 requires that both arms were actually exercised, and it
derives its subjects from the run rather than naming them, because a criterion
that names its subjects goes stale the moment the pool moves — which is the
defect recorded twice already in `specs/verification-sharing-assumptions.md`.

### Explicitly not grounds for rejection

- MSLN's risk being high, or the verdict on MSLN being unchanged. Attribution
  explains a number; it does not adjudicate it.
- A verdict resting on a single measurement. That is a finding to report, not a
  fault to fix, and the tolerance question it raises is R14's, not this one's.
- The margin being small. No bound is set on it here, because any bound would be
  a clinical policy choice and this spec is not the place to make one.

## 7. What this does not do

It does not change any risk figure, any threshold, any routing outcome or any
ranking. Nothing in §3 feeds back into scoring; the facility is a read of
intermediate state that is currently discarded. If any criterion count or hash
moves, that is a defect in the implementation, not a result.

It does not settle whether declining MSLN is correct. It makes the grounds
visible so that a clinician can settle it.

It does not address why the protein arm behaves as a presence veto. That is R14,
it is open, it is priced in `reports/staining-veto-decision.md`, and the two
interact: if MSLN's lung score turns out to come from the staining arm, then
R14's undecided tolerance is what decides MSLN. Whether that is so is a
measurement, and this spec does not presume the answer.

---

## 8. What the run reported

Implemented as specified. **T1 through T8 all clear on the first run**, and the
counts are worth reading rather than just the verdicts:

| | |
| --- | --- |
| targets carrying a risk figure, every one reconstructed to within 1e-12 | **3,400** |
| attributed organ rows, every one recomputing its own score | **65,077** |
| organ rows carrying `NOT_MEASURED` on the protein arm | **30,454** |
| targets reaching their maximum on more than one organ, reported as ties | **314** |
| winning organs decided by the transcript arm, scoring above zero | **2,769** |
| winning organs decided by the staining arm, scoring above zero | **856** |

**T5 is the one that had to be watched and it is satisfied by a wide margin.**
Both arms decide real verdicts, so the facility has been exercised on both and
is not a decomposition that only ever describes one kind of evidence.

T5 separates winning organs that carry a score from those that do not, and the
distinction is not cosmetic. A protein absent from every organ scores zero
everywhere, so every one of its organs ties for the maximum; counted raw the
transcript arm wins 4,893 organs, of which 2,769 score above zero. Reporting the
raw figure alone would claim more evidence than the run contains.

The 30,454 unmeasured protein arms are the point of §3's third state. Each is an
organ where a transcript reading exists and no staining call does, and each is
now visibly *not measured* rather than arithmetically indistinguishable from
stained-and-clean.

### The worked example, and the answer to §7's open question

MSLN's risk is **0.6366 on lung, ahead of the next organ, liver, by 0.0860**,
and it is carried by the **transcript** arm:

| organ | weighted | score | tier | arm | staining | transcript |
| --- | --- | --- | --- | --- | --- | --- |
| lung | 0.6366 | 0.6366 | 1 | `BASELINE` | bronchus, High, 0.4601 | Lung, 80.3 TPM, 0.6366 |
| liver | 0.5506 | 0.5506 | 1 | `BASELINE` | liver, Not detected, 0.0000 | Liver portal tract, 43.9 TPM, 0.5506 |
| marrow_and_blood | 0.2760 | 0.4601 | 2 | `STAINING` | tonsil, High, 0.4601 | spleen, 0.3 TPM, 0.0389 |

Twenty organs scored. A reader who wants to disagree now has something to
disagree with: a named tissue, a named quantity, and the curve that turned it
into 0.6366.

**§7 named an interaction with R14 and declined to presume it. The measurement
settles it: R14 does not decide MSLN.** The lung score comes from the transcript
arm, and with the staining arm taken alone MSLN would score 0.4601 — still above
the 0.35 terminable ceiling and far above the 0.15 persistent one. Whichever way
the staining-veto tolerance is decided, MSLN's verdict does not move.

That holds for four of the five known targets. Measured per arm, in isolation.
The table below is generated by `make_brief.py` and lands in
`reports/brief-tables.json` under `arm_isolation`, so it moves when the
measurements move rather than standing as prose nobody recomputes:

| target | risk | decided by | staining arm alone | transcript arm alone |
| --- | --- | --- | --- | --- |
| MSLN | 0.6366 | transcript | 0.4601, over 0.35 | 0.6366, over |
| CLDN18 | 0.7271 | transcript | 0.3786, over 0.35 | 0.7271, over |
| CEACAM6 | 0.8135 | transcript | 0.3786, over 0.35 | 0.8135, over |
| MUC1 | 0.9215 | transcript | 0.4601, over 0.35 | 0.9215, over |
| CEACAM5 | 0.6000 | transcript | 0.2881, **under** 0.35 | 0.6000, over |

Every one of the five is decided by the transcript arm today. CEACAM5 is the
only one whose verdict the staining arm could not reach on its own. This is
reported, not acted on: it bears on R14's tolerance question and that decision
is not taken here.

### §7 held

The stage 3 configuration hash is unchanged at `a91c696f2e1318f7`, stable across
processes. No risk figure, threshold, routing outcome or ranking moved; R14
trips exactly as before with the same grid. The only count that moved is this
specification's own eight criteria, taking stage 3 from 18 of 19 to 26 of 27.
