# P0 — evidence asymmetry

Written before the implementation and committed before any result, because the
bound and the ordering rule below decide what the run is allowed to show.

The client's four requirements are taken in order. Two are satisfiable as
written, one needs its rule stating in advance, and one cannot be satisfied by
the other three. That last point is stated here rather than discovered in the
report.

## Requirement 1 — missing protein evidence is UNKNOWN, not zero

Per organ the risk score is `max(staining, baseline)`. Where a target has no
staining call, the staining arm contributes nothing to that maximum, which is
arithmetically identical to contributing zero. The absence is therefore scored
as an absence of risk, and nothing anywhere records that it was an absence of
measurement.

**What changes.** Every target carries an explicit statement of which arms were
measured, and where the protein arm is absent the risk figure is labelled a
lower bound rather than a point value. `evidence_class` already encodes the same
fact — it is `PROTEIN_CONFIRMED` if and only if a staining call exists — so this
adds no new measurement. It makes an implication explicit that was previously
recoverable only by knowing how the scoring function works.

**What does not change.** Clearance. The gate keeps deciding admissibility on
risk, as the client requires, and this requirement is about representation.
The consequence is recorded in requirement 4.

## Requirement 2 — three separate persisted outputs

`normal_tissue_risk`, `evidence_confidence`, `evidence_class`, per target.

Audited before implementation: `evidence_class` is persisted under that name;
`evidence_confidence` is persisted as `target_confidence`; **`normal_tissue_risk`
is not persisted at all** — only `risk_organs_measured`, a count of organs. The
risk value and its peak organ reach the API from an in-memory run and are absent
from the artifact. The gap is closed here.

## Requirement 3 — the ordering constraint

The rule, stated before it is run:

> The ranking orders on the composite computed with unmeasured components scored
> at zero, rather than renormalised away. Equivalently, on
> `composite x measured_weight`.

**Derivation, from the scoring function and not from any measured value.** The
composite is a weighted mean over measured components, divided by the weight
that was measured. Dividing by the measured weight rather than the whole weight
vector imputes every missing component at the mean of the measured ones, which
is favourable whenever the measured ones are above zero. Scoring a missing
component at zero makes it uninformative rather than favourable. That is
requirement 1 applied to the attractiveness axis, and it is the same principle,
not a second one.

**What the rule commits to, in advance.** A target measured on a fraction `m` of
the evidence must beat a fully measured target by a factor of `1/m` on the part
that was measured in order to outrank it. At the evidence floor of 0.40 that
factor is 2.5x. At 0.55 it is 1.82x. This is the ordering constraint the client
asks for, expressed as a number rather than a preference, and it applies to
every target rather than to a named one.

**Separability.** This is a ranking change and touches no gate. Per-target risk
and per-target clearance are identical before and after it. What moves is the
order, and therefore the membership of the screened pool.

## Requirement 4 — rerun, and what it can and cannot show

The requirement is that evidence class no longer drives clearance. **It will not
be satisfied by requirements 1 to 3, and this is a property of the arithmetic
rather than of the implementation.**

Clearance is decided by the risk gate. Requirement 1 changes representation and
not the gate, by the client's own instruction. Requirement 3 changes the order
and not the gate. Neither moves a single target across the ceiling, so the
clearance rate by evidence class is unchanged by both.

Measured on the reference state, the gap has three separable sources:

| rule | protein-confirmed | RNA-only | ratio |
| --- | --- | --- | --- |
| `max(staining, baseline)` | 0.57% | 42.93% | 75.14x |
| staining alone | 3.17% | 0.00% | infinite |
| baseline alone | 4.31% | 42.93% | **9.96x** |

Only two things move the 75.14x. Refusing to certify a target whose protein arm
is unmeasured, which is what `r13-evidence-fix` does and which the client has
rejected as equalising by removing the comparison. Or removing the upper
envelope, which is the second half of P0 below and is a safety-tolerance
decision the client owns.

The residual 9.96x under baseline alone is not a scoring artefact. Both classes
carry that arm, no maximum is taken, and no arm is missing. It is the population
difference the R13 withdrawal documented: median breadth 51 tissues against 7.

## The second half of P0: presence reading as a veto

Requirement 1 concerns absence being read as safety. The mirror is presence
being read as a veto regardless of amount, and it lives in the same
`max(staining, baseline)` machinery.

R14 establishes it from the scoring function alone. In tier 1 the three
calibrated staining scores are 0.288, 0.379 and 0.460 and all exceed the 0.15
ceiling. In tier 2 they are 0.173, 0.227 and 0.276 and all exceed it. In tier 3
they are 0.086, 0.114 and 0.138 and none reaches it. So a Low call blocks
exactly as hard as a High call, tier 3 cannot block at any level, and any
assignment of TPM to levels yields the same cleared set provided the levels stay
on the same side of each tier threshold.

**Not fixed here.** The correction is a change to `BASELINE_TPM_SATURATION`, and
the sweep pricing every setting is already in
`reports/staining-veto-decision.md`. Section 4.4 of the client's document is
explicit that measurements get corrected and thresholds do not get relaxed, and
this is a threshold. It stays in front of him, priced, until he decides.

P0 is not complete with only the first half fixed, and this document is the
record of which half is outstanding and why.
