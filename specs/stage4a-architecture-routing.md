# Stage 4a — architecture routing

Inverts the order the pipeline currently applies. Today a fixed risk ceiling is
applied at Stage 3 before any architecture is known, so the architectures that
exist to make risky targets tolerable are only ever offered to targets that
already cleared. This stage makes the ceiling a property of the architecture.

## 1. The inversion

**Before.** One ceiling, applied blind:

```
rank -> risk <= 0.15 ? cleared : blocked -> pair -> assemble
```

199 of 200 die at that gate. Measured: NPSR1 0.0277 is the only pool member
under it, and the next is CD207 at 0.2272 — a 0.20 cliff.

**After.** The risk profile selects a candidate architecture, and the ceiling
that applies is the one declared for that architecture:

```
rank -> route(risk profile) -> architecture -> ceiling(architecture) -> admit? -> assemble
```

## 2. The routing rule

Ordered by **increasing product complexity**. The first architecture that admits
the target wins. This ordering is a stated principle, not a tuned one: a design
that needs one receptor and one product is preferred to one that needs two, and
a design that needs one manufactured product is preferred to one that needs two.

| # | condition | architecture | built? |
| --- | --- | --- | --- |
| 1 | target clears the persistent ceiling alone | `CONVENTIONAL` | yes |
| 2 | a partner exists and the pair clears the persistent ceiling | `AND_GATE` | yes |
| 3 | target does not clear alone, and clears the **terminable** ceiling | `ADAPTOR` | **this stage** |
| 4 | an exclusion marker exists for the target's risk organ | `AND_NOT` | no — §6 |
| — | none of the above | `NO_ARCHITECTURE` | — |

Rows 4 (`AND_NOT`), and the switchable, bicistronic and tandem variants, are
**not implemented**. They resolve to `NOT_IMPLEMENTED` with the reason, never to
a silent `NO_ARCHITECTURE`: a target that would route somewhere unbuilt is a
different finding from a target that routes nowhere, and collapsing the two
would understate what the missing architectures are worth.

## 3. Two ceilings, and why they are two

The current ceiling governs a **persistent** exposure: a living, self-amplifying
T cell that cannot be withdrawn once infused. An adaptor design does not remove
on-target/off-tumour risk — the adaptor still binds the antigen — but it makes
the exposure **terminable**, because activation requires a separately dosed
protein with a finite half-life.

Magnitude and reversibility are different axes, so they get different numbers
and are never blended:

```
persistent_ceiling   the risk a design may carry when it cannot be stopped
terminable_ceiling   the risk a design may carry when it can be
```

**Both are policy inputs, not measurements.** The platform cannot derive a
clinical risk tolerance from expression data, and inventing one would be the
worst kind of quiet parameter. `persistent_ceiling` is the existing
`normal_tissue_risk_ceiling`. `terminable_ceiling` is declared on the project
config beside it, and is **required** — a project that does not declare one
gets `NOT_CONFIGURED` for the adaptor row rather than a default that would
silently set clinical policy.

For this indication it is fixed at **0.35**, before the run, and §5 reports the
full sweep so its effect is visible rather than argued.

### What is scored under an adaptor design

The CAR does not bind the tumour antigen; the adaptor does. Those are two
different molecules with two different exposure profiles, and the honest
treatment is to report both rather than average them:

| | binds | persistence | risk |
| --- | --- | --- | --- |
| receptor | the tag | permanent, self-amplifying | the tag's risk |
| adaptor | the antigen | dosed, finite half-life | the target's risk, unchanged |

The tag is not a human protein, so `receptor_risk` is **not applicable** rather
than zero — the pipeline has no measurement of it and will not manufacture one.

**The target's risk number is never reduced by routing.** It is carried through
unchanged and reported. What routing changes is which ceiling it is compared
against. Scoring the receptor's risk in place of the target's would make the
safety gate vacuous, since every adaptor receptor looks identical and harmless.

## 4. Where it runs

A separate module, `car_pipeline/stages/routing.py`, called from Stage 4's
`decide()`. Stage 3 is not modified: it keeps computing risk and composite as
before, and `Ranked.cleared` keeps its existing meaning of "clears the
persistent ceiling", which several Stage 3 criteria are written against.

Cache consequence, stated before running: the routing policy enters the Stage 4
configuration hash, so `data/stage5/binders.json` is invalidated and re-fetched.
That is correct — `TargetBinders` stores `outcome` and `partner`, both of which
routing changes. Pool membership is **not** affected, because `build_pool`
ignores risk entirely, so the single-cell digest holds and the 8.3 GB matrix is
not required.

## 5. Rejection criteria

Fixed before the run.

| | criterion | trips if |
| --- | --- | --- |
| **A1** | Routing is order-independent | any target's architecture changes when the pool is shuffled |
| **A2** | The simplest admitting architecture wins | any target routed to `ADAPTOR` would also have been admitted by `CONVENTIONAL` or `AND_GATE` |
| **A3** | Risk is never reduced by routing | any routed target's reported risk differs from its Stage 3 risk |
| **A4** | The persistent ceiling still binds | any target admitted as `CONVENTIONAL` has risk above `persistent_ceiling` |
| **A5** | Positive pin — NPSR1 | NPSR1 (risk 0.0277) does not route to `CONVENTIONAL` |
| **A6** | Positive pin — MSLN | MSLN (risk 0.6366, lung) does not route to `ADAPTOR` or above the terminable ceiling |
| **A7** | Unbuilt rows are named | any target resolves to `NO_ARCHITECTURE` when an unimplemented row would have admitted it |
| **A8** | The terminable ceiling is declared, not defaulted | routing admits an `ADAPTOR` for a project with no declared `terminable_ceiling` |
| **A9** | Sensitivity is reported | the run does not emit the admitted count across the ceiling sweep |
| **A10** | Not tuned to an outcome | the declared ceiling differs from the value recorded in this spec |

A9 and A10 exist because §3 fixes a number this pipeline cannot measure. A9
makes its effect visible; A10 makes moving it after seeing output a tripped
criterion rather than an edit.

## 6. What is deliberately not built

`AND_NOT` needs a source of exclusion antigens — a normal-tissue-restricted
surface protein *absent* from tumour. Stage 4 selects partners for tumour
**co-expression**, which is the opposite relation, so no stage emits such a
list. That is a data gap, not a parts gap: every candidate inhibitory tail
(PDCD1 291 bp, CTLA4 123 bp, LILRB1 504 bp, BTLA 333 bp) slices cleanly with the
existing feature lookup. Targets that would route there report
`NOT_IMPLEMENTED: no exclusion-antigen source`.

The switchable row is additionally blocked by a molecular collision: the FKBP12
already in the build is wild-type, so the mandatory rimiducid-inducible suicide
switch and a rapamycin ON-switch would respond to the same drug. Resolving it
needs FKBP12-F36V, a point mutation, which is neither `PROTEOME` nor `SYNTHETIC`
and so breaks the provenance model that criterion K4 enforces.
