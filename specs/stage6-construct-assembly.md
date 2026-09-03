# Stage 6 — construct assembly

Written before `stages/stage6.py` exists. Short on purpose: this stage is
mechanical where 3 and 4 were not. It assembles parts whose choice is already
fixed upstream and reports whether the result fits a budget that is also fixed
upstream. Nothing here is re-derived.

Status: **implemented**. `stages/stage6.py` and its verifier exist and the criteria below were fixed before the run, not after it. The gate this line used to carry has been passed, and leaving it in place would have described the repository incorrectly.

---

## 1. Fixed upstream, not revisited

| constraint | value | source |
| --- | --- | --- |
| construct budget | **3,500 bp** | Stage 1, 4.7 kb vector less 1.2 kb backbone |
| safety switch | **mandatory** | Stage 1, conservative tolerance |
| genetic edits | at most 2 | Stage 1 |
| binder inventory | per target, two routes | Stage 5 |

**And one arithmetic result carried forward rather than rediscovered.** Stage 5
priced the architectures against the mandatory switch:

| design | binder cost | total | vs 3,500 |
| --- | --- | --- | --- |
| single `scFv` | 750 | 2,793 | fits, 707 spare |
| single `VHH` | 360 | 2,403 | fits, 1,097 spare |
| two `scFv` | 1,500 | **3,879** | **over by 379** |
| two `VHH` | 720 | 3,099 | fits, 401 spare |

So a dual design is buildable only from single-domain binders. **Measured across
Stage 5's 720 candidates, exactly one is single-domain**; 612 are Fab and 58 whole
antibody. There is no VHH inventory for this pool.

This stage therefore expects to emit **no dual construct that fits**, and that is
a finding to report, not a shortfall to engineer around. Substituting a smaller
binder that does not exist, or moving the switch out of the payload to make the
sum work, would both be answers invented to satisfy an arithmetic constraint.

## 2. What a construct record contains

Per target, in the Stage 4 order, unchanged:

- **`amino_acid_sequence`** — the translated construct end to end
- **`dna_map`** — nucleotide sequence with **domain boundaries as half-open
  intervals**, partitioning it with no gap and no overlap
- **`parts`** — one row per domain: name, source, accession and residue range
  where it came from a database, length in residues and in bases
- **`budget`** — every part's cost, the total, the headroom, and a verdict
- **`verdict`** — `BUILDABLE`, `BUDGET_EXCEEDED`, or `NO_CONSTRUCT`

The DNA is a **map, not an optimised sequence**. It is produced by reverse
translation under one fixed codon per amino acid, stated in the header, so the
boundaries are exact and reproducible. It is not a codon-optimised ordering
sequence and must not be read as one.

## 3. Domains, and where each sequence comes from

Second-generation architecture. Every human part is fetched by accession and
residue range from the pinned proteome release; nothing is transcribed from
memory.

| domain | source | accession | residues |
| --- | --- | --- | --- |
| leader | proteome | P01732 (CD8A) | signal peptide |
| binder | Stage 5 | — | VH + linker + VL |
| hinge | proteome | P01732 (CD8A) | extracellular stalk |
| transmembrane | proteome | P01732 (CD8A) | membrane segment |
| costimulatory | proteome | Q07011 (TNFRSF9) | cytoplasmic tail |
| activation | proteome | P20963 (CD247) | cytoplasmic tail |
| safety switch | proteome | P62942 + P55211 | FKBP12 and caspase-9 |
| linker, 2A | **synthetic** | — | named literal, no database |

**Two provenance classes and they are not blurred.** A `proteome` part carries an
accession, a residue range and the release pin, and can be re-derived. A
`synthetic` part — the (G4S)₃ linker, the 2A skip peptide — is a designed sequence
with no database entry; it is recorded as a literal with its name and marked
`synthetic`. A part with neither is not a part, and §5 rejects on it.

Binder sequences come from Stage 5 and are used **verbatim**. The heavy and light
variable regions are joined by the linker in that order. No humanisation, no
germlining, no mutation: Stage 5 retrieved them, this stage assembles them, and
anything else would be design presented as retrieval.

## 4. The budget, and what happens when it is exceeded

Cost is the sum of the enumerated parts, three bases per residue, plus a stop
codon. **The sum is printed with its terms**, never as a total alone.

A construct over 3,500 bp is emitted as **`BUDGET_EXCEEDED` with the full
arithmetic and the overage**. It is not silently trimmed, not re-formatted to a
smaller binder, and no part is dropped to make it fit. The stage reports what the
design costs and lets the number stand.

## 5. Targets with no binder

134 of the 200 carry no binder on either Stage 5 route. They are emitted as
**`NO_CONSTRUCT`**, carrying Stage 5's verdict as the reason.

**That is not a failure of this stage and is not counted as one.** It is a
statement about the literature, exactly as Stage 5 §1 says. A construct stage that
reported a design for a target with no binder would have invented the binder.

### 5.1 Only recommendations are built — added during implementation

A binder existing is necessary and not sufficient. **This stage builds only for
targets Stage 4 returned `SINGLE` or `DUAL` for.** A target returned `NO_DESIGN`
or `UNRESOLVED` gets `NO_CONSTRUCT` carrying that outcome as the reason, even
where a binder exists and a construct could be assembled from it.

The rule was not in the first draft of this document and is recorded here rather
than applied quietly. Without it the stage assembles a construct for a target the
pairing stage rejected on risk — MSLN among them, which has four binders and an
outcome of `NO_DESIGN`. A sequence, a DNA map and a budget verdict printed
against that target is a design, and a reader would have to carry the caveat that
upstream says it is not one. That is the wrong way round.

The cost is stated: a target whose safety verdict later changes will need this
stage re-run to obtain the construct its binders already allow.

## 6. Criteria — fixed before the run

Positive pins, not only negatives. The lesson is Stage 5's: every structure-route
check there was a negative, so all of them passed while the route was dead and
returned nothing for all 200 targets.

| id | criterion |
| --- | --- |
| **K1** | **the DNA does not translate back to the amino acid sequence, for every construct** — the round trip is the strongest single check that assembly is correct |
| **K2** | **any assembled construct does not contain its chosen binder's VH and VL verbatim, or nothing is assembled at all** — pinned to MUC16 and MUC17, the only two recommendations carrying a binder on both arms. MSLN and CLDN18 were the obvious pins and cannot serve: both carry binders, but Stage 4 returned NO_DESIGN for the first and a partner without a binder for the second, so neither yields a construct to check. This is the Stage 5 to Stage 6 join, and the equivalent join one stage earlier is the one that failed |
| K3 | domain boundaries do not partition the sequence exactly: any gap, overlap, or interval outside the sequence |
| K4 | any part carries neither an accession with a residue range nor the `synthetic` mark |
| K5 | the printed part costs do not sum to the printed total, or headroom ≠ budget − total |
| K6 | any construct is reported `BUILDABLE` without the safety switch among its parts |
| K7 | a target emits a construct without a binder; or a target with a binder, a `SINGLE` or `DUAL` outcome and — for a dual — a partner binder emits `NO_CONSTRUCT`. §5.1 is why the outcome is part of this |
| K8 | the output row count ≠ 200, or its gene set ≠ Stage 4's pool |

**K1 and K2 are the two that would have caught a dead stage.** K1 fails if any
assembly step corrupts the sequence; K2 fails if the binder never arrives.

### Explicitly not grounds for rejection

- 134 targets emitting `NO_CONSTRUCT` for want of a binder — §5
- a target with a binder emitting `NO_CONSTRUCT` because Stage 4 did not
  recommend it — §5.1
- every dual design exceeding the budget — §1, and it is the expected result
- the DNA map not being codon-optimised — §2, it is a map

---

## Build note

1. `data/domains.py` — parts by accession and range, cached under the usual
   discipline
2. K1 and K2 run before anything else is written
3. `stages/stage6.py` — assemble, cost, verdict
4. `verify_construct.py` — criteria, then the biology

---

## Amendment - K2 re-pinned

**K2's pins were chosen for a property its targets no longer have, and they are
the second stale pin found in this suite, not the first.**

**Why MUC16 and MUC17 were pinned.** Section 6 records the reasoning verbatim:
they were *"the only two recommendations carrying a binder on both arms"*. The
pin was never about those two genes. It was about the one thing K2 exists to
test - the Stage 5 to Stage 6 join across a *dual*, where two independent
binders must both arrive verbatim in one construct. MSLN and CLDN18 were named
in the same row as the obvious pins that could not serve, because neither
yielded a construct to check. MUC16 and MUC17 were what was left.

**Which decision set this is measured on.** K2 reads the persisted Stage 4
artifact through `read_decisions(allow_unusable=True)`. That artifact is written
by `verify_pairing.py:163`, which calls `decide()` with no tolerances, so every
row carries `route_reason="no tolerances supplied; routing disabled"` and the
manifest flags `usable_as_result: False`. The architecture routing therefore
never runs and no ADAPTOR row can exist in it. This is declared, not hidden, and
the verifier opts into it by name - but it bounds what K2 can see, and the bound
is recorded here because the numbers below are meaningless without it.

**Why they no longer hold.** The property has moved off them. Measured on that
artifact:

| | count |
|---|---|
| Stage 4 decisions | 200 |
| outcomes | 167 NO_DESIGN, 30 DUAL, 3 SINGLE |
| genes carrying a binder the assembler can use | 28 |
| duals carrying a binder on **both** arms | **0** |
| duals carrying a binder on their **own** arm | 4 |
| ADAPTOR rows (routing disabled) | 0 |
| constructs assembled | **0 of 200** |

MUC16 and MUC17 both still carry a usable binder on their own arm. What changed
is the arm they are joined to: MUC16 now pairs to CASR and MUC17 to PRSS21, and
neither partner carries a binder, so both rows end NO_CONSTRUCT with *"dual
design, but the partner has no binder"*. The other two duals holding an own-arm
binder are CDH17 to PRSS21 and IL22RA1 to PRSS21 - the same partner three times,
which is the partner concentration already recorded under P13.

Two things empty this stage, and the smaller one is the interesting one. With
routing disabled, only SINGLE and DUAL can build at all, so the adaptor route -
the architecture the platform actually ships for every surviving design in the
worked indication - is absent by construction rather than by outcome. Within
what remains, the 3 SINGLE rows retrieve no binder and all 30 DUAL rows have a
partner that carries none, three of them the same hub partner. Neither cause
alone would leave the stage at zero; both are needed, and only the second is a
property of the pool.

**The finding, which is not the re-pin.** There is nothing to re-pin onto. No
dual in the current pool carries a binder on both arms, so the two-arm join K2
was written to exercise is not exercised anywhere. Picking any pair that does
not test the join would restate the original mistake in fresher genes. K2 now
derives its pin set from the run, through `two_armed_duals()`, and where that
set is empty it **trips and says the join is untested** rather than clearing on
an empty set. The per-construct verbatim checks are unchanged and still run
against whatever assembled.

**Not weakened.** The old K2 also tripped today, through its `nothing was
assembled at all` clause. The amendment does not change whether K2 trips; it
changes what K2 is able to say, and stops it going quietly green the day one
dual assembles while the two-arm join is still uncovered.

**Second stale pin, not the first.** The first was the Stage 4a end-state
assertion recorded at `verify_api.py:94`, which pinned one terminal status and
would have failed once routing sent designs to an architecture that fits - a
criterion testing yesterday's answer. This is the same defect with the sign
reversed: that one would have failed on an improvement, this one would have
passed on a regression. Both come from writing a *result* into a criterion
instead of the *property* the criterion is for. The standing correction is
unchanged and now has two instances behind it: pin the property, derive the
subjects from the run.

**Third stale pin, found by running the suite against the second indication.**
A7 asserted that the top-ranked target returned by `/targets` is `CEACAM5`. That
is true of pancreatic ductal adenocarcinoma and of nothing else; pointed at
invasive breast carcinoma the criterion trips on `CD24`, which is the correct
answer for that indication. A gene symbol written into a criterion can only pass
for one indication on a platform whose whole claim is that it is
cancer-agnostic, and it would also have failed the first time a scoring change
moved the top entry for a good reason.

The property A7 exists for is that the ranked view returns a top entry carrying
a full component breakdown, and that holds for any indication. The gene clause
is dropped and the breakdown assertion is strengthened: the criterion now
derives the expected component set from `stage3.WEIGHTS` rather than counting
six, so adding or removing a scoring component moves the criterion with it
instead of leaving it asserting a stale arity.

All three instances are one shape. **A criterion holding a result rather than
the property it exists to test** fails whenever the result legitimately changes,
and passes whenever it legitimately should not. Which of the two happens is an
accident of sign, not a difference in kind.

**What K2 does not cover, and now says so.** Because routing is disabled in the
set it reads, K2 has never exercised an adaptor construct, and the adaptor is
what the platform returns for every surviving design in the worked indication.
That gap is a property of which artifact the verifier is pointed at, not of the
criterion, so it is not repaired by re-pinning. It is recorded here rather than
changed, because pointing the construct verifiers at a routed decision set
changes what the whole suite verifies and is a scope decision.

**Still open, priced separately, not fixed here.** With 0 constructs assembled,
K1, K3, K4, K5 and K6 all clear on the empty set - K4 prints *"0 parts in the
first construct"*. Their printed text is honest about the zero; their verdicts
are not. That is the same shape as this amendment, but flipping five criteria is
a criteria-design decision rather than a repair, so it is reported and left for
the tolerance call, not taken unilaterally.

---

## Amendment - the decision set is routed, and K0 added

**The call left open above has been taken.** It is specified in
`specs/stage6-routed-decision-set.md` and summarised here, because that document
changes what §6 means.

`verify_pairing` now passes the declared tolerances to `decide`, so the persisted
artifact is the decision set the platform ships rather than one with routing
switched off. Five rows move from `NO_DESIGN` to `ADAPTOR` - FER1L6, GPR35,
TMEM92, TNFSF9, BTNL8 - and five constructs assemble at 2,868 bp, `BUILDABLE`,
10 parts each. No `SINGLE` and no `DUAL` row changes, because the adaptor branch
in `decide` is reached only where neither of those applies.

**§6 gains K0, and five criteria change what they do on nothing.** K0 trips if
the manifest records no routing configuration, if any row carries the
routing-disabled reason, or if nothing assembled. K1, K3, K4, K5 and K6 each
trip on an empty population in their own right. Handed an empty set the stage
now reports 2 of 9 clear where it reported 7 of 8.

**K2 and K7 were restated, not relaxed.** Both assumed every construct's binder
is a Stage 5 sequence candidate. An adaptor's is not - it is the anti-tag part,
retrieved from a deposited structure - so run unchanged against the routed set
each produced five false failures on five correctly assembled constructs. K2 now
checks the binder by the route that supplied it, and for the anti-tag route
asserts the retrieved sequence appears verbatim, that the segment carrying it
declares that provenance and accession, and that the construct names it. It also
gained the partner VL check it was missing. K7 now asks whether every binder a
row's own architecture needs was retrieved.

**K2 still trips**, and for the one honest reason: no dual in this pool carries a
binder on both arms, so the two-arm join remains unexercised. Routing does not
change that and was never going to. What changed is that K2 no longer reports
five spurious failures beside the real one, and that the other criteria now
read 5 constructs instead of 0.
