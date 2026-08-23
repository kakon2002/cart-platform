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
