# Two gaps against AI_Pipeline.pdf

Scoped, not built. Every number here is measured from this repository or its
cache; where something could not be measured it says so rather than estimating.

Two claims I formed early were **refuted** by adversarial checking and are
corrected below. Both were the load-bearing ones.

---

## GAP 1 — the architecture table

### First correction: two rows are built, not one

| spec row | status |
| --- | --- |
| tumour-specific single antigen → conventional CAR | **built** — `stage6.py:260` emits `"single, second generation"` |
| either antigen identifies tumour → OR-gated dual | not built |
| both needed → AND-gated | **built** — `stage6.py:250` `"dual, split signal, two receptors"` |
| antigen + normal-tissue exclusion → AND-NOT / inhibitory | not built |
| serious normal-tissue expression → switchable or adaptor | not built |
| heterogeneous tumour → tandem / bicistronic / OR | not built |
| target changes over time → modular universal adaptor | not built |

The conventional CAR exists and **fits the budget comfortably**. Exactly one
target ever reaches it: `stage4.py:649` emits `SINGLE` only when a gene clears
the risk ceiling *alone*, and only NPSR1 does — which then has no binder. So the
conventional row is not blocked by the budget or by parts. It is starved of
input.

### The actual blocker is an ordering inversion

The spec selects an architecture **from** a target's risk profile. This pipeline
applies a fixed risk ceiling of 0.15 at Stage 3, **before any architecture is
known**, and only then pairs and assembles. The gate is architecture-blind.

Measured over the 200-member pool:

```
NPSR1   0.0277   <- the only one under the 0.15 ceiling
CD207   0.2272   <- next best. A 0.20 cliff.
MUC17   0.3382
MSLN    0.6366   (lung)  <- spec row five, exactly
```

| ceiling | pool members clearing |
| --- | --- |
| 0.15 (current) | 1 / 200 |
| 0.30 | 6 / 200 |
| 0.40 | 27 / 200 |
| 0.50 | 66 / 200 |

**199 of 200 are blocked at a gate that does not know what would be built for
them.** Rows four, five and seven exist precisely because they change what risk
is tolerable. Implementing them without making the gate architecture-aware
changes nothing: they would be offered only to targets that already passed.

### Budget arithmetic — corrected

The claim "adaptor fits with 677 bp of headroom" was **refuted**. 677 was an
interior point requiring a 732 bp scFv, not a measurement. Corrected:

- `BUDGET_BP = 3500` (`stage6.py:28`), derived at `stage1.py:107` as
  `vector_payload_limit_kb 4.7 − BACKBONE_OVERHEAD_KB 1.2`. It is a **payload**
  budget, already net of promoter/LTR/WPRE.
- Fixed single-receptor scaffold + mandatory switch + stop codon = **2091 bp**.
- Real scFv range across the 30 pool genes with a usable sequence-route binder:
  **705–765 bp, median 720** (not 714–765; that is only the four inside the two
  assembled duals).

| architecture | bp | vs 3500 |
| --- | --- | --- |
| **adaptor / universal (one receptor)** | **2796–2856** | **fits, 644–704 spare** |
| conventional single | 2796–2856 | fits |
| tandem OR | 3564 / 3624 | over by 64 / 124 |
| rapamycin ON-switch | 3579 / 3588 | over by 79 / 88 |
| ON-switch, chain B anchored | 3705 / 3714 | over by 205 / 214 |
| AND-gate (built) | 3834 / 3894 | over by 334 / 394 |
| AND-NOT, CTLA-4 tail | 3957 / 4017 | over by 457 / 517 |
| AND-NOT, PD-1 tail | 4125 / 4185 | over by 625 / 685 |
| bicistronic OR | 4299 / 4359 | over by 799 / 859 |

The mandatory safety switch is 1362 bp including its own T2A — **39% of the
budget**. Every two-receptor design plus that switch is over. That is the whole
shape of the problem.

**The single-domain escape hatch does not exist.** Exactly 1 of 735 retrieved
candidates is single-domain (Letolizumab, CD40LG); its gene is `NO_DESIGN` and
its sequence carries lowercase residues that `stage6.assemblable` rejects. Every
"…but with a VHH it would fit" line is hypothetical.

### Per-row cost

| row | new parts | reachable by the existing mechanism? | effort |
| --- | --- | --- | --- |
| OR tandem / bicistronic | none | n/a — a layout change | **0.5–1 d** |
| AND-NOT / iCAR | an inhibitory tail | **yes, unchanged** | **~1 week** |
| switchable ON-switch | FRB | yes, but only via a different lookup | **~1 week** |
| adaptor / universal | an anti-tag binder | **no** | **~2 weeks** |

**AND-NOT.** All four candidate tails slice cleanly with the existing
`Topological domain` + `"cytoplasmic"` lookup, verified by running
`domains._feature` against live-fetched records: PDCD1 291 bp (ITIM 221–226,
ITSM 247–251), CTLA4 123 bp (no annotated motif), LILRB1 504 bp (four ITIMs),
BTLA 333 bp (no annotated motif). The parts are the easy half. **The pipeline
has no source of NOT antigens** — Stage 4 selects partners for tumour
*co-expression*, and an iCAR needs a normal-restricted antigen *absent* from
tumour. That is a new selection problem, and it is most of the week.

**Switchable.** `_feature(MTOR, 'Domain', 'FRB')` raises `LookupError` — the
entry has no Domain named FRB. It is reachable as `Region` + `"FKBP1A"`,
residues 2012–2144, 133 aa / 399 bp. Then a hard blocker: **the FKBP12 already
in the build is wild-type**, carrying F at mature residue 36. The mandatory
iCasp9 suicide switch is therefore rapamycin-binding, and the same rapalog that
switches the CAR *on* would dimerise the suicide switch and kill the cell. The
standard fix is FKBP12-F36V — a point mutation, which is neither `PROTEOME` (it
does not match the accession) nor `SYNTHETIC` (it is not a designed literal), so
it breaks the two-provenance model and criterion K4 rejects it. That is a spec
amendment, not a code change.

**Adaptor.** The only unbuilt row that fits the budget, and the only one whose
key part cannot be retrieved by the existing mechanism. Live UniProt queries
returned **0 hits** for anti-FITC, peptide-neo-epitope and monomeric
streptavidin; `data/antibodies` contains no anti-tag entry. It needs a sequence
source outside UniProtKB. The budget it saves is paid instead as a **second
manufactured biologic** on its own CMC and regulatory path.

### The finding that governs all of the above

> Every `DUAL` target failed the 0.15 single-antigen risk ceiling. A
> one-receptor adaptor design is therefore **inadmissible for exactly the
> targets whose budget problem it solves** — unless the gate is made
> architecture-aware first.

Architecture work without the routing rule buys nothing. Routing rule first:
**~2–3 days**, and it is a spec change before it is a code change.

---

## GAP 2 — the two missing score terms

### R (shedding / soluble antigen) — NOT derivable as scoped

**Refuted claim.** I asserted that CEACAM5, the current #1, carries a shedding
signal in the cache and would be demoted. The signal I was reading is
`gpi_anchored`, which is an **anchor class, not shedding evidence**.

CEACAM5's cached row contains **zero** shedding vocabulary — no
shed/soluble/secreted/cleaved/plasma/serum anywhere; one chain with exact
bounds; `"Secreted"` absent. And it can never be detected from what is cached:

> CEACAM5 is released by cleavage of its **GPI lipid anchor**, which does not cut
> the polypeptide backbone, so **no CHAIN feature is ever created**. The spec's
> own `SHED_ANTIGEN` predicate reads chain annotations
> (`specs/stage5-binder-discovery.md:442`) and is mechanistically incapable of
> firing for it. That predicate is also spec-only — grep finds no implementation.

So the "10 of the top 20 are shed" figure is withdrawn as a shedding measure. It
counts anchor class and secreted-isoform annotations, which is a different thing.

**What is genuinely available at zero cost.** `data/hpa/proteinatlas.tsv` is
already cached (31.5 MB) and already carries *Secretome location* (col 58),
*Secretome function* (59), and blood concentration in pg/L (cols 62, 63).
`hpa.py:163` reads three columns from that file and ignores these. The HPA cache
fingerprint carries **no column list**, so adding them costs **zero re-fetch** —
unlike UniProt, where `FIELDS` sits inside the fingerprint and adding `cc_ptm`
or `keywords` would re-fetch all 20,431 entries.

Those columns discriminate where the UniProt text fails: CDH17 scores blank on
all four UniProt signals yet reads 7.5 ng/mL in blood.

**Two limits that decide the shape of the term.**

1. **Coverage.** UniProt signals reach 15.1% of the surface universe, HPA 26.9%,
   either 29.4%. **70.6% would be `NOT_MEASURED`.** Under the project's own
   evidence-floor discipline that must propagate as missing, not as zero.
2. **Wrong population.** HPA blood concentrations are *healthy-donor* plasma.
   MUC16 reads 25,000 pg/L — three orders below CEACAM5 — because CA125 is low
   in healthy people and high in patients. The quantity that would settle R,
   soluble antigen in *patient* serum, is in no connected source.

**Verdict: a three-state flag, not a scored component.** Effort **1–2 days** off
the HPA columns. Add **1 day plus a full proteome re-fetch** if you want
UniProt's `cc_ptm`/`ft_site`/`keywords`.

### A (antigen stability) — derivable, cleanly

The cached `malignant_cells_*.npz` carry a per-cell `untreated` boolean **and a
patient label** for all 64,538 malignant cells, giving patient-level resolution.
(The algebraic route from `group_means.npz` also works exactly — residual 7e-18,
`group_cells` is stored — but is unnecessary.) **11,539 treated malignant cells.**

Measured malignant-compartment means, treated vs untreated (CP10K):

| gene | untreated | treated | T/U |
| --- | --- | --- | --- |
| PSCA | 3.595 | 0.671 | **0.19×** |
| MUC17 | 0.394 | 0.096 | **0.24×** |
| MSLN | 0.670 | 0.690 | 1.03× |
| CEACAM6 | 0.997 | 1.065 | 1.07× |
| TMC5 | 6.643 | 13.209 | **1.99×** |
| MUC16 | 0.819 | 1.754 | **2.14×** |
| CEACAM5 | 0.0001 | 0.0000 | undefined (dropout) |

Effort **1–2 days**. Caveats: "treated" pools every non-untreated regimen; one
cohort; and ~2% of genes need a `NOT_MEASURED` guard where the signal is at the
dropout floor.

### Integration is the real cost — and the spec's formula cannot be used

The spec writes an additive score with negative coefficients. This pipeline
computes a **weighted mean over measured components, renormalised by their own
weight** (`stage3.py:826-834`). Adding negative weights to that breaks it.
Measured, with all six positives at 0.8:

| | measured_weight | composite |
| --- | --- | --- |
| penalties unmeasured | 1.00 | 0.8000 |
| penalties measured at **0.0** | 0.85 | **0.9412** |
| penalties measured at **1.0** (max) | 0.85 | 0.7647 |

**A zero penalty raises the score by 17.6%**, because the negative weight is
subtracted from the denominator too. Break-even is a penalty of 0.80 — below
that, a target is better off *having been measured* on the penalty. Worse, it
voids the evidence floor: two positives plus a maximum penalty reach
`measured_weight` exactly 0.40, pass, and score 0.875 — above the current #1's
neighbourhood. `measured_weight` can also reach exactly 0.0 (division by zero)
or go negative.

**So A and R must enter as positive components in inverted form (`1 − risk`), or
as gates outside the mean. That is a spec decision, and it is the first one.**

### Sensitivity, if R were added anyway

| model | effect on the #1 |
| --- | --- |
| flat subtraction from carriers | CEACAM5 holds #1 until p > **0.0716** |
| 7th positive weighted component | **w = 0.05** — the smallest weight already in the model — drops it to #2; w = 0.15 drops it to #10 |

The ranking is tightly packed: median adjacent gap in the top 20 is 0.0069, and
13 of 19 gaps are under 0.010. Only rank 1 sits on a real shoulder (0.0716).

Base rate check: 71/200 of the pool carry a signal versus 10/20 in the top 20 —
mild enrichment, and n = 20 is too small to call it real. **A shedding penalty
is a broad re-weighting, not a targeted correction.**

### Blast radius

`WEIGHTS` is inside `configuration_hash` (`stage3.py:958`), so every downstream
hash moves: s3 → s4 → s5 → s6/s9/s10/s11.

- `data/stage5/binders.json` invalidation is **enforced** → a ~5-minute,
  200-call retrieval.
- If pool *membership* changes, `malignant_cells_<digest>.npz` invalidates —
  and rebuilding it **requires the 8.3 GB single-cell matrix**, which the
  deployable cache deliberately does not carry.
- **Four of thirteen Stage 3 criteria** need re-specifying: R3, R5, R6, R12.
  R12 hard-codes arity in prose — "Seven parameters, so fourteen perturbations"
  becomes nine and eighteen.

### Fragility surfaced on the way

CEACAM5's #1 rests on `measured_weight = 0.55`. Both single-cell components are
`NOT_MEASURED` — raw malignant mean 5.31e-05, below `DROPOUT_EPSILON` 0.001 —
so its 0.8769 is a mean over 55% of the evidence rescaled as if it were the
whole score. That is independent of anything in this document and is arguably a
larger problem than either missing term.

### Spec drift found

| spec | says | measured |
| --- | --- | --- |
| stage5 §…:440 | 190 multi-chain surface proteins | **189** |
| stage6 §1 | constructs 2,793 / 3,879 bp | **3,834 / 3,894** |
| stage5 §1 | 720 candidates, 612 Fab | **735, 623** |

---

## Summary of effort

| item | effort | note |
| --- | --- | --- |
| **Architecture routing rule** | 2–3 d | prerequisite for all of GAP 1 |
| OR tandem / bicistronic | 0.5–1 d | no new parts; tandem misses by 64 bp |
| AND-NOT / iCAR | ~1 week | parts easy; no NOT-antigen source |
| Switchable ON-switch | ~1 week | FKBP collision forces a provenance amendment |
| Adaptor / universal | ~2 weeks | only row that fits; part not reachable |
| **A — antigen stability** | 1–2 d | derivable now, patient-level |
| **R — shedding flag** | 1–2 d | three-state flag off cached HPA columns |
| R as a scored component | not recommended | 70.6% unmeasured, wrong population |
| Composite integration decision | 1 d spec | negative weights are unusable |
