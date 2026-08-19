# Stage 5 — binder discovery

Written before `stages/stage5.py` exists and before any binder output exists.
Every source decision, threshold and rejection criterion below is fixed by this
document. If a criterion trips, the correction is a change to this document
followed by a re-run — never a narrative explaining why the output was acceptable
after all.

Status: **awaiting review. Nothing in `stages/stage5.py` may be written until this
is approved.**

Section references of the form "Stage 4 §6" point at the stage's own spec in this
directory.

---

## 0. Preconditions

Two upstream stages are not in a state to be built on, and this one consumes both.

| stage | state | effect here |
| --- | --- | --- |
| Stage 3 | **R13 tripped** — clearance 0.6% protein-confirmed vs 42.8% RNA-supported, ratio 75.32x against a limit of 5x | the `cleared` flag is not yet meaning what it claims |
| Stage 4 | **5 of 16 criteria tripped** — P4, P7, P8, P12, P13 | the `DUAL` recommendations this stage would take as input are not yet a result |

Stage 5 is specified now because the source question in §2 must be answered before
anything depends on it, and because writing the spec first is the rule. **It must
not be run for interpretation until both clear.**

P13 matters here specifically: one protein, **NRG3**, is the recommended partner
for 12 of Stage 4's 13 `DUAL` targets. §2.4 shows NRG3 has no structures at all,
so a binder stage run today would report "no binder available" for almost every
dual design — which would look like a Stage 5 finding and is actually P13.

---

## 1. Scope

Stage 5 answers one question per target: **is there a binding domain that can be
put in a CAR against this antigen, and what is known about it?**

A binder here is the antigen-recognition domain only — not the construct. Formats
are the ones Stage 1 admits: `scFv`, `VH_VL`, `VHH`, `ligand`.

Stage 5 produces, per target:

- **retrieved candidates**, each with format, sequence where available, source
  structure where available, and the antigen chain it was solved against
- **an ectodomain check** — does the epitope sit on the part of the protein a CAR
  can reach (§4.1)
- **a size in base pairs**, against the Stage 1 construct budget (§4.3)
- **a usability verdict per candidate**, and an explicit **`NO_BINDER`** where
  nothing was found

**Stage 5 does not design a binder, does not predict affinity, and does not dock
anything.** It retrieves and characterises what exists. Anything else would be a
prediction presented beside measurements, and the two would be read alike.

**Stage 5 does not re-rank targets.** Stage 4 §7.4 flags binder availability
rather than filtering on it, precisely so that the target list is not shaped by
which proteins happen to have been crystallised. Stage 5 must not undo that by
reordering anything: it annotates, and the annotation travels.

---

## 2. Sources — the first stage needing something not yet connected

Stage 1 declared one blocking dataset for this stage: **"SAbDab and PDB", binder
retrieval, stage 5, required**. It is the only required dataset the availability
score has never been able to satisfy — the 0.857 in Stage 2 is 6 of 7 blocking
sources, and this is the seventh.

**The answer is that half of it is reachable and half of it is not.** Measured by
probe rather than assumed:

### 2.1 What was probed

| endpoint | status | content type | verdict |
| --- | --- | --- | --- |
| RCSB search API (`search.rcsb.org/rcsbsearch/v2/query`) | 200 | `application/json` | **usable** |
| RCSB data API (`data.rcsb.org/rest/v1/core/entry/...`) | 200 | `application/json` | **usable** |
| RCSB FASTA (`rcsb.org/fasta/entry/...`) | 200 | `text/x-fasta` | **usable** |
| RCSB mmCIF (`files.rcsb.org/download/....cif`) | 200 | `chemical/x-cif` | **usable** |
| SAbDab `summary/all/` | 200 | **`text/html`** | **not usable** |
| SAbDab `summary/all/?format=tsv` | 200 | **`text/html`** | **not usable** |
| Thera-SAbDab `search/?all=true` | 200 | **`text/html`** | **not usable** |
| IMGT (`imgt.org`) | — | — | **not usable**, certificate verification fails |

### 2.2 SAbDab is not machine-readable at its documented endpoint

The summary endpoint that is documented to return a tab-separated table returns an
HTML application shell, at both the plain and the `format=tsv` URL. It answers 200,
so a connector that checked only the status code would cache an HTML page and
report the source as available — the failure is silent in exactly the way this
project keeps finding.

This is the same shape as the DepMap portal problem recorded in Stage 2: a source
that a browser can reach and an unattended run cannot. It was solved there by
going to the archived copy rather than the portal. **No equivalent has been
established here, and this document does not assume one exists.**

### 2.3 What is lost with SAbDab and what is not

SAbDab is a curated layer over the PDB. Losing it loses the curation, not the
structures:

| carried by SAbDab | available without it |
| --- | --- |
| which PDB chains are antibody heavy/light | recoverable, imperfectly, from chain descriptions and sequence |
| CDR boundaries under a numbering scheme | **not recoverable without a numbering tool** |
| which chain is the antigen | recoverable from the entity descriptions |
| species of origin | in the PDB entity taxonomy |
| affinity where curated | **not in the PDB** |
| therapeutic name and clinical stage (Thera-SAbDab) | **not recoverable** |

**The two that do not survive are CDR boundaries and affinity.** Both are load
bearing: §4.2 is about affinity and cannot be computed at all from the PDB, and a
sequence without CDR boundaries can still be retrieved and sized but cannot be
humanised, compared or grafted.

### 2.4 Structure counts, measured, and why they decide a rejection criterion

Full-text entry counts from the RCSB search API:

| target | entries | | target | entries |
| --- | --- | --- | --- | --- |
| CEACAM5 | 46 | | PTPRK | 4 |
| MUC1 | 32 | | ABCC3 | 4 |
| claudin-18 | 17 | | PTPRN2 | 3 |
| mesothelin | 15 | | PSCA | 2 |
| EPCAM | 13 | | **NRG3** | **0** |
| CEACAM6 | 6 | | VMP1 | 0 |
| EFNA5 | 5 | | NPSR1, TMC5, MUCL3 | 0 |

The known targets have structures. **Everything Stage 4 discovered has almost
none**, and NRG3 — Stage 4's recommended partner for 12 of 13 dual designs — has
none at all.

**Structure count measures how much attention a protein has had, not how good a
CAR target it is.** Any Stage 5 that scores targets by what it can retrieve will
rank the known targets top and the novel ones bottom, and will do so whatever the
biology says. That is R-null for this stage and it is criterion **B1**.

### 2.5 What this stage can honestly do today

- **Can do:** retrieve antibody-containing structures against a named target from
  the PDB; extract chain sequences; identify the antigen chain and its species;
  determine whether the antigen construct in the structure corresponds to the
  ectodomain; size a candidate in base pairs; report `NO_BINDER` honestly.
- **Cannot do:** report affinity, delineate CDRs, identify clinical-stage
  therapeutics, or say anything about humanisation state beyond species of origin.

**The spec is written to the reachable subset.** Every field that depends on
SAbDab is defined here, declared unavailable, and emitted as `NOT_CONNECTED` — not
as absent, and never as a default. A field silently missing and a field known to
be unobtainable are different states, and Stage 1's `DatasetStatus` already
distinguishes them: SAbDab is `unreachable`, not `not_configured`, because a
connector could exist and the data cannot be read.

**Stage 1's dataset list must be amended** to split "SAbDab and PDB" into two
entries with independent status, because they now have different status and one
row cannot carry both. That is a Stage 1 change, made there.

---

## 3. The pool

**Every target Stage 4 emitted, in all four outcome classes, plus every target in
the Stage 4 pool.** 200 proteins.

Not just the `DUAL` and `SINGLE` recommendations. Two reasons:

- A `NO_DESIGN` target failed on risk or coverage, not on bindability. Its binder
  status is still worth knowing, because a later change to the risk calibration
  (§0) may move it, and re-running retrieval is expensive.
- Restricting to recommendations would let Stage 4's tripped criteria propagate
  silently into what Stage 5 even looks at.

Binder retrieval is per target and independent per target, so there is no pairing
combinatorics here: 200 retrievals, not 19,900.

---

## 4. What makes a binder usable

### 4.1 The epitope must be on the part the CAR can reach

A CAR binds a cell surface from outside. A structure showing an antibody bound to
its antigen proves binding; it does not prove the epitope is **accessible on an
intact cell**.

The check: the antigen construct in the structure must correspond to residues in
the **extracellular topological domain** from the proteome, as Stage 3 §5 already
defines it.

**The same two traps apply, and for the same reasons.**

- **Unannotated is not outside.** Where the proteome gives no topological domain,
  the check is `NOT_MEASURED`, not a pass and not a fail.
- **The lipid-anchored class reports zero annotated outward residues** because it
  has no transmembrane segment to annotate around. MSLN, CEACAM5 and PSCA are all
  in this class, and all three are known CAR targets. A residue-range test would
  reject every one of them. For this class the mature chain after signal peptide
  and propeptide removal is extracellular in its entirety, and the check is
  satisfied by the anchor, exactly as Stage 3's C5 special-cases them.

Handling this wrongly rejects three of the five known targets for this indication,
which is the loudest possible signal that it was handled wrongly.

### 4.2 Affinity is not to be maximised — and cannot be measured here anyway

**Higher affinity is not better for a CAR**, and a stage that sorted by it would be
optimising the wrong direction. High-affinity binders discriminate less well
between high and low antigen density, so they engage normal tissue expressing the
antigen at low level — which is precisely the risk Stage 3 spends its whole gate
trying to bound. Affinity that is too high also promotes antigen-independent tonic
signalling and exhaustion.

The right quantity is a **window**, not a maximum, and the window depends on the
antigen density Stage 4 measured. This stage therefore states the relationship and
does not score it.

**And it cannot be measured here at all.** Affinity is not in the PDB, and the
curated source that carries it is §2.2's unreachable one. Every candidate gets
`affinity: NOT_CONNECTED`.

Recording this rather than omitting it, because "we did not rank on affinity" and
"we could not rank on affinity" are different statements, and only the second is
true.

### 4.3 Format and the construct budget — where Stage 4's dual designs get priced

Stage 1 fixes the payload at 4.7 kb with 1.2 kb of backbone overhead: **3.5 kb of
construct budget** and at most 2 genetic edits. Stage 4 §7.4 explicitly declined to
size anything and named this stage as where it happens.

Approximate coding lengths, stated as the basis of the arithmetic rather than
hidden inside it:

| format | approx. bp | note |
| --- | --- | --- |
| `scFv` | ~750 | VH + linker + VL |
| `VH_VL` | ~700 | as separate chains |
| `VHH` | ~360 | single domain |
| `ligand` | varies | sized per candidate, never assumed |

Against 3.5 kb, with hinge, transmembrane, costimulatory and CD3ζ domains also in
the budget, **a dual-antigen AND-gated construct with two scFvs is tight and with
two VHH domains is comfortable.** That is a real design consequence of Stage 4's
`DUAL` output and it belongs in the same table as the recommendation.

The stage reports the sum and the headroom. It does not choose a format.

### 4.4 Species and immunogenicity — flagged, never filtered

Species of origin is in the PDB entity taxonomy and is reported. A murine binder
carries a real anti-CAR immunogenicity risk, and most historical CAR binders are
murine.

**Not a filter.** Immunogenicity assessment is Stage 9 with IEDB, per Stage 1.
Filtering here would import a Stage 9 judgement into a Stage 5 retrieval and would
do it with a cruder instrument. Flag, carry, and let the stage that owns the
question answer it.

### 4.5 The claim, stated so it can fail

For a target `T`, Stage 5's claim is:

```
at least one candidate exists                                  (1)
its epitope lies in T's extracellular region                   (2)
its format fits the remaining construct budget                 (3)
the antigen chain it was solved against is human               (4)
```

Each is a field in the output. A target failing (1) gets `NO_BINDER`; failing (2)
gets `EPITOPE_NOT_ACCESSIBLE`; failing (4) gets `NON_HUMAN_ANTIGEN`, which is a
warning rather than a rejection since a conserved epitope may still hold.

**`NO_BINDER` is not a verdict on the target.** It is a verdict on the literature.
Stage 4 §7.4's reasoning applies unchanged and is the reason §1 forbids re-ranking.

---

## 5. Retrieval

Per target, by symbol and by UniProt accession:

1. query the RCSB search API for entries mentioning the target
2. for each entry, read entity descriptions and taxonomy from the data API
3. classify entities as antibody chain, antigen chain, or other
4. retain entries where at least one antibody chain and the target antigen are
   both present
5. pull chain sequences from the FASTA endpoint
6. map the antigen construct's residue range onto the proteome record for §4.1

**Step 3 is the weak link and is declared as such.** Without SAbDab's curation,
antibody chains are identified from entity descriptions and sequence features
rather than from a curated annotation. That is a heuristic, it will have both false
positives and false negatives, and the count of entries where classification was
ambiguous **must be reported**. A retrieval that quietly resolved every ambiguity
would be presenting a guess as a lookup.

**Cache under the same discipline as every other source**: manifest as the commit
marker, fingerprint covering the query terms and the release, atomic write. The
PDB is versioned weekly; the pin is the query date and it goes in the manifest and
the configuration hash.

---

## 6. No single number

Stage 5 emits **no binder score**. Candidates are reported with their fields and
their verdict; they are not ranked into one ordering.

The reason is §2.3 and §4.2 together: the two quantities that would dominate any
sensible ranking — affinity and CDR composition — are exactly the two that are not
reachable. A score built from what remains would be built from structure count,
species and length, and would therefore rank targets by how much attention they
have had. That is B1, and the cleanest way not to trip it is to not build the
score.

Where an ordering is needed for display, it is the Stage 4 order, carried through
unchanged.

---

## 7. Rejection criteria — fixed in advance

Prefixed `B`. Stage 3's `R` and Stage 4's `P` criteria apply to the runs feeding
this one and are not re-checked here.

### Construction invariants

| id | invariant |
| --- | --- |
| J1 | every target in the pool appears in the output exactly once |
| J2 | no candidate is emitted whose antigen chain does not map to the target |
| J3 | `NOT_CONNECTED` fields are never compared, defaulted or sorted on |
| J4 | the pool order out equals the Stage 4 order in — this stage re-ranks nothing |

### Criteria

| id | criterion |
| --- | --- |
| **B1** | **candidate count correlates above 0.9 (Spearman) with the target's total PDB entry count** — the stage is measuring attention, not bindability |
| B2 | any of MSLN, CEACAM5 or PSCA is rejected by the ectodomain check — the lipid-anchor trap of §4.1 |
| B3 | a target with a known clinical CAR binder returns `NO_BINDER` |
| B4 | more than 30% of retained entries had ambiguous antibody-chain classification (§5 step 3) |
| B5 | any candidate carries a numeric affinity — the source is not connected, so a number here was invented |
| B6 | the stage emits a combined binder score of any kind |
| B7 | a `NO_BINDER` target is moved, dropped or reordered relative to its Stage 4 position |
| B8 | fewer than 3 of the 5 known targets for this indication return at least one candidate |
| B9 | the construct budget check passes for a two-scFv dual design without the arithmetic being shown |

**B1 detail — the criterion that says this stage is doing nothing.** The null here
is not subtle: retrieval from a structure database returns more for proteins that
have been studied more. §2.4 measured the gradient — 46 entries for CEACAM5 and 0
for NRG3 — and a stage whose output is rank-equivalent to that gradient has
re-derived the literature's priorities and called it a result. B1 is computed over
the whole pool, and its value is reported whether or not it trips.

The honest expectation is that B1 **will** be high, because bindability genuinely
does correlate with study. B1 tripping is therefore not automatically a defect in
the code — it is the signal that the output must be read as a literature survey
rather than as a property of the antigens, and that reading has to be printed
alongside it.

**B3 and B8 are the known-answer checks**, the analogue of the KRT19 control in
Stage 4. Retrieval that silently returns nothing looks identical to a target with
no binders, and only a target whose answer is known separates them.

### Explicitly not grounds for rejection

- a novel target returning `NO_BINDER` — that is §4.5 and it is expected
- the known targets dominating the candidate counts
- no VHH being available for any target
- a dual design not fitting the budget with two scFvs

---

## 8. Output

Header carries: the Stage 3 and Stage 4 configuration hashes verbatim and both
their criteria outcomes; the PDB query date as the release pin; **the SAbDab
status as `unreachable`, with §2.1's probe table**; the pool size and its
provenance; every threshold and coding length in §4.3; the ambiguous-classification
count from §5; B1's measured correlation; and Stage 5's own configuration hash
covering all of it.

Per target: Stage 4 outcome and partner, candidate count, best candidate's format
and size, ectodomain verdict, species, `affinity: NOT_CONNECTED`, and for `DUAL`
targets the two-binder budget arithmetic with headroom.

**No silent caps.** All retrieved candidates are written to file.

---

## 9. Open problems

**1. SAbDab is unreachable and two load-bearing fields go with it.** §2.2. The
options are an archived copy of the summary table, a numbering tool run locally to
recover CDRs, or accepting the gap. None is chosen here; the gap is declared.

**2. Affinity cannot be measured, and it is the field that would matter most.**
§4.2. Worse, the right target is a window whose position depends on antigen
density — so even with the data, this stage would need Stage 4's per-cell
measurements to place it, and those carry the dropout floor.

**3. This stage will find nothing for most of what Stage 4 recommends.** §2.4.
That is a true statement about the literature and a misleading one about the
biology, and the output has to carry both readings at once.

**4. Stage 1's dataset row needs splitting.** §2.5. One row cannot carry two
statuses, and the availability score is currently counting a source that is half
available as one blocking gap.

---

## Build note

Once approved, and with §0 understood:

1. `data/pdb.py` — search, entry metadata, FASTA, cached under the usual discipline
2. the known-answer check (B3, B8) run before anything else is written
3. `stages/stage5.py` — retrieval, then the ectodomain check, then sizing
4. `verify_binders.py` — four invariants, then nine criteria, then the biology

Same order and same rule as Stages 3 and 4.
