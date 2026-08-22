# Stage 5 — binder discovery

Written before `stages/stage5.py` exists and before any binder output exists.
Every source decision, threshold and rejection criterion below is fixed by this
document. If a criterion trips, the correction is a change to this document
followed by a re-run — never a narrative explaining why the output was acceptable
after all.

Status: **awaiting review. Nothing in `stages/stage5.py` may be written until this
is approved.**

Revision 2. Revision 1 was reviewed and rejected; §2 was rebuilt from a re-probe,
and §10 records what changed and why, because the reasons are themselves the
calibration this project runs on.

Section references of the form "Stage 4 §6" point at the stage's own spec in this
directory.

---

## 0. Preconditions

Two upstream stages are not in a state to be built on, and this one consumes both.
Both rows below were re-measured by running the verifiers, not read from a summary.

| stage | state | effect here |
| --- | --- | --- |
| Stage 3 | **R13 tripped** — clearance 0.6% protein-confirmed vs 42.8% RNA-supported, ratio 75.32x against a limit of 5x; 12 of 13 clear | the `cleared` flag is not yet meaning what it claims |
| Stage 4 | **5 of 16 criteria tripped** — P4, P7, P8, P12, P13; the run stops | the `DUAL` recommendations this stage would take as input are not yet a result |

Stage 4's measured outcome spread is `NO_DESIGN` 176, `UNRESOLVED` 10, `DUAL` 13,
`SINGLE` 1 — 200 decisions, one per pool member.

Stage 5 is specified now because the source question in §2 must be answered before
anything depends on it, and because writing the spec first is the rule. **It must
not be run for interpretation until both clear.**

P13 matters here specifically: one protein, **NRG3**, is the recommended partner
for 12 of Stage 4's 13 `DUAL` targets. §2.6 shows NRG3 has no structures under the
correct query either, so a binder stage run today reports a missing arm for almost
every dual design — which would look like a Stage 5 finding and is actually P13.

**A stronger statement than revision 1 made, and the measured one:** every one of
Stage 4's 14 recommendations has at least one arm with no structures. The single
`SINGLE` recommendation, NPSR1, is also zero. The problem is not dual-only.

---

## 1. Scope

Stage 5 answers one question per target: **is there a binding domain that can be
put in a CAR against this antigen, and what is known about it?**

A binder here is the antigen-recognition domain only — not the construct. Formats
are the `BinderFormat` enum in `schemas/project.py`: `scFv`, `VH_VL`, `VHH`,
`ligand`. **`Fab` is deliberately not added to it.** `BinderFormat` types
`ExistingBinder.format`, which is a *project input* field, so widening it would let
a caller supply a Fab as an existing binder and would change Stage 1's contract to
solve a Stage 5 reporting problem. Instead Stage 5 emits two fields: the deposited
form, from its own `DepositedFormat` vocabulary (`Fab`, `Fv`, `scFv`, `single
domain`, `Fab+Fc`), and the CAR-converted `BinderFormat`. §4.3 sizes the second and
reports the first.

**`DesignConstraints` does not constrain binder format today.** Stage 1's
`allowed_car_formats` is the *architecture* enum — `conventional`, `dual_target`,
`logic_gated`, `switchable`, `armored` — and naming it here would make an
implementer filter `scFv` out as "not admitted". If a binder-format constraint is
wanted it is a Stage 1 change, made there.

Stage 5 produces, per target:

- **retrieved candidates**, each with format, sequence, provenance route, the
  antigen entity it was solved or raised against, and that entity's accession
- **an ectodomain check** on the *epitope*, not the construct (§4.1)
- **a size in base pairs**, against the Stage 1 construct budget (§4.3)
- **a usability verdict per candidate**, an explicit **`NO_BINDER`** where nothing
  was found, per route with a target rollup, and a per-**design** verdict for
  `DUAL` targets (§4.5)

**Stage 5 does not design a binder, does not predict affinity, and does not dock
anything.** It retrieves and characterises what exists. Anything else would be a
prediction presented beside measurements, and the two would be read alike.

**Stage 5 does not re-rank targets.** Stage 4 §7.4 flags binder availability
rather than filtering on it, precisely so that the target list is not shaped by
which proteins happen to have been crystallised. Stage 5 must not undo that by
reordering anything: it annotates, and the annotation travels.

---

## 2. Sources — measured, not assumed

Stage 1 declared one blocking dataset for this stage,
`("SAbDab and PDB", "binder retrieval", [5], True)` — the only required dataset
the availability score had never satisfied. **That row has since been split**, per
§2.5, into `("PDB", "binder structures", [5], True)` and
`("SAbDab", "antibody chain and numbering annotation", [5], True)`. Neither has a
connector yet, so both read `not_configured` and the availability score is 0.750
(6 of 8 blocking) rather than 0.857. Nothing became less available: counting two
unconnected sources as one had understated the gap.

**Revision 1 concluded that half of it was unreachable. That conclusion was
wrong**, and the way it was wrong is the thing worth recording: the probe saw an
HTML body, correctly refused to treat it as data, and stopped — without asking
what the HTML was *for*. The guard against "a 200 with an HTML body is not a
source" is right. Its mirror image is also true and was missing: **an HTML body is
not proof of absence.**

### 2.1 The probe table is not a literal

**A probe table transcribed into a spec is a probe table that cannot be wrong at
run time.** Revision 1's table was produced once, with a client whose CA bundle
dates from 2017 and which returns false certificate errors and false `000`
statuses; one row recorded a local TLS artefact as a property of a source, and the
SAbDab rows recorded a `200` at a path that returns `404`.

So §8 no longer requires a copied table. It requires **a probe performed at run
time**, recording per endpoint:

`url`, `http_status`, `content_type`, `bytes`, `sha256_prefix`, `body_shape`, and
the **HTTP client and TLS stack that made the probe**.

`body_shape` comes from a fixed vocabulary: `json`, `csv`, `tsv`, `fasta`,
`mmcif`, `server_rendered_html`, `js_app_shell`. `js_app_shell` is detected
structurally — a small HTML body with an empty root element and a module script —
and **may never stand as a terminal verdict without a recorded follow-up probe of
the shell's own API.** That rule is written because its absence is what produced
revision 1.

The client field is recorded because this project now has a documented instance of
a stale client fabricating a source verdict.

### 2.2 What is actually reachable

Re-probed through a working TLS stack. Counts are given next to the count the
source itself reports.

| source | endpoint | status / type / bytes | verdict |
| --- | --- | --- | --- |
| Structure entries | `data.rcsb.org/rest/v1/core/entry/...` | 200 · json | **usable** |
| Structure search | `search.rcsb.org/rcsbsearch/v2/query` | 200 · json | **usable** |
| Structure sequence | `rcsb.org/fasta/entry/...` | 200 · fasta | **usable** |
| Structure coordinates | `files.rcsb.org/download/*.cif` | 200 · mmcif | **usable** |
| Antibody structure summary | `…/newsabdab/api/download/all-summary` | 200 · csv · 11,694,156 B | **usable** — 21,914 rows x 45 cols |
| Antibody numbering + CDRs | `…/newsabdab/api/antibody-instances` | 200 · json | **usable** — pre-computed CDRs |
| Antibody DB release pin | `…/newsabdab/api/stats/summary` | 200 · json · 303 B | **usable** |
| Curated chain annotation | `…/newsabdab/api/rcsb-pdb-annotations` | 200 · json | **usable** |
| Therapeutic antibodies | `…/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv` | 200 · csv · 636,357 B | **usable** — 1,133 rows x 24 cols |
| Patent/literature antibodies | `…/plabdab/static/downloads/paired_sequences.csv.gz` | 200 · csv.gz · 11,392,439 B | **usable** |
| Bioactivity / mechanism | `ebi.ac.uk/chembl/api/data/...` | 200 · json | **usable** for mechanism; **empty** for affinity (§4.2) |
| Immunogenetics portal | `imgt.org` | 200 · html · 1,153,574 B | **reachable, server-rendered; no JSON route exercised — `NOT_VERIFIED` as a machine-readable source** |

The web front end of the antibody structure database serves a JavaScript
application shell (1,457 bytes) at its summary path. That observation is
reproducible and was correctly made in revision 1. The **verdict** drawn from it
was wrong: the shell's module declares a REST API, documented at
`…/newsabdab/api/openapi.json` (200, json, 27,159 B, OpenAPI 3.1.0, CC BY 4.0).

Cross-checks that must hold, and did: 21,914 summary rows against 21,914 reported
antibody instances; 11,458 distinct structure entries against 11,458 reported;
6,780 distinct antibody identifiers against 6,780 reported.

**Release pin** for the manifest: the `last_update` field from the stats endpoint
plus the API `info.version`, alongside the structure-database query date.

### 2.3 What is recoverable and what is not

Revision 1's table said four things were lost. Three of them are not.

| field | revision 1 said | measured |
| --- | --- | --- |
| antibody heavy/light chain identity | recoverable, imperfectly, from text | **curated**, from the chain-annotation endpoint, keyed to entity and asym ids |
| which chain is the antigen | recoverable from entity descriptions | **curated**, same endpoint |
| CDR boundaries under a numbering scheme | **not recoverable without a numbering tool** | **served pre-computed**, IMGT, per variable segment |
| species of origin | in the entity taxonomy | in the entity taxonomy — but see §4.4, it does not mean what §4.4 used it for |
| therapeutic name and clinical stage | **not recoverable** | **one 636 KB CSV**: name, format, isotype, highest clinical trial, status, target, companies, conditions, alternative names |
| **VH/VL sequences of clinical binders** | not contemplated | **in the same CSV** |
| affinity where curated | not in the structure DB | **still not available anywhere reachable — §4.2** |

**The two that survive as real gaps are affinity and humanisation state.** Both are
declared, and neither is imputed.

The therapeutic route also joins cleanly to the structure route: the heavy-chain
sequence for a named therapeutic is byte-identical to the sequence of its deposited
variable segment, so the two sources can be joined on sequence rather than on name.

### 2.4 What this changes about the stage

**A candidate no longer has to be a solved structure.** That is the substantive
consequence, and it changes what the stage is for. Two retrieval routes, reported
separately and never merged into one count:

- **structure route** — a deposited complex. Carries an epitope, so §4.1's check
  can run. This is the only route that can answer "where does it bind".
- **sequence route** — a named therapeutic or a patent/literature antibody with
  VH/VL and a target annotation. Carries no epitope, so §4.1 returns
  `NOT_MEASURED` for it, and it carries clinical stage, which §6 and B5'' forbid
  sorting on.

**Verdicts are per route, with a target-level rollup**, because the two routes can
disagree and a single verdict cannot carry that. CEACAM6 is the worked case: 6
structure entries, no antibody complex, one named therapeutic. Its structure-route
verdict is `NO_BINDER`; its sequence-route verdict is a candidate; its target
rollup is `BINDER_SEQUENCE_ONLY`. A target with a sequence-route binder and no
structure is therefore **not** `NO_BINDER` at target level, and the distinction is
load-bearing: the VH/VL sequence is the thing a CAR is actually built from, and
`NO_BINDER` would discard it.

`NO_BINDER` at target level means both routes are empty.

### 2.5 The status Stage 1 must carry

Revision 1 asserted that the dataset row is already `unreachable`. It is
`not_configured`, and `resolve_status` derives `UNREACHABLE` from
`source.is_cached()` — from "a connector exists but nothing is cached" — not from
reachability at all.

**Stage 1's dataset list must be amended** to split the row into a structure
source and an antibody-annotation source, each with independent status. That
touches `stage1._DOWNSTREAM_DATASETS` **and** `data/availability.py:CONNECTORS`,
and `availability._ORPHANED` raises at import time if a connector is registered
under a name Stage 1 does not emit.

**The split is one atomic change, not a sequence.** The guard is one-directional:
connector-without-row raises at import and is caught for you; row-without-connector
raises nothing and silently changes the availability denominator. There is no safe
intermediate state to stop in, so neither ordering is prescribed — both files move
together.

**Two hardcoded expectations move with them**, and neither is in the files above:
`verify_schema.py` asserts fixed dataset and blocking counts in **both** discovery
and validation mode — four literals, not two — and `availability.py` carries the expected score `0.857` as a literal. The
row split breaks both, and a spec that named only the two obvious files would have
left a verifier failing for a reason that looks unrelated.

State the resulting availability score in the same change. The prior framing — six
of seven blocking sources connected, this being the seventh — stops being true the
moment the row splits, and the score must be restated rather than left to be
re-derived from a sentence that no longer exists.

### 2.6 Structure counts, measured with the right instrument

Revision 1 counted **full-text hits on a name string**. Those counts reproduce
exactly — and the quantity is wrong. Full text counts documents matching a string,
not structures of a protein. Measured, same session: `AMN`, an actual Stage 4
`DUAL` target, returns 369 full-text hits whose top entries include a bacterial
RNA chaperone, a sulfur transferase and a photosystem supercomplex. Not one is
amnionless. The table also mixed conventions, using gene symbols for most rows and
protein names for two — and the two differ by up to 19x on the same API.

The correct instrument is an accession-anchored query: `exact_match` on
`rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession`
with `database_name` UniProt. Counts under it:

| target | accession | entries | | target | accession | entries |
| --- | --- | --- | --- | --- | --- | --- |
| MUC1 | P15941 | 23 | | EPCAM | P16422 | 2 |
| MSLN | Q13421 | 12 | | PSCA | O43653 | 1 |
| CEACAM5 | P06731 | 6 | | **NRG3** | P56975 | **0** |
| CEACAM6 | P40199 | 6 | | NPSR1 | Q6W5P4 | 0 |
| EFNA5 | P52803 | 6 | | TMC5 | Q6UXY8 | 0 |
| CLDN18 | P56856 | 3 | | VMP1 | Q96GC9 | 0 |
| PTPRN2 | Q92932 | 3 | | | | |

The error is not a uniform inflation, and saying so would be its own overstatement.
It runs from 7.7x too high (CEACAM5, 46 against 6) through 6.5x (EPCAM), 5.7x
(CLDN18), 1.4x (MUC1) and 1.25x (MSLN) to exactly right (CEACAM6, 6 against 6) and
one that was too *low* (EFNA5, 5 against 6). The distortion is target-specific,
which is worse than a scale error: it reorders. CEACAM5 falls from first to joint
third and MUC1 rises to first.

**The qualitative conclusion survives and the gradient does not.** "46 for CEACAM5
and 0 for NRG3" was an order of magnitude overstated; the honest figures are 6 and
0. NRG3 is genuinely 0 under the correct instrument, so §0's reasoning holds — but
it held for a reason that would not have survived re-derivation, which is exactly
the kind of thing this project writes down.

**Structure count measures how much attention a protein has had, not how good a
CAR target it is.** Any Stage 5 that scores targets by what it can retrieve will
rank the known targets top and the novel ones bottom whatever the biology says.
That is the null for this stage and it is criterion **B1**.

### 2.7 The zero-result contract

A zero-hit structure search returns **HTTP 204 with a zero-byte body**, not 200
with a count of zero. `json.loads(b"")` raises, and the natural repair is an
exception handler that yields "0 candidates" — after which *"this protein has no
structures"* and *"the query broke"* are byte-identical outcomes.

Every zero in §2.6 is a 204. NRG3 = 0 is the single number §0 rests on, and it is
exactly the number a silent failure would fabricate. The contract is fixed here and
enforced by **J5** and **B12**.

---

## 3. The pool

**The 200 `Decision` records Stage 4 emits, one per pool member, in all four
outcome classes.** `decide()` iterates the pool and appends exactly one decision
per member; revision 1 described this single set as two ("plus every target in the
Stage 4 pool"), which reads as a union and invites a pool larger than 200.

The pool is Stage 3's top 200 by composite with risk ignored entirely,
deduplicated to one accession per gene.

Not just the `DUAL` and `SINGLE` recommendations. A `NO_DESIGN` target failed on
risk or coverage, not on bindability, and re-running retrieval is expensive.

Binder retrieval is per target and independent per target, so there is no pairing
combinatorics here: 200 retrievals, not 19,900.

### 3.1 Stage 4 must persist its decisions

**Stage 4 writes nothing to disk.** Neither `stage3.py` nor `stage4.py` emits an
artifact; `Decision` carries no accession; and `verify_pairing.main()` returns
before the decisions are ever printed when a criterion trips — which it does. There
is nothing for Stage 5 to read.

Stage 4 must gain a writer emitting, per decision: `gene`, `accession`, `outcome`,
`partner`, `partner_accession`, `pool_index` — under the same cache discipline as
every source. Without it Stage 5 must re-run Stages 3 and 4 in process to re-derive
the very numbers §0 says are unsettled, and **J4** is unverifiable against anything
on disk. This is a Stage 4 change and it is a precondition, not a nicety.

---

## 4. What makes a binder usable

### 4.1 The epitope must be on the part the CAR can reach

A CAR binds a cell surface from outside. A structure showing an antibody bound to
its antigen proves binding; it does not prove the epitope is **accessible on an
intact cell**.

**The check is on the epitope, not on the construct.** Revision 1 tested whether
the antigen construct's residue range fell in an extracellular topological domain.
A multi-pass protein can only be crystallised whole, so that test rejects
CLDN18 — one of the two validation cases, and an approved drug's target — because
the construct spans 1–208 against an extracellular annotation covering 84 of 261
residues. The test must be: **the antigen residues within contact distance of the
antibody chains, computed from the deposited coordinates, must fall in an
extracellular topological domain.** Where contacts cannot be computed, the verdict
is `NOT_MEASURED`.

**Contact is fixed here as any heavy atom of the antigen within 4.5 Å of any heavy
atom of an antibody chain**, with at least 3 such residues required before any
verdict other than `NOT_MEASURED` is emitted. Both numbers are fixed before output
exists, are printed in §8's header, and are not to be revisited once candidate
lists are visible — a cutoff chosen after seeing which targets pass is exactly the
tuning this project forbids, and B2, B11 and B18 all resolve through it.

**The proteome record does not currently carry the ranges this needs.**
`ProteinRecord` exposes `extracellular_residues: int | None` — a scalar produced by
`_count_extracellular_residues`, which parses `lo` and `hi` from each topological
segment and then discards them. A sum cannot answer a membership question. The
parser already computes the ranges; `extracellular_ranges: list[tuple[int, int]]`
must be retained on the record. That much is free: `ft_topo_dom` is already in
`uniprot.FIELDS`.

**The rest is not free, and the cost is a proteome re-fetch.** The chain-selection
rule below needs `ft_chain`, which is not in
`FIELDS`, and `FIELDS` is inside the cache fingerprint — so adding them invalidates
the proteome cache and re-fetches all 20,431 entries.

**Isoform sequences are not a column at all**, and this is the sharper problem.
They are not obtainable by adding a field; they require requesting isoforms in the
query, which **changes the entry set itself** — so the 20,431 reviewed entries and
the 3,466 surface proteins derived from them stop being the same numbers, and every
count this project reports against them moves. Isoform sequences must therefore be
fetched as a **separate, additively-cached entry keyed on the accessions that need
them**, not by widening the proteome query. Which accessions qualify is **not** circular: the count of isoforms per entry and
their differing ranges are metadata on the canonical record and can be enumerated
without pulling any isoform sequence. The enumeration query is run first, its result
count is reported, and only the qualifying accessions have sequences fetched. Fetching
for those is cheap and leaves the reference counts intact.
The alternative — deriving chain boundaries from anything other than the
annotation — is the imputation this project forbids.

Revision 1 cited "Stage 3 §5" for this check. Stage 3 §5 defines a **saturating
size term for a score**, not a membership test, and `_score_c5` returns 1.0 or 0.6
— a pass/fail check has no analogue of 0.6. That citation is withdrawn.

**Seven traps, each of which silently returns a wrong answer.**

- **Unannotated is not outside.** Where the proteome gives no topological domain,
  the check is `NOT_MEASURED`, not a pass and not a fail.

- **The lipid-anchored class reports `None`, never zero.**
  `_count_extracellular_residues` returns `total if measured else None`, and
  `measured` is only set inside a branch that has already added a positive length —
  the function is structurally incapable of returning 0. Stage 3 keys its own
  special case on `None`. An implementer writing revision 1's wording (`== 0`) gets
  a branch that never fires, and the anchored targets fall through to
  `NOT_MEASURED`, which §4.1 itself defines as neither pass nor fail. The trap is
  sprung in a form the criterion cannot see.

- **The anchored class needs a pass rule, and revision 1's was wrong rather than
  unnecessary.** A GPI-anchored protein has no transmembrane segment to annotate
  around, so it carries no extracellular topological domain and the check above
  returns `NOT_MEASURED` for every one of them — which B2 must treat as a failure.
  The replacement rule: **for a GPI-anchored protein, the chain carrying the
  lipidation site is extracellular in its entirety, and epitope contacts falling
  inside that chain satisfy the check.** That is narrower than revision 1's "the
  mature chain after signal peptide and propeptide removal", which is what admitted
  mesothelin's secreted chain below.

- **The class is computed, not listed.** Revision 1 named MSLN, CEACAM5 and PSCA.
  The predicate already exists: `uniprot.py` sets `gpi` from the lipidation
  annotation and assigns `MembraneClass.GPI_ANCHORED`. Drive the special case off
  that. Measured GPI status among the five known targets: MSLN, CEACAM5 and
  CEACAM6 are anchored; MUC1 and CLDN18 are not. **PSCA is anchored but is not one
  of the five**, and CEACAM6 — which is — was omitted. EFNA5 is also anchored and
  is in the pool. Note also that the predicate matches a GPI anchor only:
  myristoyl, palmitoyl and prenyl anchors are not captured, so "the lipid-anchored
  class" is broader in prose than in code. Say GPI.

- **Mesothelin's mature chain is not extracellular in its entirety.** Revision 1's
  rule — "the mature chain after signal peptide and propeptide removal" — admits
  residues 37–598. Residues 37–286 are megakaryocyte-potentiating factor, which the
  proteome annotates as **secreted**: a soluble serum protein, not on the cell. 250
  of 562 residues, 45%, are declared surface-accessible and are not. A binder
  solved against that region would pass and be reported as a usable CAR binder for
  a target it can only meet in plasma.

  **The rule is the smallest chain containing the anchor position**, and the
  qualifier is not pedantry — it is what the annotation forced. Now that chain
  boundaries are fetched, MSLN measures as **three** chains, not two:

  ```
   37..598   Mesothelin                          spans the whole precursor
   37..286   Megakaryocyte-potentiating factor   secreted
  296..598   Mesothelin, cleaved form            anchored
  ```

  Two of the three contain the anchor, so "the chain carrying the anchor" is
  ambiguous and a first-match or longest-match reading returns 37–598 — the whole
  precursor, which is the original defect restored. Taking the smallest chain that
  contains the anchor gives 296–598. Emit the chosen chain id so the choice is
  auditable, and report the count of proteins where more than one chain contained
  the anchor, because that count is how often this rule was load-bearing rather
  than incidental. Measured: 190 of the 3,466 surface proteins carry more than one
  chain. Any candidate whose antigen range overlaps a chain annotated secreted is
  flagged `SHED_ANTIGEN`, never passed.

- **The deposited numbering is not the proteome's numbering.** Measured on the
  amatuximab complex 4F3F: the antigen entity's own sequence position 2 aligns to
  proteome position **302**, length 58 — the epitope is precursor residues
  302–359, and the entry is titled "Msln7-64" in mature-chain coordinates. An
  implementer taking the deposited residue numbers, or parsing the title, places
  the epitope at 7–64 — inside the signal peptide — and gets a confident wrong
  answer with nothing raised. **The mapping comes from the entity's alignment
  block (`rcsb_polymer_entity_align.aligned_regions`, entity offset to reference
  offset), never from deposited numbering and never from a description string.**
  302–359 falls inside the anchored chain 296–598, so the chain-selection rule
  above and this one agree, which is the check that they are both right.

  The anchor position this needs is also not on the record: `parse_row` reduces
  `ft_lipid` to a boolean and keeps no coordinate, so `lipid_anchor_position: int |
  None` must be retained alongside `extracellular_ranges`. Same parser, same field
  already in `FIELDS`, same reason — the value is computed and then discarded.

  And `parse_row` unpacks a **fixed seven columns**, padding and truncating to that
  width. Adding `ft_chain` without widening the unpack silently drops it, and the
  drop is only visible after the 20,431-entry re-fetch has already been paid for.
  The unpack width and `FIELDS` must change in the same edit.

- **Tandem repeats make coordinates meaningless.** MUC1's extracellular domain is
  24–1158 and its epitope peptide occurs roughly 40 times, so two structures of the
  identical epitope receive residue ranges hundreds of residues apart and any
  overlap, dedup or per-residue mapping is noise. Where the antigen construct
  sequence occurs more than once in the target, emit
  `REPEAT_COORDINATE_AMBIGUOUS` with the occurrence count instead of a range.

Handling the anchored class wrongly rejects three of the five known targets for
this indication, which is the loudest possible signal that it was handled wrongly.

### 4.1.1 Isoform is part of the target's identity

**CLDN18 is not one protein for this purpose, and nothing in the pipeline knows
that.** Every source in Stages 2–4 is keyed on gene symbol; the word "isoform"
appears nowhere in the codebase or in any prior spec. CLDN18 has two mutually
exclusive first exons: **CLDN18.1 is the lung isoform and CLDN18.2 is the gastric
one and the therapeutic target.** The two are the same length with identical
topology annotation, so every residue-range check succeeds numerically while being
computed against the wrong sequence, with no error raised anywhere. The differences
concentrate in the epitope: 8 of the 21 differing residues sit in the first
extracellular loop, which is where the approved binder binds.

**And it is not confined to CLDN18.** Measured: the antigen in 7UED, one of the
pinned mesothelin entries, is deposited as "**Isoform 4** of Mesothelin". So the
second validation case carries an isoform assignment too, stated in the entity
description and absent from the accession. A stage keyed on accession alone
records both as `Q13421` and cannot tell which sequence the epitope was mapped
against.

This is not confined to Stage 5. Measured on the current pipeline, gene-level
CLDN18 carries 150.96 TPM in normal lung, and `_baseline_score(150.96)` is 0.7271
— **exactly CLDN18's risk, whose peak organ is therefore the lung, on the wrong
isoform.** The AND-gate consequence is direct:

| MSLN + CLDN18, conservative gate | risk | peak organ |
| --- | --- | --- |
| as measured today, gene-level | 0.6366 | lung |
| with CLDN18 lung measured and absent, i.e. isoform-resolved | 0.2277 | gi_tract |

Isoform resolution moves the number 2.8x and relocates the peak organ from a
tissue the target isoform is not in to the stomach — which is the approved
antibody's actual dose-limiting toxicity. **It does not rescue the pair**; both
figures are above the 0.15 ceiling. It is a correction, not a fix, and it belongs
to Stage 3 as much as to Stage 5.

For Stage 5 the rule is: **isoform is an explicit key on the target record.** Any
target whose proteome entry has more than one isoform differing inside an
extracellular topological domain must carry an explicit isoform assignment or emit
`ISOFORM_UNRESOLVED`, and a candidate may not be reported usable with an unresolved
isoform. Assignment is by aligning the deposited antigen sequence to each isoform
— not by accession, which does not carry it, and not by parsing a description
string. Enforced by **B11**.

### 4.2 Affinity is not to be maximised — and is not obtainable

**Higher affinity is not better for a CAR**, and a stage that sorted by it would be
optimising the wrong direction. High-affinity binders discriminate less well
between high and low antigen density, so they engage normal tissue expressing the
antigen at low level — precisely the risk Stage 3 spends its whole gate bounding.

The right quantity is a **window**, not a maximum. This stage states the
relationship and does not score it.

Two corrections to revision 1, both of which would have propagated into a wrong
action:

- **Tonic signalling is misattributed.** Antigen-independent signalling does not
  involve antigen binding, so the antigen dissociation constant cannot drive it.
  It is driven by binder framework self-association and aggregation propensity, and
  by hinge and transmembrane choice. Acting on revision 1's wording would
  deprioritise a high-affinity non-aggregating binder while passing a low-affinity
  aggregation-prone one. That matters here specifically, because the candidates are
  murine Fabs being converted to scFv, where aggregation propensity is the
  unmeasured hazard.

- **The window cannot be placed from anything this pipeline has.** Revision 1 said
  it "depends on the antigen density Stage 4 measured". Stage 4 measured a per-cell
  *detection fraction* from single-cell RNA — a dropout-limited probability of a
  non-zero count, not surface protein copies per cell. Placing a dissociation
  constant window from it is a unit error. **Surface copy number per cell is not
  available anywhere in this pipeline**, and that is stated rather than worked
  around.

**And affinity itself is not obtainable — for a different reason than revision 1
gave.** Revision 1 said the curated source carrying it is unreachable. The source
is reachable and **no longer carries the field**: affinity, free-energy, method and
temperature columns are absent from the current release's 45-column summary and
from its API schema. Independently, the bioactivity route returns
`total_count = 0` for both validation molecules and both validation targets.

Every candidate gets `affinity: NOT_CONNECTED`. Recording the corrected reason
matters, because with revision 1's reason left standing the next reader would read
"the source is reachable now" as "affinity is available", switch the field on, and
populate it from a source that does not have it — or from a prediction.

The retired release did carry those columns. Recovering them from an archive is
**not assumed to be possible**: no live URL for that file was reproduced.

### 4.3 Format and the construct budget — with the arithmetic shown

Stage 1 fixes the payload at 4.7 kb with 1.2 kb of backbone overhead:
`max_construct_kb = 3.5`, `max_genetic_edits = 2`, and — the term revision 1 omitted
entirely — **`require_safety_switch`**.

That last is **not a constant**: `TOLERANCE_RULES` derives it from the project's
safety tolerance, and it is `True` only for `CONSERVATIVE`. This project is
conservative, so it is `True` here and every sum below includes the switch. A
moderate or permissive project drops 1,308 bp and the conclusion inverts — two
scFvs then fit with 929 bp spare. **The sums are therefore conditional on the
tolerance, and the stage must read the flag rather than assume it**, or B9 trips a
correct implementation on the first non-conservative project.

Coding lengths are fixed here, in advance, so the sum is checkable. Base pairs are
three times the residue count; the residue counts are the stated basis of the
arithmetic rather than hidden inside it.

| component | aa | bp |
| --- | --- | --- |
| signal peptide | 21 | 63 |
| `scFv` (VH + linker + VL) | 250 | 750 |
| `VH_VL` as separate chains | 233 | 699 + 129 |
| `VHH` single domain | 120 | 360 |
| `Fab` as deposited | 469 | 1,407 |
| `ligand` | varies | sized per candidate, never assumed |
| hinge | 45 | 135 |
| transmembrane | 24 | 72 |
| costimulatory domain | 42 | 126 |
| CD3ζ | 112 | 336 |
| ribosomal skip peptide | 22 | 66 |
| safety switch | 414 | 1,242 |
| stop codon | — | 3 |

`VH_VL` carries a `+ 129` because two separate chains need a second signal peptide
(63) and a second skip peptide (66) that the single-chain formats do not. It is
written into the row rather than left to the sums below, which cost only the
single-chain formats — a design costed at a flat 699 would be short by exactly that
amount, per arm, and B9 would not see it.

**An AND gate is two receptors, not two binders.** Pricing it as two binders is
what let revision 1 assert a conclusion it could not derive. Split-signal
architecture, both receptors plus the mandatory switch:

```
receptor 1   signal 63 + binder + hinge 135 + TM 72 + CD3ζ 336
receptor 2   signal 63 + binder + hinge 135 + TM 72 + costim 126
             + skip peptide 66 x2 + stop 3 + safety switch 1,242
```

| design | binder cost | construct total | vs 3,500 |
| --- | --- | --- | --- |
| two `scFv` | 1,500 | **3,879** | **over by 379 — does not fit** |
| two `VHH` | 720 | **3,099** | fits, 401 spare |
| two `Fab`, unconverted | 2,814 | 5,193 | far over — see below |

**The single-antigen case is costed too**, because 187 of the 200 pool members are
not dual and `BUDGET_EXCEEDED` and B9 apply to every one of them. One receptor plus
the mandatory switch:

```
signal 63 + binder + hinge 135 + TM 72 + costim 126 + CD3ζ 336
       + skip peptide 66 + safety switch 1,242 + stop 3
```

| design | binder cost | construct total | vs 3,500 |
| --- | --- | --- | --- |
| one `scFv` | 750 | 2,793 | fits, 707 spare |
| one `VHH` | 360 | 2,403 | fits, 1,097 spare |
| one `Fab`, unconverted | 1,407 | 3,450 | fits, 50 spare — which is why §4.3 sizes from variable domains and not from deposited chains |

**So revision 1's headline was wrong.** "Two scFvs is tight and two VHH domains is
comfortable" becomes, with the mandatory switch priced: two scFvs **do not fit**,
and two VHH domains do, with 401 bp of headroom. The conclusion is now a computed
result rather than an assertion, which is what **B9** exists to require.

Two further statements the arithmetic forces:

- **The deposited form is a Fab and the budget is not.** Every retrieved antibody
  in the validation set is a Fab. Sizing a candidate from its deposited chain
  length gives 1,407 bp where the CAR domain is 750 — a false "does not fit".
  Candidates are sized from **variable-domain boundaries**, with expression tags
  stripped, and the deposited format and the CAR-converted format are reported as
  **two separate fields** so the conversion is visible.
- **A split-signal design is not a true AND gate.** The CD3ζ-only receptor drives
  partial activation on antigen A alone, so the safety property Stage 4's `DUAL`
  recommendation was purchased for is not fully delivered by the construct that
  was priced. A synNotch AND gate does not fit the budget at all. That is recorded
  here rather than discovered later.

Whether `max_genetic_edits = 2` binds a dual design — one bicistronic cassette
versus two constructs — is stated as out of scope for this stage and left to the
construct stage, rather than mentioned once and never used.

The stage reports the enumerated components, the sum and the headroom. It does not
choose a format.

### 4.4 Species — two fields, and not a humanisation state

The entity taxonomy reports **the source organism of the deposited construct**. It
is not a humanisation state and cannot distinguish murine from chimeric from
humanised. The approved CLDN18.2 antibody is chimeric and its variable chains are
deposited as mouse; revision 1's rule would flag it as a murine binder carrying "a
real anti-CAR immunogenicity risk", which is the opposite of the clinical record
for that molecule.

Three separate fields, because revision 1 collapsed them into one and a reader
could not tell which molecule a value referred to:

- `antigen_species` — from the **UniProt cross-reference** of the mapped accession,
  never from the deposited source organism, which describes the expression system
- `binder_source_organism` — the deposited construct's source organism
- `humanisation_state` — `NOT_CONNECTED`

For the approved CLDN18.2 complex the honest answer is both, differently: antigen
human, binder chains mouse.

**Not a filter.** Immunogenicity assessment is Stage 9, per Stage 1. Filtering here
would import a Stage 9 judgement into a Stage 5 retrieval with a cruder instrument.

### 4.5 The claim, stated so it can fail

For a target `T`, Stage 5's claim is:

```
at least one candidate exists                                  (1)
its epitope lies in T's extracellular region                   (2)
it fits the remaining construct budget                         (3)
the antigen it was raised against is human                     (4)
the isoform it binds is T's isoform                            (5)
```

Each is a field. Failing (1) gives `NO_BINDER`; (2) `EPITOPE_NOT_ACCESSIBLE`;
**(3) `BUDGET_EXCEEDED`**; (4) `NON_HUMAN_ANTIGEN`, a warning rather than a
rejection since a conserved epitope may still hold; (5) `ISOFORM_UNRESOLVED`.

Revision 1 gave clause (3) no verdict at all and then listed "a dual design not
fitting the budget" under *not grounds for rejection* — so a claim "stated so it
can fail" was stated so it could not. A `DUAL` target whose two binders do not fit
received the same verdict as one whose binders do. `BUDGET_EXCEEDED` fixes that.
Given §4.3, it will fire on every two-scFv dual design, which is the point.

**A design verdict, distinct from the target verdict.** §1's "one question per
target" cannot describe a dual design, which has two arms. Measured, 13 of 13 dual
designs have at least one arm with zero structures, and for 12 of them the empty
arm is the partner. Several dual targets do have entries on the target arm — the
per-target counts across all 13 are measured in the run, not asserted here, and
§2.6 deliberately does not carry them because entries are not candidates: CEACAM6
and PSCA both have entries and are pinned at zero candidates. Wherever the target
arm does yield a candidate, the row would otherwise show it beside a headroom
figure computed from one binder, reading as a feasible dual design when the partner
arm has nothing.

So: `DESIGN_NO_BINDER` where either arm has no candidate, and the budget row is
emitted as `NOT_COMPUTABLE` naming the missing arm — **never as a headroom number
derived from one binder.**

**Keyed on `outcome == DUAL`, which is the design's arity — not on `partner is not
None`, which is not.** `decide()` sets a partner on the `SINGLE` branch too
whenever an admissible pair exists, so a non-null partner there is an artefact of
how the decision was reached and not a second arm. Costing it as one would price a
two-binder construct for a single-antigen design.

The artefact must still be visible rather than dropped: a `SINGLE` carrying a
non-null partner is **reported with that partner and explicitly not costed**. It
does not arise in this run only because the one `SINGLE`, NPSR1, happens to carry
`None` — which is luck, not a guarantee, and is why the rule is written down.

**`NO_BINDER` is not a verdict on the target.** It is a verdict on the literature.
Stage 4 §7.4's reasoning applies unchanged and is the reason §1 forbids re-ranking.

---

## 5. Retrieval

Per target, keyed on **UniProt accession**, with symbol full text used only as a
recall supplement whose extra hits are counted separately and must each pass the
accession test before becoming candidates.

**Structure route**

1. query the structure search API with `exact_match` on the accession attribute
   named in §2.6
2. read entity descriptions, taxonomy and UniProt cross-references from the data API
3. classify chains from the **curated chain-annotation endpoint** — not from
   description text
4. retain entries where an antibody chain and the target antigen are both present
5. pull chain sequences, and verify each against the entity record's own length and
   leading residues (**B15**)
6. compute antibody–antigen contact residues from the coordinates, translate them
   through the entity alignment block into proteome coordinates, and map them onto
   the proteome record for §4.1; assign the isoform by aligning the deposited antigen sequence to each isoform
   sequence, per §4.1.1 — an isoform named in the entity description (7UED reads
   "Isoform 4 of Mesothelin") is corroboration and a mismatch trips, but it is
   never the assignment itself

**Sequence route**

1. match the therapeutic table on its target field by **tokenising first, then
   matching exactly**. The field is compound and synonym-laden — measured values
   include `CEACAM5/CD66e`, `MUC1/PEM/EMA`, `CLDN18;CD3E`, `IAP/CD47;CLDN18` — so
   the rule is: split on `;` to separate the antigens of a bispecific, split each
   on `/` to separate synonyms, strip whitespace, then require an exact match on a
   resulting token. Neither raw substring nor exact match on the whole field works:
   the whole field misses `CEACAM5/CD66e` entirely, and substring pulls MUC16 and
   MUC18 into the MUC1 bucket **and matches `CLDN1` inside `CLDN18`**, which the
   table contains as separate targets. Under this rule the pinned counts hold:
   CLDN18 11, CEACAM5 5, MSLN 4, MUC1 3, CEACAM6 1.
2. take VH/VL, format, isotype, highest clinical trial, status, companies,
   conditions, alternative names
3. join to the structure route on sequence identity where a structure exists
4. supplement from the patent/literature table for targets with no named therapeutic

**Step 3 of the structure route is no longer the weak link.** Revision 1 declared
chain classification a heuristic over description text and required its ambiguity
rate to be reported. The curated annotation supplies heavy, light and antigen chain
identity keyed to entity and asym ids, so the heuristic is replaced rather than
measured. Where an entry is absent from the curated annotation the heuristic
applies, and **those** entries are counted and reported (**B4**).

**Cache under the same discipline as every other source**: manifest as the commit
marker, fingerprint covering query terms and release pins, atomic write. A cached
empty must record the status it was born from (§2.7).

---

## 6. No single number

Stage 5 emits **no binder score**.

The reasoning has changed and the conclusion has not. Revision 1 justified it by
unreachability: the two quantities that would dominate a ranking were the two that
could not be fetched. One of those — CDR composition — is now available. The
justification is therefore no longer "we cannot", it is **"we must not"**: a score
built from structure count, clinical stage and length ranks targets by how much
attention they have had, and that is B1.

Clinical stage is newly available and is the most tempting thing in the output to
sort on. It is reported and **may not be compared, sorted, thresholded, or used in
candidate selection** (**B5''**).

Where an ordering is needed for display, it is the pool order **`stage4.build_pool`
produces**, which is `sorted(key=(-composite, gene, accession))` over Stage 3's
composites, deduplicated to one accession per gene and truncated at 200. Revision 1
called it "the Stage 4 order", which is right about where the sort happens and
silent about whose numbers it sorts; the values are Stage 3's and the ordering
operation is Stage 4's. Naming only one of the two is how J4 gets checked against
the wrong artifact. **That order is void across any Stage 3 recalibration**, and R13 (§0) is
exactly such a recalibration pending. The header records the hashes; this document
records that the order dies with them.

---

## 7. Rejection criteria — fixed in advance

Prefixed `B`. Stage 3's `R` and Stage 4's `P` criteria apply to the runs feeding
this one and are not re-checked here.

Revision 1's set had a structural flaw worth naming: it measured the stage against
numbers the stage itself produced, pinned no known answer to a retrievable
identifier, and routed every identity question — right protein, right isoform,
right chain, right species, request actually succeeded — through invariants, which
stop the run rather than being reported and are checkable only against the stage's
own mapping. A self-consistent wrong mapping satisfied all of them.

### Construction invariants

| id | invariant |
| --- | --- |
| J1 | no target appears in the output more than once |
| J2 | no candidate is emitted whose antigen chain does not map to the target |
| J4 | the pool order out equals the order in — this stage re-ranks nothing |
| J5 | every retrieval **records** its HTTP status and byte count, and no zero is ever inferred from an exception. This is a recording invariant: a non-2xx is legitimate and transient, so it does not halt the run — B12 judges what may be emitted under one |

J3 of revision 1 — "`NOT_CONNECTED` fields are never compared, defaulted or sorted
on" — was a statement about code paths that no verifier can evaluate from output.
It becomes **B5'** and **B6'**, which are checkable.

### The known-answer table

Fixed here, in advance, pinned to identifiers. Revision 1 referenced "a target with
a known clinical CAR binder" and never enumerated one, which meant the list would
be chosen from whatever the run returned.

The five known targets for this indication are taken verbatim from
`verify_ranking.py:KNOWN_TARGETS` — **CEACAM5, CEACAM6, CLDN18, MSLN, MUC1** —
rather than restated. Revision 1's §4.1 named a different trio including PSCA;
PSCA is GPI-anchored and belongs in the anchor criterion, but it is **not** one of
the five and has no named therapeutic.

| target | expected, structure route | expected, sequence route |
| --- | --- | --- |
| MSLN | amatuximab in 4F3F / 7UED / 8CXC; anetumab in 8CZ8 | amatuximab, anetumab, misitatug, inezetamab |
| CLDN18 | zolbetuximab in 9V32; osemitamab in 9V2U / 9V31 | 11 therapeutics incl. zolbetuximab (approved) |
| CEACAM5 | tusamitamab in 8BW0 | 5 therapeutics incl. tusamitamab |
| MUC1 | SM3 in 1SM3 / 6FZQ / 6FZR / 6TGG; AR20.5 in 5T6P / 5T78 | 3 therapeutics |
| CEACAM6 | **none — 6 entries, none with an antibody** | 1 therapeutic (tinurilimab) |
| PSCA | **none — 1 entry, no complex** | **none** |
| NRG3 | **none — 0 entries, HTTP 204** | **none** |

The negatives are as load-bearing as the positives; without them the checks are
one-sided. Two facts recorded so a correct implementation is not failed by them:
**the SS1 scFv — the actual clinical mesothelin CAR binder — is not in the structure
database** (five query phrasings, all 204), which is why the positive checks are
phrased at entry level rather than binder-name level.

The second is a worked example of why **B10 keys antigen identity on the accession
cross-reference and not on any name string**. CEACAM5 has 6 accession-keyed entries
and exactly one antibody complex, 8BW0 — verified as "Structure of CEACAM5 A3-B3
domain in Complex with Tusamitamab Fab". But that row's curated `antigen_name`
records a carbohydrate, so a stage matching antigens by name either drops the one
real candidate or keeps it for the wrong reason. CEACAM5 must return exactly one
candidate, from 8BW0, and must not return one from its five non-antibody entries.
The curated annotation is authoritative for *which chain is which*; the accession
cross-reference is authoritative for *which protein it is*. Conflating the two is
how a name string gets to decide identity.

### Criteria

| id | criterion |
| --- | --- |
| — | **B1 is deleted. See the note below: B3 already carries its failing state, and a criterion that cannot fail independently is the double-counting this document deletes B7 for.** |
| B2 | for any **structure-route** candidate of a GPI-anchored target, the ectodomain verdict is anything other than PASS — `FAIL` trips, and so does `NOT_MEASURED` **where contacts were computable** — that is the anchor trap. Where contacts could not be computed at all, or fewer than 3 were found, the candidate is outside this criterion and its count is reported instead; otherwise B2 and B3 contradict each other on the pinned MSLN entries. Sequence-route candidates carry no epitope and are outside it by construction (§2.4); a candidate flagged `SHED_ANTIGEN` trips, which is the intent |
| B3 | any pinned entry above is absent from its target's candidate list; or any documented-negative target returns a structure-route candidate; or **CEACAM5 returns other than exactly one structure-route candidate** — 6 accession-keyed entries, 1 complex, so an entry-count-echoing stage returns 6 and fails here. That last clause is what carries the null B1 was deleted for, and it is stated as a criterion rather than left in prose |
| B4 | entries needing the description-text fallback exceed those resolved by the curated annotation |
| B5' | any affinity-typed field is neither the literal `NOT_CONNECTED` nor a record carrying value, relation, units, assay id, source and source release |
| B5'' | any ordering, threshold or selection whose key includes clinical stage or therapeutic name |
| B6 | **the stage emits a combined binder score of any kind** — any per-candidate or per-target number that fuses two or more of the reported fields into one ordering quantity. Carried forward from revision 1 unchanged: §6 is the prohibition and this is the only criterion that enforces it |
| B6' | a candidate is singled out in any per-target field without its selection rule printed verbatim in the header as a deterministic total order over named measured fields, or that rule references a count, a search rank, clinical stage, or a `NOT_CONNECTED` field |
| B8 | fewer than 5 of the 5 known targets return at least one candidate across both routes |
| B9 | the printed construct components do not sum to the printed total, or any **non-binder** component in §4.3's table is absent from the sum, or headroom ≠ budget − total. The binder rows are mutually exclusive alternatives, so exactly one applies per arm |
| B10 | any candidate's antigen entity maps to a **different human** accession than the target's Stage 3 accession. An antigen mapping to a non-human **orthologue of the target** is flagged `NON_HUMAN_ANTIGEN` per §4.5(4) and does **not** trip; an antigen mapping to a non-human accession whose cross-referenced **gene symbol differs from the target's** trips, because "wrong protein" does not become a warning by being in another species. Gene symbol equality is the instrument — orthologues share it, and it is a cross-reference field, not a name string parsed out of a description — otherwise B10 would halt the run on the case B16 exists to warn about. An antigen entity with no cross-reference is `ANTIGEN_NOT_MAPPED`, never usable, and its count is reported |
| B11 | any **structure-route** candidate for a multi-isoform target is reported usable with an unresolved isoform; specifically, if CLDN18's structure-route candidates do not resolve to the .2 isoform. Sequence-route candidates carry no deposited antigen to align against, and the therapeutic table's target field records the **gene** (`CLDN18`), not the isoform — measured. They are therefore `ISOFORM_UNRESOLVED` unless the source names an isoform explicitly: reportable, and never `usable` |
| B12 | any target emits `NO_BINDER` while any of its requests returned a non-2xx status; the count of targets with any non-2xx must be 0 |
| B13 | output row count ≠ 200, or the output gene set ≠ Stage 4's pool gene set; the symmetric difference is printed in both directions |
| B14 | a candidate's binder chain carries the target's own accession, or its binder-chain classification evidence is `none` |
| B15 | **any** candidate's retrieved sequence disagrees with its entity record's own length or leading residues — checked for every candidate, not only the pinned ones, since a join fault on the other ~196 targets would otherwise feed wrong lengths into §4.3 undetected. The pinned entries are additionally hand-checkable, so any disagreement there is a defect rather than a rate |
| B16 | any candidate whose antigen maps to a human accession is flagged `NON_HUMAN_ANTIGEN`, **or any candidate whose antigen maps to a non-human accession is not so flagged** — the criterion is two-sided, because a missing warning is the direction that ships |
| B17 | a candidate is drawn from an entry whose method or title indicates a computed model |
| B18 | any candidate's antigen residue range was not derived from the entity alignment block; specifically, if 4F3F's epitope does not map to proteome residues 302–359 |

**B1 is no longer a correlation.** Revision 1 set it as Spearman above 0.9 against
"the target's total PDB entry count", then wrote, in the same document, that B1
would probably trip and that tripping "is not automatically a defect" — a narrative
explaining the output in advance, which the preamble forbids. A criterion with a
pre-written excuse cannot stop anything: a stage that counts entries and re-emits
them as candidates trips B1, invokes the excuse, and ships.

It also had no well-defined instrument. Full-text and accession-keyed counts differ
by up to 8x and reorder the pool, so whether the old B1 tripped was decided by an
unstated choice of query string. And the repo has **two** Spearman helpers that
disagree: `verify_ranking.py` averages ranks across ties, `verify_pairing.py`
assigns ordinal ranks with no tie correction. Over a vector that is mostly zeros —
which this one is — they give different answers, so a criterion phrased as a
correlation must also name which helper computes it.

So the correlation becomes **a mandatory reported header number**, computed both
ways, with tied ranks, over targets with at least one entry, printed beside the size
of the zero block and beside the pool-wide retention ratio — entries examined
against candidates kept.

**And B1 is deleted rather than reformulated.** Three attempts were made to give it
an independent failing state, and each collapsed the same way: every formulation of
"the filters are inert" is already decided by B3's pinned answers. CEACAM6 has 6
entries and no antibody complex, so it must come out `entries > 0` with zero
structure-route candidates; CEACAM5 has 6 entries and exactly 1 complex, so it must
come out with 1 candidate and not 6. **An inert stage fails B3 on both.** A separate
B1 restating that is the double-counting this document deletes B7 for, and keeping
it would inflate the "N of M clear" headline with a check that cannot fail.

The null B1 was reaching for is real and it does not go away — it is simply not a
rejection criterion. It is the reported retention ratio and the two correlations,
and §2.6's warning that structure count measures attention. The reader is told the
number; the run is stopped by B3.

**B8 moves from 3 of 5 to 5 of 5.** The old floor tolerated a total retrieval
failure on the flagship target: a stage whose MSLN retrieval was entirely broken
still cleared at exactly 3/5. With the sequence route connected, all five clear
with sequences — which is itself the sharpest demonstration that revision 1's
"therapeutic name and clinical stage are not recoverable" was load-bearing.

**B7 of revision 1 is deleted.** It restated J4, and a criteria list that
double-counts inflates the "N of M clear" headline that a handoff is read on.

### Explicitly not grounds for rejection

- a novel target returning `NO_BINDER` — that is §4.5 and it is expected
- the known targets dominating the candidate counts
- no VHH being available for any target
- PSCA returning nothing on either route — it has no named therapeutic, and this is
  carved out explicitly rather than left to fire and be argued away

---

## 8. Output

Header carries: the Stage 3 and Stage 4 configuration hashes verbatim and both
their criteria outcomes; the structure-database query date, the antibody-database
`last_update` and API version, and the therapeutic-table release as pins;
**§2.1's run-time probe record, including the client and TLS stack**; the pool size
and its provenance; **every threshold this document fixes** — the 4.5 Å contact
cutoff, the 3-residue minimum, and every coding length in §4.3 — with the full
component sum; the
fallback-classification counts from §5; both correlation values from B1's note with
the zero-block size; the selection rule if any; and Stage 5's own configuration
hash covering all of it.

Per target: Stage 4 outcome and partner, **each carrying an inline marker that P12
and P13 are tripped for them** — a header note does not travel with a column —
candidate counts **per route**, ectodomain verdict, `antigen_species`,
`binder_source_organism`, `humanisation_state: NOT_CONNECTED`,
`affinity: NOT_CONNECTED`, isoform assignment, the single-receptor budget
arithmetic, and for `DUAL` targets the two-arm arithmetic or `NOT_COMPUTABLE`
naming the missing arm. A `SINGLE` carrying a non-null partner reports the partner
and is costed as one arm.

**Which candidate the budget row sums is fixed here**, because otherwise B9 cannot
be evaluated from the output: it is the candidate with the **smallest CAR-converted
binder length**, ties broken by entry identifier ascending, and **the chosen
candidate is named in the row**. That is a deterministic total order over a
measured field, so it satisfies B6' rather than smuggling a ranking in. Choosing
the smallest is stated as a bound — it reports the most favourable budget the
retrieved set allows, which is the honest direction for a feasibility check to err
in, and it is labelled as such rather than read as a recommendation.

**No silent caps.** All retrieved candidates are written to file.

---

## 9. Open problems

**1. Affinity is genuinely gone, and it is the field that would matter most.** §4.2.
The current release dropped the columns the retired one carried, and no reachable
substitute has them. Worse, the right target is a window whose position depends on
surface copy number per cell, which this pipeline does not have anywhere — so even
with the data the window could not be placed.

**2. Humanisation state is not recoverable** from any source connected here, and
the deposited source organism is not a substitute. §4.4.

**3. This stage will find no structure for most of what Stage 4 recommends.** §2.6.
Every one of the 14 recommendations has at least one arm with no structures. That
is a true statement about the literature and a misleading one about the biology,
and the output has to carry both readings at once. The sequence route narrows this
for the known targets and not for the novel ones.

**4. Three upstream changes block this stage, and none of them is a Stage 5 change.**
Stage 1's dataset row must split (§2.5); Stage 4 must persist its decisions (§3.1);
and `uniprot.FIELDS` must gain `ft_chain` (§4.1),
which changes the proteome cache fingerprint and forces a re-fetch of all 20,431
entries. The last is the expensive one and should be sequenced first, because every
other change is cheap to redo and that one is not.

Isoform sequences are **not** part of this: they are fetched as a separate
accession-keyed cache entry precisely so the entry set and the 20,431 / 3,466
counts do not move. And `BinderFormat` is **not** widened — see §1.

**One upstream spec also needs a correction, not just code.** Stage 3 §5 states
that the anchored class reports "zero annotated outward residues". It reports
`None`; §4.1 gives the measurement. Stage 3's implementation already keys on
`None`, so this is a wording defect in the document rather than a bug in the
stage — but it is the wording revision 1 copied, and leaving it in place guarantees
the next stage to read Stage 3 copies it again.

**5. The isoform problem is upstream.** §4.1.1. Stage 5 can refuse to attach a
binder to an unresolved isoform, but the risk score that made CLDN18 look like a
lung-risk protein is Stage 3's, and correcting it requires isoform-level expression
evidence the pipeline does not currently ingest.

---

## 10. What changed from revision 1, and why

Recorded because the failure mode is the calibration.

| # | revision 1 | measured | consequence |
| --- | --- | --- | --- |
| 1 | probe table transcribed into the spec | client's CA bundle fabricated one verdict; one path 404s | §2.1 — probe at run time, record the client |
| 2 | antibody structure DB "not machine-readable" | documented REST API, 11.7 MB CSV, 21,914 rows | §2.2 — the stage is not built around a gap |
| 3 | CDR boundaries "not recoverable without a numbering tool" | served pre-computed under IMGT numbering | §2.3 |
| 4 | therapeutic name and clinical stage "not recoverable" | 636 KB CSV, 1,133 rows, **with VH/VL** | §2.4 — a second retrieval route |
| 5 | affinity unavailable *because the source is unreachable* | source reachable; **field retired**; bioactivity route returns 0 | §4.2 — right verdict, wrong reason, corrected so it is not undone |
| 6 | structure counts from full-text name search | distorted target-specifically, 7.7x too high to slightly too low; one target's count was 369 unrelated entries | §2.6 — accession-anchored |
| 7 | ectodomain check on the antigen **construct** | rejects an approved drug's target | §4.1 — check the epitope |
| 8 | lipid-anchor class as a three-name list "reports zero" | returns `None`; list omits CEACAM6, includes a non-member | §4.1 |
| 9 | mesothelin's mature chain "extracellular in its entirety" | 45% of it is a secreted chain | §4.1 |
| 10 | isoform not mentioned | CLDN18 risk peaks on the lung, on the wrong isoform | §4.1.1, B11 |
| 11 | budget: "two scFvs is tight" | mandatory safety switch omitted; two scFvs **do not fit** | §4.3 |
| 12 | B1 a correlation with a pre-written excuse | instrument undefined, failure pre-excused | §7 — inertness check |
| 13 | B3 referenced an unenumerated list | list would be chosen after the run | §7 — pinned table |
| 14 | identity checks as invariants | a self-consistent wrong mapping passes | B10–B16 |
| 15 | residue ranges taken as deposited | 4F3F numbers its antigen 7–64 for proteome 302–359 | §4.1, B18 |

Revision 2 was itself reviewed before being offered, and the following were found
in it and fixed rather than shipped — recorded because a revision that presents
itself as the corrected one has the most to gain from hiding its own defects:

| # | revision 2 draft | why wrong | fixed in |
| --- | --- | --- | --- |
| 16 | GPI pass rule deleted with the wrong wording that carried it | every anchored target would return `NOT_MEASURED` and trip B2 | §4.1 |
| 18 | B10 tripping on any accession mismatch | halts on the cross-species case B16 exists to warn about | B10 |
| 19 | "overstated 2x to 8x" | one count was exact and one was too low | §2.6 |
| 20 | isoform sequences as a proteome field | they change the entry set, moving 20,431 and 3,466 | §4.1 |
| 21 | only the dual construct costed | 187 of 200 targets are single and B9 applies to them | §4.3 |
| 22 | `Fab` added to `BinderFormat` | that enum types a project *input*, not a Stage 5 output | §1 |
| 23 | design verdict keyed on `partner is not None` | prices two binders for a single-antigen design | §4.5 |
| 24 | B1 kept at all, across three reformulations | every phrasing of it is already decided by B3's pinned answers | B1 deleted entirely |
| 25 | B11 with no sequence-route carve-out | CLDN18's 11 therapeutics would trip a correct run | B11 |
| 26 | GPI rule needing an anchor position the record drops | `ft_lipid` is reduced to a boolean | §4.1 |
| 27 | B16 one-sided | an unflagged non-human antigen passed | B16 |
| 28 | §2.5's ordering inverted | the guard fires on the *safe* order, not the unsafe one | §2.5 |
| 29 | `require_safety_switch` stated as a constant | it is derived from safety tolerance; the budget conclusion inverts without it | §4.3 |
| 30 | J5 halting on any non-2xx | made B12 unfailable, the same inflation B1 and B7 were deleted for | J5, B12 |
| 31 | J1 relaxed to "no more than once" | dropped targets became reportable rather than halting; B13 now carries completeness | J1, B13 |
| 32 | B4's threshold moved from 30% ambiguous to a fallback ratio | the denominator changed with it and the move was unlogged | B4 |

---

## Build note

Once approved, and with §0, §2.5 and §3.1 understood:

1. `uniprot.FIELDS` gains `ft_chain`, with the proteome
   re-fetch that fingerprint change forces; then Stage 1's dataset row and its connector
   move together as one change (§2.5); then Stage 4 gains its decision writer. Isoform sequences are a separate accession-keyed cache entry, and
   `BinderFormat` is not touched
2. `data/structures.py` and `data/antibodies.py` — search, entity metadata,
   sequences, coordinates, curated chain annotation, therapeutic table; cached
   under the usual discipline
3. the known-answer table (B3, B8) run before anything else is written
4. `stages/stage5.py` — retrieval by both routes, then isoform assignment, then the
   epitope check, then sizing
5. `verify_binders.py` — four invariants, then the criteria, then the biology

Same order and same rule as Stages 3 and 4.
