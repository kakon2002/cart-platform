# Stage 12 — the final candidate package

Written before `stages/stage12.py` exists. The criteria in §7 are fixed here and
committed before any result.

This is the thing a reader actually receives. Every input already exists: the
package assembles what stages 1 through 11 produced and adds nothing measured.
Its one new claim is the **gaps section**, and that claim is the reason the
document is careful.

---

## 1. What the reference document asks, and what can be answered

Section 11 lists twelve deliverables. Assessed against what the document asks
each one to *contain*:

| | deliverable | state |
| --- | --- | --- |
| 1 | Top 3–5 CAR-T constructs | partial |
| 2 | One conservative backup design | produced, as a refusal |
| 3 | One advanced or innovative design | produced |
| 4 | Complete sequence and domain map | partial |
| 5 | Target and binder evidence report | partial |
| 6 | Safety-risk matrix | partial |
| 7 | Structural report | **absent** |
| 8 | Functional predictions | **absent** |
| 9 | Manufacturability assessment | partial |
| 10 | Experimental validation plan | produced |
| 11 | Model confidence and applicability domain | produced |
| 12 | Full evidence and decision audit trail | partial, closable here |

Eight of the twelve have something real to carry. The package carries those and
**names what the other four are missing**, per deliverable, rather than shipping
eight sections and letting a reader assume twelve.

### 1.1 Correcting the record on Stages 7 and 8

The project handoff states that Stages 7 and 8 are *"schema only by decision"*.
**That is wrong and this specification corrects it.** There is no schema. There
is no module in `car_pipeline/stages/`, no dataclass, no field, and no stub
anywhere in `car_pipeline/schemas/`. The numbering gap is recorded in prose —
`README.md` says the two are *"absent from the pipeline, not from the
numbering"* — and nothing is stubbed behind that sentence.

The distinction matters for exactly the reason this repository keeps
rediscovering: a schema with no implementation is a shape a reader can inspect
and a caller can fail against, and its absence is visible. Prose describing a
schema that does not exist is a claim nobody can check. Two of the twelve
deliverables rest on these stages, and a reader deciding what to build next is
entitled to know that the ground is bare rather than half-laid.

The two are not equivalent. **Stage 7 is buildable from what is connected** —
structure prediction over sequences the construct stage already emits. **Stage 8
is not**: the reference document itself names partner-generated experimental
data among its required training inputs, and no such data is connected.

### 1.2 The cheapest closable gap is inside deliverable 6

Stage 9's required safety output includes a genomic and construct-safety arm:
recombination-prone regions, cryptic splice sites, unwanted open reading frames,
sequence repeats, and editing-related risks for allogeneic products. None of it
exists.

**None of it needs a new stage.** All but the editing item are sequence analysis
over the DNA map `stage6` already emits — a construct whose nucleotide sequence,
domain boundaries and provenance are all in hand. No new data source, no model,
no external call. Of everything missing across the six partial deliverables,
this is the shortest path from absent to produced, and the gaps section says so
in those terms rather than listing it flat beside items that need a data source
nobody has.

This specification does not build it. It records that the cost is small and the
blocker is nobody's decision, so the next reader can price it correctly.

### 1.3 The conservative backup is a refusal, not an omission

`design_class` fixes what conservative means against the architecture table —
the conventional single-antigen receptor carrying a binder with clinical
precedent — and fixes it independently of what survived, which is the decision
recorded at `validation.py:219`. In the worked indication no design meets it.

The package therefore carries the section with a **counted refusal**: no
conservative design exists in this pool, three single-antigen targets were
recommended and none assembles for want of a binder, no dual design assembles at
all because every dual recommendation names a partner that retrieves no binder.

**A blank section and a refusal are different claims** and the difference is the
whole point of the deliverable. Blank reads as "not filled in yet". The refusal
reads as "we looked, here is what we found, and here is the count behind it".
Criterion Q7 exists so the section can never silently become the first.

## 2. Where the package lives

`car_pipeline/stages/stage12.py`, assembling from the objects
`car_pipeline/api/pipeline.py` already returns. It computes no new measurement.
Every number in a package is carried from the stage that measured it, and the
criteria in §7 assert that carrying is lossless rather than assuming it.

Served at `GET /projects/{id}/package` for the set and
`GET /projects/{id}/package/{gene}` for one. Rendered to
`reports/packages/<GENE>.md`, one file per candidate, by `make_package.py` —
the reader-facing artifact.

## 3. What one package contains

Per candidate, eight sections carrying what exists:

| section | source | deliverable |
| --- | --- | --- |
| `ranking` | Stage 11 `Ranked`: position, the four objectives, `on_front` | 1 |
| `design_class` | `validation.design_class` | 2, 3 |
| `construct` | Stage 6: amino-acid sequence, DNA map, domain boundaries, per-part provenance, budget arithmetic | 4 |
| `target_evidence` | Stage 3 `Ranked` plus the per-organ `risk_attribution` | 5, 11 |
| `binders` | Stage 5 `TargetBinders`, reported by route and never summed | 5 |
| `safety` | Stage 9 `SafetyRecord` | 6 |
| `developability` | Stage 10 `Developability` rows, flags listed and never scored | 9 |
| `validation_plan` | `validation.plan` | 10 |

Plus two that belong to the package rather than to any one stage:

| section | contents | deliverable |
| --- | --- | --- |
| `provenance` | every connected source with its release pin, and the configuration-hash chain from Stage 3 to Stage 11 | 12 |
| `gaps` | §4 | 7, 8, and the missing halves of 1, 4, 5, 6, 9, 12 |

**No section is emitted for a stage that does not exist.** There is no
`structural_report: null` and no `functional_predictions: []`. A null field reads
as computed-and-found-nothing; the absence belongs in `gaps`, where it is stated
once with its reason. Criterion Q9 enforces this.

### 3.1 Provenance, which closes most of deliverable 12

The document defines the audit trail as a versioned evidence graph linking every
recommendation to its dataset, publication, release date and confidence level.
Today the release pins exist as module constants and inside cache manifests, and
the configuration hashes chain stage to stage — and **no endpoint serves any of
them**. The trail a reader can reach names its stages and not its data.

The package emits both: each source with the release it was pinned to, and the
hash chain, so a package identifies the exact configuration that produced it.
That does not make an evidence graph and does not link publications; those stay
in `gaps` under deliverable 12. It does mean a candidate can be traced to the
data release behind it, which is the part that was missing and reachable.

## 4. The gaps section, and why it is the risky part

Everything else in the package is carried from a stage that measured it. The
gaps section is the one place the package makes an assertion of its own — *"the
platform does not produce X"* — and an assertion nobody recomputes is how the
narrative in `validation.py` came to claim three single-antigen targets were one,
and how a live endpoint came to report a criteria count that had been wrong for
weeks.

So each gap entry carries a **probe**: a declarative statement of what would have
to exist for the gap to be closed, in a form the verifier executes.

| probe kind | closes when |
| --- | --- |
| `module` | the named module exists under `car_pipeline/stages/` |
| `field` | the named dataclass gains the named field |
| `key` | the named key appears in a package section |

A gap whose claim is a judgement rather than a presence — deliverable 5's
recommended-status vocabulary differing from the document's five classes — carries
**no probe and says so**. Criterion Q6 checks every probe that exists and reports
how many entries were checkable, so the split between verified and asserted is
visible rather than implied.

**This is the mechanism that stops the gaps table becoming the next stale
sentence.** Build Stage 7 and forget to update the table, and Q6 trips.

## 5. What the package does not do

- **It computes nothing.** No new score, no new threshold, no re-derivation. If
  any figure in a package differs from the stage that produced it, that is a
  defect in the assembly, and Q4 and Q5 are there to catch it.
- **It does not rank.** Order is Stage 11's, carried through.
- **It does not soften a gap.** A missing element is named, not described as a
  future enhancement or a limitation of scope.
- **It does not fill the conservative backup.** §1.3.
- **It changes no existing stage.** The only edits outside `stage12.py` are the
  hash chain and release pins added to the run dict, the two views, and the
  renderer. If any existing criterion count or hash moves, that is a defect in
  this change, not a result.

## 6. Order of work

1. This document, committed before any code.
2. `stages/stage12.py`: the package, the gap table with its probes.
3. `api/pipeline.py`: the hash chain and release pins onto the run dict.
4. `api/server.py`: the two views.
5. `verify_package.py`: the criteria below, then the biology.
6. `make_package.py`: one markdown artifact per candidate.
7. `run_all.py` gains the stage.
8. Full suite, both indications, every count reported beside the count predicted.

## 7. Rejection criteria — fixed before the run

Identities and positive pins. None is a threshold fitted to an observation.

| id | trips when |
| --- | --- |
| **Q1** | the package set is empty, or its gene set differs from Stage 11's survivor set. Stated so an empty set **fails**: a run with no survivors must report that as a status, as Stage 11 does, never as an empty list of packages |
| **Q2** | a package exists for a candidate that did not reach the end, or a candidate that reached the end has no package. Both directions, because either alone passes on a stage that emits nothing |
| **Q3** | any package omits one of the eight carried sections, or carries one that is empty while the stage behind it produced something for that gene |
| **Q4** | for any package, the DNA in the package does not translate to the amino-acid sequence in the package, or the domain boundaries do not partition it exactly. The construct stage already asserts this; Q4 asserts it on the **copy**, which is what a reader receives |
| **Q5** | for any package, the per-organ attribution does not reconstruct the packaged risk to within 1e-12, or the packaged risk differs from Stage 3's |
| **Q6** | any declared gap carrying a probe is not actually open — the module exists, the field exists, the key is present — or any gap entry claims an element the package in fact carries. Positive pin: trips if **no** probe executed, so a table of unprobed assertions fails |
| **Q7** | the conservative-backup section is absent or empty, or reports a backup that `design_class` did not label conservative, or reports a refusal without the counts behind it |
| **Q8** | any connected source is missing its release pin, or the configuration-hash chain has a gap between Stage 3 and Stage 11 |
| **Q9** | the package emits a section, key or placeholder value for Stage 7 or Stage 8 rather than recording their absence in `gaps` |

**Q6 is the one that would catch this document going stale**, and Q1 is the one
that would catch the stage doing nothing. Q4 and Q5 are the two that would catch
a lossy assembly, which is the only kind of defect an assembler can have.

### Explicitly not grounds for rejection

- The conservative backup being a refusal. §1.3.
- Deliverables 7 and 8 being absent. That is the finding, recorded in `gaps`.
- A package carrying fewer than twelve deliverables' worth of content. Eight is
  what exists; claiming twelve would be the defect.
- The gaps section being long. It is the honest length.
