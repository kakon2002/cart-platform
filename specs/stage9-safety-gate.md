# Stage 9 — safety gate

Written before `stages/stage9.py` exists. Short, like Stage 6: this stage
aggregates measurements other stages made and adds two of its own. It decides
nothing that Stage 3 already decided.

Status: **implemented**. `stages/stage9.py` and its verifier exist and the criteria below were fixed before the run, not after it. The gate this line used to carry has been passed, and leaving it in place would have described the repository incorrectly.

---

## 1. What this gate is for, and what it is not

Stage 3 already gates on normal-tissue risk and Stage 4 on the pair. **This stage
does not re-open either.** It asks the questions those stages could not, because
they had no binder and no construct to ask about:

1. **On-target, off-tumour** — carried from Stage 3, not recomputed. A gate that
   re-derived it would be a second place for the tissue-mapping bugs to live.
2. **Anti-CAR immunogenicity** — a property of the binder, which only exists from
   Stage 5.
3. **Prior clinical outcomes** — what has already happened to people dosed against
   this antigen.

A target failing here has a construct that should not be built. A target passing
here has not been shown safe; it has failed to be shown unsafe by three specific
questions, which is a weaker statement and is reported as one.

## 2. Sources

| question | source | verdict |
| --- | --- | --- |
| off-tumour risk | Stage 3 | carried |
| binder origin | the therapeutic's own name | **derived, see §3** |
| prior outcomes | trial registry, v2 API | connected |
| epitope-level immunogenicity | epitope database | **`NOT_CONNECTED`, see §4** |

## 3. Binder origin is read from the name, and that is a convention not a measurement

International non-proprietary names encode the source species in a stem:
`-o-` murine, `-xi-` chimeric, `-zu-` humanised, `-u-` human. Amatuximab and
zolbetuximab are both `-xi-`, and a chimeric binder carries a real anti-CAR
immunogenicity risk that a human one does not.

**This is a naming convention, not a sequence measurement**, and the output says
so on every row. It can be wrong in both directions: a name assigned before a
molecule was re-engineered keeps its original stem, and a binder with no INN has
no stem at all. Stage 5 already records `humanisation_state: NOT_CONNECTED` from
the deposited construct's taxonomy; this is a second, independent and equally
indirect reading, reported beside it rather than replacing it.

Where no stem is recognisable the value is `ORIGIN_UNKNOWN`. It is never guessed
from the sequence.

## 4. Epitope-level immunogenicity is not connected, and the reason is not reachability

The epitope database is reachable and queryable — its search endpoint answers on a
linear sequence. **The question this stage would need answered is not a lookup.**
Establishing whether a binder contains known immunogenic epitopes means scanning
every k-mer of a ~240-residue variable region against the epitope table, which is
a bulk download and a scan, not a query.

So: **the source is connected and the question is not answered here.** Every
binder carries `epitope_immunogenicity: NOT_CONNECTED`. Recording it this way
because "we did not check" and "we could not check" are different statements, and
here the first is true — which is the direction that must not be dressed up as the
second.

## 5. Prior outcomes

Per target, from the trial registry: how many interventional trials name the
antigen, their phases, and — the part that matters — **how many were terminated,
withdrawn or suspended**. A terminated trial is not proof of a safety problem and
is not reported as one; it is a signal that something is worth reading before
dosing another person.

Matching is on the antigen name in the query and is **inexact by construction**.
The count is reported as trials *mentioning* the antigen, never as trials *of a
binder against* it, because the registry's free text cannot support the second.

## 6. The verdict

Per target, one of:

- **`BLOCKED`** — Stage 3's risk exceeds the project ceiling. Carried, not
  re-decided.
- **`FLAGGED`** — passes Stage 3, but carries a chimeric or murine binder, or has
  terminated trials against the antigen. Both are reported with their reason.
- **`NO_GATE`** — no binder, so there is nothing to gate. Not a pass.
- **`PASSES_STATED_CHECKS`** — passes all three. Named at length deliberately:
  `SAFE` would be a claim this stage cannot make.

## 7. Criteria — fixed before the run

Positive pins first. Stage 5's lesson was that a route can be dead while every
negative check passes.

| id | criterion |
| --- | --- |
| **S1** | **the registry returns no trials for MSLN or CLDN18** — both have many; zero means the route is dead, not that the antigen is untried |
| **S2** | **amatuximab or zolbetuximab is not classified chimeric** — both are `-xi-`, and this pins the origin rule to a known answer |
| S3 | any target with a Stage 5 binder carries no Stage 3 risk value |
| S4 | any target whose Stage 3 risk exceeds the ceiling is not `BLOCKED` — this gate must never contradict Stage 3 |
| S5 | any binder carries a numeric or non-`NOT_CONNECTED` epitope immunogenicity value |
| S6 | output row count ≠ 200, or its gene set ≠ Stage 4's pool |
| S7 | a target with terminated trials is not flagged |

### Explicitly not grounds for rejection

- most targets returning `NO_GATE` — there is no binder for 134 of them
- a chimeric binder being flagged rather than blocked — immunogenicity is a
  managed risk, not a disqualification, and Stage 5 §4.4 already refuses to filter
  on it

---

## Build note

1. `data/trials.py` — registry counts per antigen, cached
2. S1 and S2 run before anything else
3. `stages/stage9.py`, then `verify_safety.py`
