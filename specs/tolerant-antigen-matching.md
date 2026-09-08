# Specification: tolerant antigen-name matching

**Status: specified, not implemented. No code has been changed.** Criteria are
fixed here before anything can be checked against an observation, and the
predicted outcome is written down before the run that would produce it.

This is a separate change from the binder-count correction in
`binder-count-correction.md`, deliberately. That correction was reviewed and
shipped on the figures the current matching rule produces. Folding a second
correction into it would have made neither reviewable.

## The finding

The target-match check derives a set of names for a target from the reference
record — the recommended name, each parenthesised synonym, the note on each
mature chain, and the gene symbol — and compares each recorded antigen against
that set by **exact string equality after lowercasing**.

Exact equality is defeated by three label shapes that do not indicate a
different molecule. Measured across the pancreatic pool, **roughly 26 of the 83
flagged structural rows** are the target itself under a label the comparison
could not see:

| shape | recorded antigen | target | rows |
| --- | --- | --- | ---: |
| isoform label | `Isoform A2 of Claudin-18` | CLDN18 | 4 |
| isoform label | `5'-nucleotidase, ecto (CD73), isoform CRA_a` | NT5E | 5 |
| word-order variant | `Receptor protein-tyrosine kinase erbB-2` | ERBB2 | 3 |
| fusion-partner label | `Ubiquitin-like protein SMT3,Cadherin-1` | CDH1 | 1 |
| fusion-partner label | `Maltose/maltodextrin-binding periplasmic protein,Mucin-16` | MUC16 | 2 |
| fragment label | `MUC1 Peptide Fragment` | MUC1 | 2 |
| construct label | `Inactive TMPRSS2 construct` | TMPRSS2 | 1 |

Nine targets are affected: **CDH1, CLDN18, ERBB2, F2RL1, MUC1, MUC16, NT5E,
TMPRSS2, TREM2.**

The CLDN18 case is the one that matters most. Isoform A2 is claudin-18.2 — the
clinically validated form of that target and the antigen of an approved
therapeutic. The check calls four binders against it wrong-antigen.

**The error is safe in direction and still wrong.** It withholds credit from a
real binder, which can only make the platform more conservative. It cannot
advance a design on evidence that does not exist. No design that reached the
ranking is affected, and the binder-count correction stands unchanged.

## What this must not become

This is the sixth string-matching failure on record in this project. The other
five all ran the opposite way — a match that was too permissive:

- localisation phrases matched as substrings of free-text notes, admitting
  *"not detected at the cell membrane"* as evidence of membrane localisation;
- a compound antibody-target field matched as a substring, putting MUC16 and
  MUC18 into the MUC1 bucket and matching CLDN1 inside CLDN18;
- tissue names mapped to organs by substring, putting the adrenal at kidney
  criticality and folding the kidney into the brain;
- a free-text name query returning several hundred entries topped by a bacterial
  RNA chaperone, a sulfur transferase and a photosystem supercomplex;
- a blinding guard enumerating modules by name, behind which thirteen modules
  went unexamined.

The standing rule those produced — match on exact values and explicit token
boundaries, never on substrings — is what caused this sixth one. **The fix is
not a retreat to substring matching.** Any implementation that would
reintroduce any of the five above is rejected regardless of how well it handles
the seven rows in the table.

### The latent hazard that must be closed first

The derived name set contains degenerate entries. Measured over the full surface
set of 3,466 records:

```
distinct degenerate names (<=3 characters, or containing no letter):  380
name-set entries affected:                                            612
targets carrying at least one:                                        486
```

The most common is `'+'`, which appears in the name set of **115 targets** —
extracted by the synonym rule from constructions such as
`Na(+)/K(+)-transporting ATPase`. Others are bare `'A'` (20 targets), `'TM'`,
`'SU'`, `'2+'`, `'-'`.

Under the current exact-equality rule these are harmless: a recorded antigen is
essentially never the single character `+`. **Under any containment or
subsequence rule they are catastrophic** — `'+'` would match every antigen
string containing a plus sign, and `'A'` would match nearly everything.

That immunity is an accident of the current rule, not a property of the data.
Excluding degenerate names is therefore a **precondition** of this change, not a
detail of it, and it carries its own criterion.

## What a tolerant match must accept

The three shapes above, stated as transformations of a reference name:

1. **An isoform or variant qualifier** wrapped around the reference name —
   `Isoform A2 of X`, `X, isoform CRA_a`.
2. **A fusion or expression partner** listed alongside the reference name,
   comma-delimited within a single antigen element.
3. **A word-order variant** of the reference name using the same words —
   `Receptor protein-tyrosine kinase erbB-2` against the reference's
   `Receptor tyrosine-protein kinase erbB-2`.

Fragment and construct qualifiers (`Peptide Fragment`, `Inactive … construct`)
fall under the same treatment as (1).

## What a tolerant match must refuse

The prohibition is the substance of this specification. Every case below is real
and present in the current pool, and every one of them is correctly `FAIL`
today. **A tolerant match that turns any of them into a `PASS` is worse than the
error it fixes**, because it would manufacture evidence rather than withhold it.

| target | recorded antigen | why it must stay flagged |
| --- | --- | --- |
| MET | `Hepatocyte growth factor beta chain` | the **ligand**, where the reference name is `Hepatocyte growth factor receptor`. Every word of the antigen appears in the reference name |
| ERBB3 | `Receptor tyrosine-protein kinase erbB-2` | a sibling receptor. The same string must match ERBB2 and refuse ERBB3 |
| ITGB6 | `Integrin alpha-V heavy chain` | the partner chain of the same heterodimer |
| ITGA6 | `Integrin beta-1` | the partner chain of the same heterodimer |
| CLDN18 | `Claudin-1` | a different claudin, carried separately in the same pool |
| MUC1 | `Mucin-16` | a different mucin, carried separately in the same pool |
| GPR35 | `Guanine nucleotide-binding protein subunit alpha-13` | the stabilising G-protein, not the receptor |

The MET case is the sharpest and should be treated as the design constraint
rather than as one test among seven. A rule that discards generic words to
improve recall — dropping `receptor`, `protein`, `chain`, `subunit` — matches
the hepatocyte growth factor **ligand** to the hepatocyte growth factor
**receptor**. Those are different molecules, and an antibody against one is not
evidence about the other. **Relational and structural words carry meaning and
must not be discarded.**

Note that the screening heuristic used to produce the figure of roughly 57 in
the case study *did* discard `protein` and `isoform`. That was acceptable for
after-the-fact triage of an already-flagged set. **It is not acceptable as an
implementation**, and it is recorded here so the two are not confused.

## Criteria

Fixed before implementation. Positive pins wherever a criterion could otherwise
clear by finding nothing.

| id | criterion |
| --- | --- |
| **M1** | All seven accept-cases in *What a tolerant match must accept* move from `FAIL` to `PASS`. **Positive pin:** all seven must move. A partial rule that fixes isoform labels and not word order is not this change. |
| **M2** | All seven refuse-cases stay `FAIL`, MET included. Asserted case by case, not as an aggregate rate. |
| **M3** | No reference name shorter than four characters, and none containing no letter, is used for tolerant matching. **Positive pin:** the exclusion must be shown to fire — at least one target in the pool must have a name excluded by it, or the guard is untested. |
| **M4** | For every pair of pool targets whose reference names differ only in a trailing number, an antigen naming one does not match the other. **Positive pin:** the pool must contain at least one such pair, and the criterion names the pairs it found. |
| **M5** | An antigen consisting of a reference name plus a relational or structural word — `receptor`, `ligand`, `binding`, `associated`, `subunit`, `chain`, `heavy`, `light`, `partner` — does not match, unless the reference name itself contains that word. MET is the live instance. |
| **M6** | Every row that reads `PASS` under exact matching still reads `PASS`. Tolerance may only add matches, never remove them. Monotonicity makes the change reviewable as a strict relaxation. |
| **M7** | The flagged count moves and the movement is reported against its prediction. **Positive pin:** if the flagged count is unchanged, the rule did not take effect and this trips rather than clears. |
| **M8** | If any candidate's ranking decision changes, both decisions are carried on the record in the same four-link chain the binder-count correction established. A decision that moves silently is a defect regardless of which direction it moves. |

## Predicted outcome, fixed before the run

- Flagged structural rows fall from **83 of 315 to approximately 57 of 315**.
- Affected candidates fall from **25 of 200 to approximately 21 of 200**.
- The nine named targets gain counted binders. CLDN18 goes from 11 counted of 15
  to 15 of 15, and its counted set stops being entirely sequence-route.
- **No design that reached the ranking changes.** None of the five is affected,
  GPR35's flagged binder is a genuine mismatch, and the other four retrieved no
  binder at all. The Pareto front stays `['FER1L6']` and GPR35 stays `BACKUP`.
- The Stage 11 configuration hash changes, because the counted objective changes
  for nine targets. This must be confirmed by computing both, not by inspection
  — the assertion that a hash would move has already been wrong once in this
  project, in `decision-binder-count-objective.md`.

**If the run lands anywhere else — and in particular if any ranking decision
moves — that is reported before anything else, and this specification is amended
and the run repeated rather than the result explained.**
