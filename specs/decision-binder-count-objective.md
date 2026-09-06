# Decision required: what `binder_count` should count

**Status: open. Awaiting a decision. No code has been changed.**

This changes which design the platform tells a reader to advance, so it is put
here rather than patched.

## The finding

Stage 5 retrieves binders by searching a structural database on the target's
accession. A hit means the entry *contains* the target. It does not mean the
antibody in that entry was raised against it: an entry may carry an antibody
against a different chain of the same complex.

The target-match check now measures how often that happens. Across the current
pancreatic pool, **83 of 315 structure-route binders — 26% — are annotated
against a different protein than the target they were retrieved for.** The
pattern is systematic rather than scattered: G-protein-coupled receptor entries
where the search finds the receptor–G-protein complex and the antibody in it is
annotated against the G-protein subunits. ADGRF1, ADGRG6, CASR and others show
the identical shape.

**The retrieval count did not change. What it means did.** Every structure-route
binder count this platform has reported was a count of database hits.

## Why it reaches a decision rather than a bug fix

`binder_count` is one of the four objectives the Pareto front is computed from,
alongside tumour attractiveness, safety margin and binder cleanliness. It counts
hits.

Among the five designs that ship, four retrieved no binder at all. The fifth is
GPR35, which retrieved exactly one: entry `8H8J`, whose recorded antigens are
*Guanine nucleotide-binding protein subunit alpha-13* and *Guanine
nucleotide-binding protein G(I)/G(S)/G(T) subunit beta-1*. GPR35 itself is not
listed as an antigen anywhere in its own entry.

| design | attractiveness | safety margin | binder_count | cleanliness | front | decision |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| FER1L6 | **0.5023** | **−0.1381** | 0 | 0 | yes | ADVANCE |
| GPR35 | 0.5022 | −0.1627 | **1** | 0 | yes | ADVANCE |
| TMEM92 | 0.4959 | −0.1381 | 0 | 0 | no | BACKUP |
| TNFSF9 | 0.4114 | −0.1391 | 0 | 0 | no | BACKUP |
| BTNL8 | 0.4059 | −0.1747 | 0 | 0 | no | BACKUP |

**FER1L6 beats GPR35 on attractiveness, beats it on safety margin, and ties on
cleanliness.** GPR35 is non-dominated on exactly one objective: a binder count
of one. That one binder is annotated against G-protein alpha-13.

Recomputing the front counting only target-matched binders:

```
front as computed                        ['FER1L6', 'GPR35']
front counting matches rather than hits  ['FER1L6']
```

GPR35 becomes dominated, leaves the front, and **its decision moves from
ADVANCE to BACKUP**. Nothing else in the pool moves.

## The options

**A — count matches.** `binder_count` counts only binders whose recorded antigen
names the target. GPR35 drops to zero, leaves the front, and reads BACKUP. The
platform advances one design rather than two.

**B — count hits, as today.** `binder_count` stays a retrieval statistic. GPR35
continues to advance on the strength of a binder that is annotated against
something else, and the package says so in its gaps section.

**C — carry both.** Report `binder_count` and `matched_binder_count` separately
and leave the objective on hits. Honest, and it defers the question: the front
is still computed from the number that measures retrieval luck.

## Recommendation: A, count matches

**A scoring objective should measure a property of the design, not a property of
the search that found it.** The other three objectives all do. Tumour
attractiveness measures the target. Safety margin measures residual risk below
the ceiling the design was judged against. Cleanliness measures sequence
liabilities in the binder itself. `binder_count` alone measures how many rows a
database query returned, and a query that returns a crystallisation chaperone
counts it the same as a therapeutic antibody.

Under option B the platform advances GPR35 for a reason that does not survive
being stated out loud: *this design is not beaten on every axis because one
database row exists, and that row is an antibody against a different protein.*
That is retrieval luck presented as evidence.

Option A is also the smaller claim. It does not assert that GPR35 is a bad
target — GPR35 keeps its attractiveness, its safety margin and its BACKUP
position. It asserts only that a binder annotated against G-protein alpha-13 is
not evidence about GPR35, which is the same thing the target-match check already
says everywhere else in the output.

**What A does not fix.** A target-match PASS is a claim about annotation, not
about a measured interface. Counting matches makes the objective mean *binders
whose recorded antigen names this target*, which is better than *rows returned*
and still weaker than *binders shown to bind*. Settling that needs the
coordinates that are not connected.

## If A is chosen

The change is small and its consequences are known in advance, which is why they
are written here before it is made:

- `binder_count` counts rows whose target-match verdict is `PASS`. Rows reading
  `UNKNOWN` do not count, because unknown is not evidence.
- Sequence-route binders read `UNKNOWN` for want of deposited antigen
  annotation, so under a strict reading they would stop counting too. **They
  should keep counting**: they are named therapeutics recorded against the
  target, and their UNKNOWN is about a check that does not apply rather than
  about doubt. The rule is therefore *count a binder unless its recorded antigen
  names something else*, which is one clause and not two.
- The Stage 11 configuration hash changes, so no cached run compares equal to a
  run under the old objective. That is correct and intended.
- The predicted outcome, fixed before the run: the front becomes `['FER1L6']`,
  GPR35 moves to BACKUP, and no other design changes. If the run lands elsewhere
  that is reported first.
