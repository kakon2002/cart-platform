# Verification that shares an assumption with the thing it verifies

Read this before writing a criterion. It has been found six times in this
repository. The fifth was found by suspecting the fourth, and the sixth was
introduced by a change made after the first five were written down.

## The pattern

A check is written to test something. It draws on the same value, the same
constant, the same list or the same fitting procedure as the thing under test.
It then reports success no matter what the thing does, because the two move
together. Nothing in the output looks wrong. The number is often unusually
clean, and cleanliness is the only visible symptom.

**The check that would have caught every instance below is one question: what
would this criterion report if the thing it tests were broken?** If the answer
is "the same as it reports now", the criterion is not a test. Ask it before the
criterion is written, not after it passes.

## The six

**1. The summary that was a literal.** Each verifier printed
`{9 - len(tripped)}/9 criteria clear`. The denominator was a constant, not a
count of the criteria that ran. Adding a tenth criterion left the summary saying
nine, and the runner preferred the printed summary to the criteria it had
parsed, so the new criterion was invisible in the total while sitting in the log
two lines above it. Nine of ten verifiers carried the same literal.

*Broken-thing question:* if a criterion silently stopped running, the summary
would report the same nine.

**2. The criterion cited as a rationale and never written.** R13's accepted
rationale in `run_all.py` said it had been "Replaced by R13-prime, which
clears." No criterion with that identifier existed anywhere in the code. The
exemption stood on a replacement that was never implemented, and the citation
was load-bearing: it is why R13 was allowed to keep failing.

*Broken-thing question:* if R13-prime did not exist, the rationale would read
exactly as it did.

**3. The check sharing a skip list with its subject.** The comment strip skipped
directories named in `SKIP_DIRS`, which contains `data` for the cache. The test
was every path component, so the source package `car_pipeline/data/` matched and
was skipped: fifteen modules, three hundred comments, untouched. The
verification that confirmed "zero comments remain" reused the same `SKIP_DIRS`
logic, so it confirmed its own blind spot and reported a clean sweep.

*Broken-thing question:* if the strip skipped a whole package, the check would
report zero comments remaining.

**4. The statistic that was the fitting procedure's definition.** R13-prime
tested whether staining agrees with transcript, reporting a median difference of
exactly +0.0000 and a directional split of exactly 50.0%. `calibrate_atlas_levels`
sets each staining level's value to the **median** of the paired baseline values
at that level, over the same organ-pairing construction the criterion then
measures. The baseline score is monotonic in TPM, so "staining exceeds
transcript" is exactly "the median exceeds the observation", which is 50% by
definition. Measured per level: 2,228/2,228, 4,035/4,035, 2,294/2,294. Exact
ties were 2 of 17,116, so ties were not the mechanism; the median fit was.

*Broken-thing question:* if the calibration were badly mis-centred and then
refit, the criterion would still report exactly 50%.

This one was invisible from the numbers. A median of +0.0000 and a split of
50.0% look like strong agreement and are a tautology.

**5. The limit that arrived with the criterion.** R13's 5x limit and R13 itself
first appear in commit `c21794c`, which rewrote 1,931 lines of `stage3.py`,
1,029 of the stage 3 spec and 934 of the verifier. No independent derivation of
5x exists anywhere in the history. The same is true of R13-prime's bounds: the
±0.05, the 35–65%, the observed `30,906` and the observed `50.0%` all first
appear together in commit `89ea896`, a single-file change, with the spec
asserting in that same commit that the bounds were "set from what the scoring
function does rather than from the observed value".

*Broken-thing question:* a bound chosen after seeing the value it bounds will
always be satisfied by that value.

**6. The origin check blind to the class it most needs to see.** Stage 9 flags
non-human sequence content in a binder, and it determines species from the INN
name stem: `-xi-` chimeric, `-o-` murine, `-u-` human. When the adaptor
receptor gained a binder retrieved from a deposited structure, that binder had
no INN name, because it is not a named therapeutic. The check therefore could
not fire on it. Worse, the usable-binder test looked only at stage 5 records,
which an adaptor design has none of, so seven of the eight surviving designs
returned `NO_GATE` with the reason "no binder, so there is nothing to gate" —
on designs whose receptor carried a murine scFv.

*Broken-thing question:* if a construct carried a wholly non-human binder from
a source with no naming convention, the check would report `ORIGIN_UNKNOWN` and
the gate would report that there was nothing to gate.

This is the same family as instance 3. There, the check shared a skip list with
the code it checked. Here, the check shares a naming assumption with the only
kind of binder it had ever been given. Both are blind spots aligned exactly
with the case that matters, and in both the output looks clean: `NO_GATE` reads
as "not applicable", not as "I cannot see this".

The fix reads the species from the deposition rather than from a name, and
refuses an entity that declares no source organism, so the value cannot be
quietly defaulted. The criterion that verifies it, A12, had to be moved above
the summary line before it counted — as first written it ran after the tally,
so the verifier would have printed 10/11 while reporting twelve criteria. That
is instance 1 reappearing as a placement rather than a constant, in the same
change that closed instance 6, and it was caught by reading the output rather
than by the suite.

## What to do instead

Write the bound before the run and commit it separately from the result. If the
spec commit and the result commit are the same commit, the exercise has already
failed. Where no defensible bound can be derived in advance, say so and report
the statistic ungated, stating that the property is uncharacterised. That is a
weaker claim honestly made, and it is worth more than a threshold fitted to an
observation.

Derive the bound from what the scoring function does, not from what the data
did. R14 is the current example: it asks whether any criticality tier places two
staining levels on opposite sides of the gate. That mentions no observation and
could have been written before the first run. It fails, and it fails
informatively.

## Structural notes from this line of work

**The staining axis takes four values.** Level 0, plus three calibrated points
at 6.3168, 12.6750 and 23.0102 TPM-equivalents. It is a coarse ordinal mapped
onto a continuous scale. Any future test of that axis is testing a three-point
ordinal, and should not demand agreement finer than three points can express.
It is part of why the staining-versus-transcript interdecile disagreement is
0.549 while the centre is pinned at zero.

**The arm-switch statistic cannot separate information from noise.** Of
protein-confirmed targets clearing on the transcript arm alone, 72 of 83 are
blocked once staining is added. That set is predicted exactly by presence
alone — whether the target carries any positive staining call in a tier-1 or
tier-2 organ — with 72 predicted against 72 observed and no discrepancy in
either direction. The statistic does not depend on staining magnitude, so a
bound on it would be instance six. It is reported ungated for that reason.
