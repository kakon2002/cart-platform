# Verification that shares an assumption with the thing it verifies

Read this before writing a criterion. It has been found nine times in this
repository. The fifth was found by suspecting the fourth, the sixth was
introduced by a change made after the first five were written down, the seventh
is the root cause of the third, found only when the third was fixed, the eighth
was found by reading a criterion's own output rather than its verdict, and the
ninth was found while fixing the eighth, in a verifier nobody was looking at.

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

## The nine

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
was skipped: thirteen modules, three hundred comments, untouched. The
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

**7. The name-based match that reached past its subject.** Instance 3 recorded
the effect: the strip skipped `car_pipeline/data/` and the check reused the same
predicate. Fixing it exposed the cause, which is not about caches or comments at
all. `SKIP_DIRS` held the *name* `data`, and the test asked whether any path
component equalled it. A name is not an address. One directory was meant - the
cache at `data/` - and the rule selected every directory in the tree that
happened to share the label, which is how a package of thirteen modules ended up
inside a skip list written for a cache.

*Broken-thing question:* if the rule selected the wrong directories, the sweep
would report a clean pass over the ones it did visit, which is exactly what it
reported.

This is the same shape as two matches already recorded elsewhere in this
repository: `renal` matching **adrenal gland**, and `cortex` matching
**Kidney_Cortex**. In all three the string is a real name for the intended
subject, the match is textual rather than structural, and the surplus is silent
- an adrenal gland is not a kidney, `Kidney_Cortex` is not every cortex, and
`car_pipeline/data/` is not the cache. Nothing errors, because over-matching
produces a *larger* answer that still typechecks. The fix in each case is the
same: match the thing by what identifies it - a resolved path under a known
root, an identifier, a column key - not by a word that appears in its name.

The reason this belongs in a note about verification is that a name-based rule
is what a verifier is most likely to reuse. It is short, it reads correctly, and
copying it into the check is the obvious way to keep the two consistent. The two
are then consistent about the wrong set.

**8. The five that cleared on an empty set.** The construct stage's K1, K3, K4,
K5 and K6 each iterate the constructs that were assembled. The artifact they
read was written by a call that passed no tolerances, so routing was disabled,
no adaptor row could exist, and nothing assembled. Five loops over an empty list
produced five empty failure lists and five clear verdicts. K4 printed *"0 parts
in the first construct"* and cleared on that sentence.

*Broken-thing question:* if the construct stage assembled nothing at all — which
is exactly what it was doing — every one of the five would report success.

The number was in the output the whole time. `0` is a count, it typechecks, and
a criterion phrased as *"no construct fails X"* is satisfied by there being no
construct. Only K2 tripped, and only because a previous amendment had already
given it a positive clause; the other five had none. **A criterion phrased over
a population must say what it does when the population is empty**, and the
answer must be that it fails.

The fix is in `specs/stage6-routed-decision-set.md`: the artifact is now written
routed, so there is something to read; K0 asserts the set is routed and non-empty
before the others run; and each of the five trips on an empty population in its
own right, saying so in its own words. Both halves are wanted. K0 makes the
failure legible in one line, the five clauses make it true even if K0 is ever
removed. Handed an empty set the stage now reports 2 of 9 clear where it used to
report 7 of 8.

This is the third of a sub-family: a criterion that never executes the path it
claims to cover. The first was M4 and M5, which asserted a field of an object
they had constructed themselves and grepped source text, both green over 55 lines
of dead code. The second was instance 1's summary literal, under which a tenth
criterion ran, passed, and was invisible in a total that said nine. **A criterion
that does not execute the path it claims to cover is not a criterion**, and none
of the three was found by the suite. All three were found by reading output.

**9. The cache read that was checked against the wrong property.** Stage 11's
verifier loaded its binder records with `stage5.read_binders()`, and stage 10's
with the same call. That function validates the payload against its own recorded
digest and against nothing else — not the Stage 4 hash, not the gene set, not
the indication. Its four sibling verifiers all go through
`load_or_retrieve(decisions, source, manifest["stage4_hash"])`, which refuses a
cache belonging to another configuration.

The binder cache is a single shared slot. The last thing to write it in a full
run is the **deliberately degraded** run inside the multi-indication verifier —
the one that screens with an unresolvable dependency lineage to prove the
platform names its missing sources. Its pool differs from the real one by two
genes and by an order swap at ranks 8 and 9.

Run against that state, stage 11's verifier reported **6 of 6 criteria clear**
and output identical to the logged run. Nothing in it could see that its binders
came from a run whose entire purpose was that a source was missing. Only the
order of the suite kept it from mattering, and `run_all.py` reuses derived
artifacts unless `--fresh` is passed, so the next ordinary run reads it.

*Broken-thing question:* if the binder set belonged to a different screen
entirely, every criterion would report exactly what it reported.

**This is the worse shape, and it is worth separating from instance 8.** An
empty set at least prints its zero: K4 said *"0 parts in the first construct"*,
and the number was visible to anyone who read the line. A wrong-but-populated
set prints 107 scored sequences and a full attrition table. It produces
confident, well-formed, plausible output with no number anywhere in it that a
reader could challenge, because every number is internally consistent — they
were all computed from the same wrong input.

The digest check is what makes it feel safe, and it is a real check. It proves
the file is not truncated or corrupt. It says nothing whatever about whether the
file belongs to this run, and the two questions are easy to conflate because
both are answered by the word "valid". **A check that passes on the wrong
property reads exactly like a check that passed.** Stage 10's own D5 has the
same shape from the other side: *"107 rows against 107 binders carrying a
sequence"* compares the scored rows to the records they were scored from, which
any binder set satisfies.

The fix is not a better digest. It is that an artifact must be admitted by
something that knows what run it belongs to — the configuration hash the
producing stage recorded — and every reader of a shared slot must go through the
same door.

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
