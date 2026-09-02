# The construct verifiers read a routed decision set

Written before any implementation. It follows the amendment at the end of
`specs/stage6-construct-assembly.md`, which recorded the problem and declined to
fix it because the fix is a scope decision rather than a repair. This document
takes that decision, states what it costs, and fixes the criteria before the run
that tests them.

Status: **implemented**. §1–§8 were written and committed before any code
changed; §9 records what the run then reported, against the numbers §2 and §6
predicted.

---

## 1. What the construct verifiers read, and why it cannot exercise them

`verify_construct.py:39` reads the persisted Stage 4 artifact:

```python
decisions, manifest = stage4.read_decisions(allow_unusable=True)
```

That artifact is written once, by `verify_pairing.py:390`, from the decisions
produced at `verify_pairing.py:163`:

```python
decisions = stage4.decide(pool, pairs, tumour_tpm, per_organ=per_organ)
```

`tolerances` is not passed. `stage4.decide` treats that as a positive
instruction — decision 32 in `specs/design-decisions.md` — and disables routing
rather than inventing a ceiling. Every row therefore carries
`architecture="NOT_CONFIGURED"` and
`route_reason="no tolerances supplied; routing disabled"`.

Three consequences follow, and only the first is visible in the log.

**No adaptor row can exist.** `stage4.decide` emits `ADAPTOR` only from the
branch guarded by `decided is not None and decided.architecture == routing.ADAPTOR`.
With routing off, `decided` is `None`, so the branch is unreachable. The adaptor
is the architecture the platform returns for **every** surviving design in the
worked indication.

**Stage 6 can build almost nothing.** Only `SINGLE`, `DUAL` and `ADAPTOR` are
build-eligible. Of what remains, the 3 `SINGLE` rows retrieve no binder and all
30 `DUAL` rows have a partner that carries none. Zero constructs are assembled.

**Five criteria then report success over nothing.** K1, K3, K4, K5 and K6 all
iterate the assembled set. On an empty set every loop body is skipped, every
failure list is empty, and every criterion clears. K4 prints *"0 parts in the
first construct"* and clears on that sentence.

## 2. The measurement

Both decision sets were built from the same pool, the same pairs and the same
binder records, differing only in whether `tolerances` was passed.

| | routing off | routing on |
| --- | --- | --- |
| rows | 200 | 200 |
| `NO_DESIGN` | 167 | 162 |
| `DUAL` | 30 | 30 |
| `SINGLE` | 3 | 3 |
| `ADAPTOR` | 0 | **5** |
| constructs assembled | **0** | **5** |

The five rows that move are FER1L6, GPR35, TMEM92, TNFSF9 and BTNL8, each
`NO_DESIGN` → `ADAPTOR`. No `DUAL` and no `SINGLE` row changes. The routing-off
set reproduces the persisted artifact field for field, so this is a measurement
of the change and not of a drifting pool.

Each of the five assembles at 2,868 bp against the 3,500 bp budget, `BUILDABLE`,
10 parts — the same five the architecture-routing verifier already reports, and
the same five the API already serves.

### What the current criteria say once there is something to look at

The existing K1–K8 logic, run unchanged against the routed set:

| | routing off | routing on | |
| --- | --- | --- | --- |
| K1 | clear over **0** | clear over **5** | now executes |
| K2 | TRIPPED, nothing to pin on | TRIPPED, **5 false failures** | see §4 |
| K3 | clear over **0** | clear over **5** | now executes |
| K4 | clear, *"0 parts"* | clear, **10 parts** | now executes |
| K5 | clear over **0** | clear over **5** | now executes |
| K6 | clear over **0** | clear over **5** | now executes |
| K7 | clear | **TRIPPED, 5 false failures** | see §4 |
| K8 | clear, 200 rows | clear, 200 rows | unchanged |

The five failures under K2 read *"chosen binder 'anti-tag binder, peptide
neo-epitope, …' not in Stage 5"*, and the five under K7 read *"construct without
a usable binder"*. Both are the same mistake: **both criteria assume every
construct's binder is a Stage 5 sequence candidate, and an adaptor construct's
binder is not.** It is the anti-tag part, retrieved from a deposited structure.
The constructs are correct; the criteria cannot see the route that supplied
them. Fixing that is §4, and it is a restatement, not a relaxation — each gains
a positive assertion it does not have today.

## 3. Where the routed set comes from

Two ways to give the construct verifier a routed set. They are not equivalent.

**A. `verify_pairing` writes the set the platform actually uses**, by passing
the tolerances it already has in hand — it builds `spec` at
`verify_pairing.py:63` and needs only
`spec.design_constraints.terminable_risk_ceiling` beside the persistent ceiling
it already reads.

**B. `verify_construct` derives a routed set for itself**, either by re-running
the pipeline or by re-routing the persisted rows.

**A is chosen.** B makes the verifier construct its own input, which is the
failure recorded as instances M4 and M5 — a check asserting a property of an
object it built. It would also leave two different objects in the suite both
called "the Stage 4 decisions", which is worse than one wrong one.

The deeper reason is that the artifact is mislabelled today. It is written to
`data/stage4/decisions.json`, read by five verifiers as *the* Stage 4 result,
and it is not the result the platform produces. A is the only option that makes
the file's name true.

### What A costs, measured

Routing the shared artifact changes what stages 5, 6, 9, 10 and 11 verify. The
whole of that change was measured before this document was written:

| stage | effect |
| --- | --- |
| 5 binders | 5 records change `outcome` only. No B criterion moves. |
| 6 constructs | 0 → 5 assembled. K1/K3/K4/K5/K6 begin to execute. K2 and K7 need §4. |
| 9 safety | 197 `BLOCKED` → 192 `BLOCKED` + 5 `FLAGGED`. **S4 trips.** |
| 10 developability | nothing; it scores binder sequences and reads no decision. |
| 11 ranking | `NO_DESIGN_REACHES_THE_END` → `RANKED`, 5 survivors. N3, N5 and N6 stay clear. |

Stage 11's new numbers are the API's existing numbers: 192 blocked, 3 without a
binder, 5 reaching the end. The suite currently reports a different result from
the platform it verifies, and this is what closes that gap.

**S4 tripping is correct and is not to be silenced.** S4 reads

```python
g.risk is not None and g.risk > ceiling and g.verdict != stage9.BLOCKED
```

against the single persistent ceiling, in a design that has two and holds them
apart on purpose: persistent 0.15 for exposure that cannot be withdrawn,
terminable 0.35 for exposure that can. The five adaptor targets sit above 0.15
and below 0.35 and are admitted against the terminable ceiling, which is the
whole point of the row. S4 has never seen a target admitted that way, because
until now no such target existed in the set it reads. It must be restated
against **the ceiling that was applied to each row**, with a second clause
asserting that a row may only be admitted above the persistent ceiling when its
own route declares the exposure terminable. That is strictly more than S4 says
today, and it is written in §5.

## 4. What each criterion tests once the set is non-empty

K1, K3, K5, K6 and K8 are unchanged in substance. K2 and K7 are restated. K0 is
new.

**K0 — new. The set is routed, and something was built.** Trips if the manifest
records no routing configuration, if any row carries the routing-disabled
reason, or if no construct assembled. Stated so that an empty or unrouted set
**fails**; it cannot be satisfied by absence. This is the criterion whose want
let five others clear on nothing.

**K1 — the DNA round trip.** Every assembled construct's DNA must translate back
to its own amino acid sequence. Trips also when nothing was assembled.

**K2 — the binder arrives verbatim, by the route that supplied it.** Three
clauses, one per route, and a construct is checked under the route its
architecture uses:

- *Stage 5 route* (`SINGLE`, `DUAL`): the chosen candidate must exist in the
  Stage 5 record under the name the construct claims, and its VH and VL must
  both appear verbatim in the sequence. For a dual, the partner's VH **and VL**
  as well — today only the partner's VH is checked, which is a gap this closes.
- *Anti-tag route* (`ADAPTOR`): the retrieved anti-tag part's sequence must
  appear verbatim; the segment carrying it must declare provenance `structure`
  with the accession the retrieval returned; and the construct's `binder_name`
  must be that part's name. A construct whose binder is not attributable to
  either route fails.
- *The two-arm join*: pins derived from the run through `two_armed_duals()`.
  Where that set is empty, K2 trips and says the join is untested, as the
  amendment already has it.

**K3 — the boundaries partition the sequence.** Unchanged. Trips also when
nothing was assembled.

**K4 — every part names its source.** Unchanged. Trips also when nothing was
assembled. The routed run gives it 10 parts to look at.

**K5 — the printed arithmetic agrees.** Unchanged. Trips also when nothing was
assembled.

**K6 — a buildable construct carries the safety switch.** Unchanged. Trips also
when nothing was assembled.

**K7 — owed and not owed, per architecture.** A row is owed a construct when its
outcome is build-eligible and every binder its architecture needs is available:
its own Stage 5 binder for a single, both arms' for a dual, and a supplied
anti-tag binder for an adaptor. A construct must not exist without the binder
its own architecture needs — for an adaptor that is the anti-tag part, not a
Stage 5 record. §5.1 of the assembly spec is unaffected: a target with a binder
and no recommendation is still counted and reported, not tripped on.

**K8 — the row count and gene set.** Unchanged.

### The rule behind K0

**A criterion must not be able to clear on an empty population.** Each of K1,
K3, K4, K5 and K6 examines the assembled constructs; each therefore trips when
that set is empty and says so, rather than reporting an empty loop as success.
K0 states the same thing once at the top, so the failure is legible without
reading six detail lines. Both are wanted: K0 makes it obvious, the individual
clauses make it true.

## 5. Five defects found alongside, and what happens to each

None of these was known when the work started. The first three were found while
measuring, the last two by the review of the change itself. Each is a
consequence of the same root — an artifact, a manifest or a criterion describing
a configuration that is not the one in front of it. All five are fixed here,
because each is either created or made live by routing the artifact, and leaving
one behind would mean shipping a change that breaks something it can see.

**5.1 The Stage 4 hash does not cover the routing configuration.**
`write_decisions` computes `configuration_hash(stage3_hash, pool_genes)` with no
tolerances, while `pipeline.py:173` computes it with them. Decision 13 against
`routing.py:177` says the declared ceilings belong in that hash, *"otherwise two
runs with opposite verdicts are indistinguishable by hash and a cached result is
reused across tolerances"*. Today the omission is harmless only because the
artifact is unrouted. **The moment it is routed, a routed and an unrouted
artifact hash identically.** Fixed in this change: `write_decisions` takes the
tolerances and threads them through the hash, and records them in the manifest
so K0 has something structural to read rather than a reason string to match.

**5.2 The manifest's outcome tally has no `ADAPTOR` row.**
`write_decisions` tallies `SINGLE`, `DUAL`, `NO_DESIGN` and `UNRESOLVED` only.
On the routed set that reports 195 of 200 rows and silently drops five. Nothing
reads the field today, which is why it has been wrong without effect. Fixed by
tallying every outcome the stage can emit, so a new outcome cannot go missing
from the count again.

**5.3 Stage 11's verifier reads the binder cache with nothing blessing it.**
`verify_ranking_final.py:26` calls `stage5.read_binders()`, which validates the
payload against its own digest and nothing else — not the Stage 4 hash, not the
gene set. Its four sibling verifiers all go through
`load_or_retrieve(decisions, source, manifest["stage4_hash"])`, which refuses a
cache from another configuration.

This is live, not hypothetical. `data/stage5/binders.json` at the time of
writing holds the binder set from the **deliberately degraded** run inside
`verify_indications.py` — the one that screens with an unresolvable dependency
lineage to prove the platform names its missing sources. That run is the last
writer of the single shared cache slot. Its pool differs from the persisted
decision pool by two genes (`AQP3`, `CLSTN1` for `LRP5`, `SLC16A4`) and by an
order swap at ranks 8 and 9.

Run against that state, `verify_ranking_final.py` reports **6/6 criteria clear**
and output identical to the logged run. Nothing in it can see that its binders
came from a run whose whole purpose was that a source was missing. Only the
order of the suite — stage 11 at 05:03, the degraded run at 05:19 — kept it from
mattering, and `run_all.py` reuses derived artifacts unless `--fresh` is given,
so the next non-fresh run reads it.

`verify_developability.py:15` does the same, and its own D5 — *"107 rows against
107 binders carrying a sequence"* — compares the scored rows to the very records
it scored, so a wrong binder set is invisible there too.

Fixed by routing both reads through the same blessed path as their siblings.

**5.4 The safety gate is blind to a structure-derived binder in both verifiers.**
`stage9.gate` learns that a receptor carries a structure-derived binder only from
the `constructs` argument. `pipeline.py:188` passes it. `verify_safety.py:65` and
`verify_ranking_final.py:57` do not, so `structure_binders()` returns nothing and
an adaptor row falls through to `NO_GATE` with the reason *"no binder, so there
is nothing to gate"*.

That sentence is instance 6 of the verification-assumptions note, printed on a
design whose receptor carries a murine binder. It was recorded as fixed. It is
fixed in the platform and was never fixed in the two verifiers, and it has been
invisible because no adaptor row existed in the set they read. Routing that set
makes it visible in the same change, so it is corrected in the same change: both
verifiers pass the constructs they have already built. Criterion C7 pins the
result — 5 `FLAGGED`, not 5 more `NO_GATE`.

**5.5 The pairing report prints a Stage 4 hash that omits the routing.**
`_report_biology` in `verify_pairing.py` recomputes the hash for display with
`configuration_hash(s3_hash, pool_genes)` and no tolerances, so once the artifact
is routed the printed hash and the manifest's `stage4_hash` disagree — a reader
comparing them would conclude the file came from another configuration. It has
never printed, because it runs only when nothing tripped and pairing trips on
five criteria, so this is a latent divergence rather than a wrong number anyone
has read. It is 5.1 in a second place and is fixed with it.

## 6. Criteria for this change, fixed before it runs

Positive pins, and each names a specific observation so that absence cannot
satisfy it.

| id | criterion |
| --- | --- |
| **C1** | The persisted decision set carries a routing configuration in its manifest, and no row carries `route_reason="no tolerances supplied; routing disabled"`. |
| **C2** | The persisted set contains **exactly 5** `ADAPTOR` rows, and they are FER1L6, GPR35, TMEM92, TNFSF9 and BTNL8. A positive pin: the set must not merely be routed, it must route to the architecture that was measured. |
| **C3** | `SINGLE` stays at 3 and `DUAL` at 30. Routing must convert `NO_DESIGN` rows only; if it moves a single or a dual, the change did more than it claims. |
| **C4** | The construct stage assembles **exactly 5** constructs, each 2,868 bp, `BUILDABLE`, 10 parts. |
| **C5** | The Stage 4 configuration hash **changes** from `8cb155b103141a36`. A routed set that hashed the same as the unrouted one would be 5.1 unfixed. |
| **C6** | K0 trips when handed an unrouted set, and trips when handed a routed set that assembles nothing. Both halves tested directly, because a gate that cannot fail is not a gate. |
| **C7** | Stage 9 reports 192 `BLOCKED`, 5 `FLAGGED`, 3 `NO_GATE`, and the applied ceilings in play are exactly {0.15, 0.35}. |
| **C8** | Stage 11 reports `RANKED` with 5 survivors and attrition 192 / 0 / 3 / 0 / 0, matching what the API already returns for this indication. |

C2, C4, C7 and C8 are the pins that a dead or half-applied change would fail.
C6 is the one that tests the new criterion rather than the new code.

## 7. What this does not do

- **It does not touch `stage4.decide`, `routing.route` or `stage6.build`.** No
  scoring, no threshold, no architecture rule changes. The routed decisions
  measured here come from the existing code, called the way the platform calls
  it.
- **It does not re-pin K2's two-arm clause.** No dual in this pool carries a
  binder on both arms; routing does not change that, because routing converts
  `NO_DESIGN` rows only. **K2 is expected to keep tripping after this change**,
  now for the single honest reason rather than for that reason plus five
  spurious ones. That is a finding about the pool, priced in the amendment to
  the assembly spec, and it is not repaired here.
- **It does not resolve R14 or P17.** Both are open decisions with their own
  documents.
- **It does not change the terminable ceiling**, or what routes to it. A6 stays
  on the accepted list for the reason recorded there.
- **It does not claim the five adaptor designs are good.** It claims the suite
  should verify them rather than verify their absence. Every caveat the API
  already prints about them — the second biologic, the murine binder, the
  crystallisation artifacts, the unassessed immunogenicity — stands untouched.

## 8. Order of work

1. This document, reviewed, committed before any code.
2. `write_decisions` / `read_decisions`: tolerances into the hash and the
   manifest (5.1), every outcome into the tally (5.2).
3. `verify_pairing`: pass the declared tolerances to `decide`.
4. `verify_construct`: K0 added; K2 and K7 restated; K1, K3, K4, K5, K6 made
   unable to clear on an empty population.
5. `verify_safety`: S4 restated against the applied ceiling, with the second
   clause on terminable exposure.
6. `verify_ranking_final`: the binder read routed through the blessed path (5.3).
7. The pattern file gains this as an instance, alongside M4/M5 and the summary
   literal.
8. Full suite, both indications, with every count reported beside the count
   this document predicts.

---

## 9. What the run reported

Every criterion in §6 landed on the number written before it ran.

| | predicted | observed | |
| --- | --- | --- | --- |
| C1 | routing recorded, 0 rows disabled | routing recorded, **0 of 200** | clear |
| C2 | exactly 5 `ADAPTOR`, the named five | **5**: FER1L6, GPR35, TMEM92, TNFSF9, BTNL8 | clear |
| C3 | `SINGLE` 3, `DUAL` 30 | **3**, **30** | clear |
| C4 | 5 constructs, 2,868 bp, `BUILDABLE`, 10 parts | **5**, 2,868 bp, `BUILDABLE`, 10 parts (50 across all five) | clear |
| C5 | hash moves off `8cb155b103141a36` | **`5d097e05887e5b28`** | clear |
| C6 | K0 trips on unrouted, and on empty | trips on **both**, each naming its own clause | clear |
| C7 | 192 `BLOCKED`, 5 `FLAGGED`, 3 `NO_GATE`, ceilings {0.15, 0.35} | **192 / 5 / 3**, both ceilings in play | clear |
| C8 | `RANKED`, 5 survivors, attrition 192 / 0 / 3 / 0 / 0 | **`RANKED`**, 5 survivors, 192 / 0 / 3 / 0 / 0 | clear |

### Stage 6

**8 of 9 clear, against 7 of 8 before.** K0 is the new criterion. K1, K3, K4,
K5 and K6 now read five constructs where they read none: K4 reports *"10 parts
in the first construct, 50 across all 5"* where it reported *"0 parts"*. K7
clears — the five false failures are gone with the per-route restatement.

**K2 still trips**, with the five spurious failures gone and the honest reason
left: *"no dual carries a binder on both arms, so the two-arm join is not
exercised anywhere in this decision set (5 of 200 rows assembled, 5 of them by
the anti-tag route, whose binder is verified above)"*. That is the finding the
assembly spec's amendment already priced, and this change does not repair it.
Routing converts `NO_DESIGN` rows only, so it was never going to.

### C6 in full, because it is the criterion that tests the criterion

Handed a decision set with routing switched off but constructs still present,
K0 trips alone and names the routing clause. Handed a routed set that assembles
nothing, the stage reports **2 of 9 clear** — K0, K1, K2, K3, K4, K5 and K6 all
trip, each in its own words: *"no construct to translate"*, *"no construct to
partition"*, *"no construct, so no part was examined"*, *"no construct to
cost"*, *"no buildable construct to check the switch on"*. K7 and K8 clear, and
correctly: nothing was owed and 200 rows exist. They read the decision set, not
the assembled set.

Under the old criteria the same empty set reported **7 of 8 clear**. That
difference is the whole of this change.

A third case, from the review: K2's clause for *"adaptor rows exist but no
anti-tag sequence was retrieved"* was written inside the loop over assembled
constructs, where it could never run — `_assemble` returns an empty protein when
any part is unsupplied, so such a construct never reaches that loop. A criterion
that cannot execute its own path, in the change whose subject is criteria that
cannot execute their own path. It was moved out to the decision set, where it is
reachable, and tested by substituting the unsupplied part: the stage then reports
**2 of 9 clear**, K2 naming the cause — *"5 adaptor row(s), but no anti-tag
sequence was retrieved, so the anti-tag join cannot be verified on any of them"*
— while K7 correctly clears, because nothing was owed. That is not a hypothetical
state; it is the state this repository was in before the anti-tag binder was
retrieved.

### Stage 9

S4 clears in its restated form and now says which branch admitted what: *"3
admitted against the persistent 0.15, 5 against the terminable 0.35, each on a
route declaring the exposure stoppable"*. The gate's own report line moves from
*"reached the immunogenicity and trials questions 0 of 200"* to **5 of 200** —
the first run in which this gate's own logic decided anything, rather than every
row being settled by risk carried from Stage 3.

### Stage 4, 5, 10 and 11

Pairing reports 10 of 15 with the same five tripped as before — P4, P8, P13,
P15, P17 — so routing changed no pairing verdict, which is C3 stated another
way. Binder discovery 7 of 7 and developability 6 of 6, both unchanged. Final
ranking 6 of 6, now over `RANKED` with five survivors and a two-design Pareto
front, where it previously reported that no design reaches the end.

That last line is the point of the exercise. The service has been returning five
designs for this indication throughout; the suite has been verifying that there
were none.
