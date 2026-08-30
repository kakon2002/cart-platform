# Full run

Derived artifacts deleted and rebuilt; raw sources read from `data/` unchanged.

**119/125 criteria clear across 12 stages**, 34.9 minutes.

| Stage | | Criteria | | Time |
| --- | --- | --- | --- | --- |
| 1 | Design spec | 31/31 | clear | 1s |
| 2 | Surface proteome | 2/2 | clear | 1s |
| 3 | Target discovery | 12/13 | **TRIPPED** | 23s |
| 4 | Target pairing | 10/14 | **TRIPPED** | 11s |
| 4a | Architecture routing | 10/11 | **TRIPPED** | 366s |
| 5 | Binder discovery | 7/7 | clear | 298s |
| 6 | Construct assembly | 8/8 | clear | 0s |
| 9 | Safety gate | 7/7 | clear | 4s |
| 10 | Developability | 6/6 | clear | 1s |
| 11 | Final ranking | 6/6 | clear | 4s |
| API | HTTP surface | 10/10 | clear | 306s |
| MULTI | Multi-indication | 10/10 | clear | 1080s |

## Every criterion

### Stage 1 — Design spec (31/31)

- clear `check 1` — discovery_mode — got 'B', expected 'B'
- clear `check 2` — target_antigen — got None, expected None
- clear `check 3` — cancer_type — got 'Pancreatic Ductal Adenocarcinoma', expected 'Pancreatic Ductal Adenocarcinoma'
- clear `check 4` — malignancy_type — got 'solid', expected 'solid'
- clear `check 5` — product_type — got 'autologous', expected 'autologous'
- clear `check 6` — car_format — got 'auto', expected 'auto'
- clear `check 7` — safety_tolerance — got 'conservative', expected 'conservative'
- clear `check 8` — vector_payload_limit_kb — got 4.7, expected 4.7
- clear `check 9` — max_genetic_edits — got 2, expected 2
- clear `check 10` — pancreas override tier — got 2, expected 2
- clear `check 11` — mistyped field rejected — got True, expected True
- clear `check 12` — blank antigen -> None — got None, expected None
- clear `check 13` — blank antigen -> mode B — got <DiscoveryMode.DISCOVER: 'B'>, expected <DiscoveryMode.DISCOVER: 'B'>
- clear `check 14` — supplied antigen -> mode A — got 'A', expected 'A'
- clear `check 15` — blank cancer_type rejected — got True, expected True
- clear `check 16` — override without rationale rejected — got True, expected True
- clear `check 17` — spec discovery_mode — got 'B', expected 'B'
- clear `check 18` — datasets — got 10, expected 10
- clear `check 19` — blocking datasets — got 8, expected 8
- clear `check 20` — construct budget kb — got 3.5, expected 3.5
- clear `check 21` — safety switch required — got True, expected True
- clear `check 22` — risk ceiling — got 0.15, expected 0.15
- clear `check 23` — allowed formats — got 5, expected 5
- clear `check 24` — auto excluded — got True, expected True
- clear `check 25` — spec target_antigen — got None, expected None
- clear `check 26` — unresolved availability score — got 0.0, expected 0.0
- clear `check 27` — resolved availability score — got 0.75, expected 0.75
- clear `check 28` — validation datasets — got 7, expected 7
- clear `check 29` — validation blocking — got 5, expected 5
- clear `check 30` — input not mutated by build — got None, expected None
- clear `check 31` — project id unique — got True, expected True

### Stage 2 — Surface proteome (2/2)

- clear `validation sets` — pass
- clear `count drift` — yes

### Stage 3 — Target discovery (12/13)

- clear `R1` — outside top decile: none
- clear `R2` — cleared the ceiling: none (across 4 accessions)
- clear `R3` — CEACAM5 composite=0.8769 measured_weight=0.55 tier_rank=1
- clear `R4` — non-surface entries present: 0
- clear `R5` — highest |rho| tumour_vs_normal=0.815 (over measured targets only)
- clear `R6` — worst retention 94% at malignant_expression x0.8
- clear `R7` — cleared 646 of 3,466
- clear `R8` — most repeated composite 0.1273 occurs 17x (0.50%)
- clear `R9` — best unresolved None vs best protein-confirmed 0.8769
- clear `R10` — 0 of 100 reached only after a symbol failed
- clear `R11` — 22 of 25 depart from the systematic offset (3.6x) by more than 2x; all listed, so none unread
- clear `R12` — worst retention 82% at c3_fold x0.5
- **TRIPPED** `R14` — no criticality tier places two staining levels on opposite sides of the 0.15 ceiling, so the arm gates on presence only — tier1 0.288/0.379/0.460(all above)  tier2 0.173/0.227/0.276(all above)  tier3 0.086/0.114/0.138(all below)

### Stage 4 — Target pairing (10/14)

- clear `P1` — combined risk vs min of members, rho=0.7926 (limit 0.95)
- clear `P2` — 6,328 of 19,900 (31.80%) beat the better member by more than 0.05 (limit 1%)
- clear `P3` — 99 blocked targets rescued by some pair: ACSL5, ADAM9, ADGRF1, AMIGO2, AMN, ANO1, ANTXR1, AQP5
- **TRIPPED** `P4` — f_AB vs f_A x f_B over 14,535 measured pairs, rho=0.9934 (limit 0.98)
- clear `P5` — 0 pairs marked cleared on the optimistic arm
- clear `P6` — 0 recommended pairs are unmeasured
- **TRIPPED** `P7` — 1 cleared pairs contain a ubiquitous immune protein (in pool: ['HLA-A']) ['NPSR1+HLA-A']
- **TRIPPED** `P8` — 108 of 222 cleared pairs (48.6%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
- clear `P10` — 0 targets recommended dual despite clearing alone
- clear `P11` — outcome spread {'NO_DESIGN': 185, 'DUAL': 13, 'UNRESOLVED': 1, 'SINGLE': 1}, largest 92.5% (limit 95%)
- clear `P12` — 4 of 13 dual recommendations change at 2 counts (30.8%, limit 50%)
- clear `P13` — most common partner takes 38.5% of dual recommendations (LAMP5)
- clear `P14` — top pair is NPSR1+PTPRN2, top two singles are CEACAM5+TMC5
- **TRIPPED** `P15` — pool halved to 100: 5 of 7 shared dual targets change partner (71.4%, limit 50%)

### Stage 4a — Architecture routing (10/11)

- clear `A1` — architecture is order-independent; 0 differ on a reversed pool
- clear `A2` — 0 targets routed ADAPTOR that CONVENTIONAL would have admitted
- clear `A3` — 0 routed risks differ from the Stage 3 risk
- clear `A4` — 0 CONVENTIONAL targets sit above the persistent ceiling
- clear `A5` — NPSR1 (risk 0.0277) routes CONVENTIONAL
- **TRIPPED** `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
- clear `A7` — 0 targets resolve to NO_ARCHITECTURE with no reason
- clear `A8` — with no declared terminable ceiling, 0 targets still route ADAPTOR
- clear `A9` — adaptor admissions across the ceiling sweep: 0.15->0, 0.2->0, 0.25->1, 0.3->5, 0.35->9, 0.4->26, 0.5->65, 0.6->109, 0.7->150
- clear `A10` — declared ceiling 0.35 matches the spec value 0.35
- clear `A11` — 8 adaptor constructs; 0 emit a sequence despite an unsupplied binder

### Stage 5 — Binder discovery (7/7)

- clear `B1` — 73 targets have entries but no antibody among them (an inert stage would have none); 6 of 127 with entries echo their count
- clear `B3` — known answers hold: 5 sequence targets, 2 pinned structure targets, 2 documented negatives, and CEACAM5 at 1 structure entries of 6 (expected 1)
- clear `B8` — 5 of 5 known targets return a binder on some route
- clear `B5` — 0 candidates carry an affinity value (the source does not have one)
- clear `B7` — pool order carried through unchanged
- clear `B13` — 200 records out of 200 decisions in; 0 genes dropped, 0 added
- clear `B11` — 0 candidates claim a resolved isoform (neither route can determine one)

### Stage 6 — Construct assembly (8/8)

- clear `K1` — 2 constructs translate back to their own sequence
- clear `K2` — 2 constructs carry their binders verbatim, including the pinned MUC16, MUC17
- clear `K3` — domain boundaries partition every construct exactly
- clear `K4` — every part of every construct names its source (19 parts in the first construct)
- clear `K5` — part costs sum to the printed total for every construct
- clear `K6` — every buildable construct carries the mandatory safety switch
- clear `K7` — every owed construct was built and none was built without a binder; 26 targets have a binder but no recommendation (§5.1)
- clear `K8` — 200 rows and 200 distinct genes against the 200 the Stage 4 manifest records

### Stage 9 — Safety gate (7/7)

- clear `S1` — registry returns trials for both pins (MSLN=169, CLDN18=3)
- clear `S2` — pinned names classify as expected (Amatuximab=chimeric, Zolbetuximab=chimeric)
- clear `S3` — every target with a binder carries a Stage 3 risk
- clear `S4` — no target over the ceiling escapes BLOCKED
- clear `S5` — epitope immunogenicity is NOT_CONNECTED on every row
- clear `S6` — 200 rows and 200 distinct genes against the 200 the Stage 4 manifest records
- clear `S7` — every target with stopped trials is flagged or blocked

### Stage 10 — Developability (6/6)

- clear `D1` — poly-K pI 11.463, poly-E pI 2.623 — the charge model orders the two known answers correctly
- clear `D2` — NST yields 1 sequon, NPT yields 0 — the proline exclusion holds
- clear `D3` — a 3-cysteine control reports parity 'odd', and never 'unpaired: 0'
- clear `D4` — every scored binder has a pI in 1..14 and 0..5 flags
- clear `D5` — 113 rows against 113 binders carrying a sequence
- clear `D6` — no liability is summed into a single score; flags are counted and listed

### Stage 11 — Final ranking (6/6)

- clear `N1` — a dominated point is excluded from the front
- clear `N2` — both non-dominated points are on the front
- clear `N3` — attrition accounts for 200 of 200
- clear `N4` — no weighted or summed score across objectives is emitted
- clear `N5` — status NO_DESIGN_REACHES_THE_END matches the survivor count 0
- clear `N6` — 200 rows against the 200 the Stage 4 manifest records

### Stage API — HTTP surface (10/10)

- clear `A1` — project created (201), target_antigen None and discovery mode B
- clear `A2` — a view before any run answers 409 RUN_NOT_COMPLETE with instructions, not an empty list
- clear `A3` — a run returns 202 with job 21a6b0af5c61 rather than blocking
- clear `A4` — job finished complete after stages ['sources', 'pairing', 'binders', 'ranking']
- clear `A5` — 200 BUILDABLE_AWAITING_BINDER: 8 buildable = 0 complete + 8 awaiting a binder; 2 over budget, 1 reasons
- clear `A6` — end state RANKED_AWAITING_BINDER, attrition accounts for 192 + 8 of 200; 8 reached = 0 complete + 8 awaiting
- clear `A7` — top target CEACAM5 with a 6-component breakdown
- clear `A8` — pairs carry the span percentile beside the raw fraction (0.006321856890514115 at percentile 0.0785)
- clear `A10` — an unknown project answers 404 NOT_FOUND and one that exists without a finished run answers 409 RUN_NOT_COMPLETE: a client can tell a bad id from a run in progress
- clear `A9` — evidence trail for MSLN spans 7 stages: stage3, stage4, stage5, stage6, stage9, stage10, stage11

### Stage MULTI — Multi-indication (10/10)

- clear `M1` — 38 indication-tagged artifacts; 0 changed, 0 disappeared after running both
- clear `M2` — shared sources carry no per-indication copy (0 found)
- clear `M3` — indication-specific module constants remaining: none
- clear `M4` — an atlas-less indication returns NOT_USABLE with no ranking
- clear `M5` — the refusal names malignant_vs_stroma as the missing discriminator, not just a lost weight
- clear `M6` — Mode A on CD19 returns NOT_ASSESSED with 3 reasons
- clear `M7` — CD19 ranks 1303 of 3400 -- outside the top 20, so the verdict is not self-agreement
- clear `M8` — Mode A and Mode B report the same evidence for CD19 (risk 0.5353, composite 0.227)
- clear `M9` — reference unchanged: top3 ['CEACAM5', 'TMC5', 'MUCL3'], pool 200, hash a91c696f2e1318f7, outcomes {'NO_DESIGN': 177, 'DUAL': 13, 'ADAPTOR': 8, 'UNRESOLVED': 1, 'SINGLE': 1}
- clear `M10` — a degraded indication names its missing source: ["dependency lineage 'NoSuchLineage': ValueError: need at least one array to concatenate"]

## What the platform returns for this indication

```
    GET /constructs -> BUILDABLE_AWAITING_BINDER
      - 8 design(s) fit the 3500 bp budget but carry no binder sequence: the adaptor receptor binds a tag, and no anti-tag binder exists in the connected sources, so its size is declared and its sequence is not invented.
    GET /result     -> RANKED_AWAITING_BINDER
      blocked on normal tissue risk      - 191     9 remain
      no design recommended              -   0     9 remain
      no binder retrieved                -   1     8 remain
      no construct assembled             -   0     8 remain
      construct over budget              -   0     8 remain
```

## Tripped

- Stage 3 `R14` — no criticality tier places two staining levels on opposite sides of the 0.15 ceiling, so the arm gates on presence only — tier1 0.288/0.379/0.460(all above)  tier2 0.173/0.227/0.276(all above)  tier3 0.086/0.114/0.138(all below)
  - **Not on the accepted list.** Open decision, not a regression.
  - The staining arm vetoes on presence rather than amount: every grade blocks in tiers 1 and 2, none reaches the ceiling in tier 3. Under a conservative tolerance that may be the right design. It is a decision about tolerance, priced in reports/staining-veto-decision.md, and it is not on the accepted list because nobody has taken it yet. The non-zero exit is this decision being open, not a broken build.
- Stage 4 `P4` — f_AB vs f_A x f_B over 14,535 measured pairs, rho=0.9934 (limit 0.98)
  - Accepted: Coverage is span-confounded: f_AB tracks genomic span (+0.68) more than expression (+0.20). Reported beside a span-matched percentile and removed from partner selection.
- Stage 4 `P7` — 1 cleared pairs contain a ubiquitous immune protein (in pool: ['HLA-A']) ['NPSR1+HLA-A']
  - Accepted: One cleared pair contains HLA-A. Recorded rather than filtered, because the pool is not curated by hand.
- Stage 4 `P8` — 108 of 222 cleared pairs (48.6%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
  - Accepted: 48.6% of cleared pairs stop clearing if an unmeasured antigen saturates its organ. This is the cost of treating missing as a third state instead of imputing it.
- Stage 4 `P15` — pool halved to 100: 5 of 7 shared dual targets change partner (71.4%, limit 50%)
  - Accepted: Partner choice is unstable under pool halving (71.4%). The pairing stage is complete-with-limitations by decision.
- Stage 4a `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
  - Accepted: A positive pin written before the run. It expected MSLN to route to an adaptor because it matches that row's condition in words: serious normal-tissue expression. It does not, because its measured risk 0.6366 is nearly twice the declared terminable ceiling of 0.35. Admitting it needs a ceiling near 0.65, which also admits about 120 others - a clinical policy decision, not a code change. The ceiling stays where the spec pinned it and A9 reports the whole sweep so the trade is visible.
