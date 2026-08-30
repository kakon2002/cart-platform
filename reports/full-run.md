# Full run

Derived artifacts deleted and rebuilt; raw sources read from `data/` unchanged.

**118/126 criteria clear across 12 stages**, 32.4 minutes.

| Stage | | Criteria | | Time |
| --- | --- | --- | --- | --- |
| 1 | Design spec | 31/31 | clear | 1s |
| 2 | Surface proteome | 2/2 | clear | 1s |
| 3 | Target discovery | 12/13 | **TRIPPED** | 20s |
| 4 | Target pairing | 8/14 | **TRIPPED** | 9s |
| 4a | Architecture routing | 11/12 | **TRIPPED** | 344s |
| 5 | Binder discovery | 7/7 | clear | 269s |
| 6 | Construct assembly | 8/8 | clear | 1s |
| 9 | Safety gate | 7/7 | clear | 4s |
| 10 | Developability | 6/6 | clear | 1s |
| 11 | Final ranking | 6/6 | clear | 4s |
| API | HTTP surface | 10/10 | clear | 277s |
| MULTI | Multi-indication | 10/10 | clear | 1014s |

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
- clear `R3` — CEACAM5 composite=0.8769 measured_weight=0.55 tier_rank=100
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

### Stage 4 — Target pairing (8/14)

- clear `P1` — combined risk vs min of members, rho=0.7929 (limit 0.95)
- clear `P2` — 6,337 of 19,900 (31.84%) beat the better member by more than 0.05 (limit 1%)
- clear `P3` — 102 blocked targets rescued by some pair: ABCC3, ACSL5, ADAM12, ADAM19, ADAM9, ADGRF1, AMIGO2, AMN
- **TRIPPED** `P4` — f_AB vs f_A x f_B over 19,503 measured pairs, rho=0.9924 (limit 0.98)
- clear `P5` — 0 pairs marked cleared on the optimistic arm
- clear `P6` — 0 recommended pairs are unmeasured
- **TRIPPED** `P7` — 2 cleared pairs contain a ubiquitous immune protein (in pool: ['HLA-A']) ['HLA-A+LRRC15', 'NPSR1+HLA-A']
- **TRIPPED** `P8` — 142 of 272 cleared pairs (52.2%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
- clear `P10` — 0 targets recommended dual despite clearing alone
- clear `P11` — outcome spread {'DUAL': 100, 'UNRESOLVED': 96, 'SINGLE': 2, 'NO_DESIGN': 2}, largest 50.0% (limit 95%)
- **TRIPPED** `P12` — 73 of 100 dual recommendations change at 2 counts (73.0%, limit 50%)
- **TRIPPED** `P13` — most common partner takes 72.0% of dual recommendations (LRRC15)
- clear `P14` — top pair is NPSR1+PTPRN2, top two singles are TMC5+ITGB6
- **TRIPPED** `P15` — pool halved to 100: 51 of 53 shared dual targets change partner (96.2%, limit 50%)

### Stage 4a — Architecture routing (11/12)

- clear `A1` — architecture is order-independent; 0 differ on a reversed pool
- clear `A2` — 0 targets routed ADAPTOR that CONVENTIONAL would have admitted
- clear `A3` — 0 routed risks differ from the Stage 3 risk
- clear `A4` — 0 CONVENTIONAL targets sit above the persistent ceiling
- clear `A5` — NPSR1 (risk 0.0277) routes CONVENTIONAL
- **TRIPPED** `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
- clear `A7` — 0 targets resolve to NO_ARCHITECTURE with no reason
- clear `A8` — with no declared terminable ceiling, 0 targets still route ADAPTOR
- clear `A9` — adaptor admissions across the ceiling sweep: 0.15->0, 0.2->1, 0.25->1, 0.3->4, 0.35->7, 0.4->22, 0.5->61, 0.6->103, 0.7->145
- clear `A10` — declared ceiling 0.35 matches the spec value 0.35
- clear `A11` — 1 adaptor constructs; 0 emit a sequence despite an unsupplied binder
- clear `A12` — 1 of 1 adaptor constructs carry a structure-derived binder the origin check can see; epitope immunogenicity stays NOT_CONNECTED on all of them, so the species gap and the immunogenicity gap remain separate

### Stage 5 — Binder discovery (7/7)

- clear `B1` — 77 targets have entries but no antibody among them (an inert stage would have none); 4 of 130 with entries echo their count
- clear `B3` — known answers hold: 5 sequence targets, 2 pinned structure targets, 2 documented negatives, and CEACAM5 at 1 structure entries of 6 (expected 1)
- clear `B8` — 5 of 5 known targets return a binder on some route
- clear `B5` — 0 candidates carry an affinity value (the source does not have one)
- clear `B7` — pool order carried through unchanged
- clear `B13` — 200 records out of 200 decisions in; 0 genes dropped, 0 added
- clear `B11` — 0 candidates claim a resolved isoform (neither route can determine one)

### Stage 6 — Construct assembly (8/8)

- clear `K1` — 12 constructs translate back to their own sequence
- clear `K2` — 12 constructs carry their binders verbatim, including the pinned MUC16, MUC17
- clear `K3` — domain boundaries partition every construct exactly
- clear `K4` — every part of every construct names its source (19 parts in the first construct)
- clear `K5` — part costs sum to the printed total for every construct
- clear `K6` — every buildable construct carries the mandatory safety switch
- clear `K7` — every owed construct was built and none was built without a binder; 16 targets have a binder but no recommendation (§5.1)
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
- clear `D5` — 112 rows against 112 binders carrying a sequence
- clear `D6` — no liability is summed into a single score; flags are counted and listed

### Stage 11 — Final ranking (6/6)

- clear `N1` — a dominated point is excluded from the front
- clear `N2` — both non-dominated points are on the front
- clear `N3` — attrition accounts for 200 of 200
- clear `N4` — no weighted or summed score across objectives is emitted
- clear `N5` — status RANKED matches the survivor count 1
- clear `N6` — 200 rows against the 200 the Stage 4 manifest records

### Stage API — HTTP surface (10/10)

- clear `A1` — project created (201), target_antigen None and discovery mode B
- clear `A2` — a view before any run answers 409 RUN_NOT_COMPLETE with instructions, not an empty list
- clear `A3` — a run returns 202 with job c4a33892d699 rather than blocking
- clear `A4` — job finished complete after stages ['sources', 'pairing', 'binders', 'ranking']
- clear `A5` — 200 BUILDABLE: 2 buildable = 2 complete + 0 awaiting a binder; 11 over budget, 5 reasons
- clear `A6` — end state RANKED, attrition accounts for 198 + 2 of 200; 2 reached = 2 complete + 0 awaiting
- clear `A7` — top target CEACAM5 with a 6-component breakdown
- clear `A8` — pairs carry the span percentile beside the raw fraction (0.006321856890514115 at percentile 0.1291)
- clear `A10` — an unknown project answers 404 NOT_FOUND and one that exists without a finished run answers 409 RUN_NOT_COMPLETE: a client can tell a bad id from a run in progress
- clear `A9` — evidence trail for MSLN spans 7 stages: stage3, stage4, stage5, stage6, stage9, stage10, stage11

### Stage MULTI — Multi-indication (10/10)

- clear `M1` — 44 indication-tagged artifacts; 0 changed, 0 disappeared after running both
- clear `M2` — shared sources carry no per-indication copy (0 found)
- clear `M3` — indication-specific module constants remaining: none
- clear `M4` — an atlas-less indication returns NOT_USABLE with no ranking
- clear `M5` — the refusal names malignant_vs_stroma as the missing discriminator, not just a lost weight
- clear `M6` — Mode A on CD19 returns NOT_ASSESSED with 3 reasons
- clear `M7` — CD19 ranks 1303 of 3400 -- outside the top 20, so the verdict is not self-agreement
- clear `M8` — Mode A and Mode B report the same evidence for CD19 (risk 0.5353, composite 0.227)
- clear `M9` — reference unchanged: top3 ['CEACAM5', 'TMC5', 'MUCL3'], pool 200, hash a91c696f2e1318f7, outcomes {'DUAL': 100, 'UNRESOLVED': 95, 'SINGLE': 2, 'ADAPTOR': 1, 'NO_DESIGN': 2}
- clear `M10` — a degraded indication names its missing source: ["dependency lineage 'NoSuchLineage': ValueError: need at least one array to concatenate"]

## What the platform returns for this indication

```
    GET /constructs -> BUILDABLE
      - Every surviving design routes to an adaptor architecture. That is two manufactured biologics, not one: the receptor and, separately, the tagged adaptor antibody that gives it its specificity. The second carries its own CMC package and its own regulatory path, and the payload budget the adaptor route saves is paid there instead.
      - The anti-tag binder is a murine anti-GCN4 single-chain Fv retrieved from PDB 1P4B entities 1+2 at revision 1.4. PDB does not record the identifier 52SR4 anywhere in that entry. The identification rests on an exact match between the deposited CDRs and those quoted for 52SR4 in the Calibr/Scripps patent family, together with the shared Zahnd 2004 primary citation. That is an inference drawn by this pipeline, not a fact taken from the source, and a reader is entitled to disagree with it.
      - The retrieved binder is murine. The clinical construct in this tag system is humanized and its sequence is not established, so what is built here is the crystallised murine scFv, not the clinical one. Non-human sequence content is an explicit Stage 9 immunogenicity question and that arm is empty: epitope-level immunogenicity reports NOT_CONNECTED on every row because no epitope source is connected, and the origin check reads INN name stems, which a structure-derived binder does not carry. Nothing in this pipeline has assessed the immunogenicity of this binder.
      - The binder is emitted as deposited, including its crystallisation artifacts, because trimming them is a design decision this pipeline does not take silently. Each construct therefore carries MADYADA at residues 22-28, expression leader carried on the light-chain entity; and ASGADHHHHHH at residues 270-280, purification tag carried on the heavy-chain entity. As emitted these are not manufacturable: the first is a second leader sitting inside the mature protein, the second a His tag between the binder and the hinge. Removing them is a wet-lab step that has not been taken here.
      - 1 advanced design(s) are available, all of them adaptor receptors, which is the architecture row the spec lists for serious normal-tissue expression.
    GET /result     -> RANKED
      blocked on normal tissue risk      - 197     3 remain
      no design recommended              -   0     3 remain
      no binder retrieved                -   1     2 remain
      no construct assembled             -   0     2 remain
      construct over budget              -   0     2 remain
```

## Tripped

- Stage 3 `R14` — no criticality tier places two staining levels on opposite sides of the 0.15 ceiling, so the arm gates on presence only — tier1 0.288/0.379/0.460(all above)  tier2 0.173/0.227/0.276(all above)  tier3 0.086/0.114/0.138(all below)
  - **Not on the accepted list.** Open decision, not a regression.
  - The staining arm vetoes on presence rather than amount: every grade blocks in tiers 1 and 2, none reaches the ceiling in tier 3. Under a conservative tolerance that may be the right design. It is a decision about tolerance, priced in reports/staining-veto-decision.md, and it is not on the accepted list because nobody has taken it yet. The non-zero exit is this decision being open, not a broken build.
- Stage 4 `P4` — f_AB vs f_A x f_B over 19,503 measured pairs, rho=0.9924 (limit 0.98)
  - Accepted: Coverage is span-confounded: f_AB tracks genomic span (+0.68) more than expression (+0.20). Reported beside a span-matched percentile and removed from partner selection.
- Stage 4 `P7` — 2 cleared pairs contain a ubiquitous immune protein (in pool: ['HLA-A']) ['HLA-A+LRRC15', 'NPSR1+HLA-A']
  - Accepted: One cleared pair contains HLA-A. Recorded rather than filtered, because the pool is not curated by hand.
- Stage 4 `P8` — 142 of 272 cleared pairs (52.2%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
  - Accepted: 48.6% of cleared pairs stop clearing if an unmeasured antigen saturates its organ. This is the cost of treating missing as a third state instead of imputing it.
- Stage 4 `P12` — 73 of 100 dual recommendations change at 2 counts (73.0%, limit 50%)
  - **Not on the accepted list. This is new.**
- Stage 4 `P13` — most common partner takes 72.0% of dual recommendations (LRRC15)
  - **Not on the accepted list. This is new.**
- Stage 4 `P15` — pool halved to 100: 51 of 53 shared dual targets change partner (96.2%, limit 50%)
  - Accepted: Partner choice is unstable under pool halving (71.4%). The pairing stage is complete-with-limitations by decision.
- Stage 4a `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
  - Accepted: A positive pin written before the run. It expected MSLN to route to an adaptor because it matches that row's condition in words: serious normal-tissue expression. It does not, because its measured risk 0.6366 is nearly twice the declared terminable ceiling of 0.35. Admitting it needs a ceiling near 0.65, which also admits about 120 others - a clinical policy decision, not a code change. The ceiling stays where the spec pinned it and A9 reports the whole sweep so the trade is visible.
