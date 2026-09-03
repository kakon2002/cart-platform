# Full run

Derived artifacts deleted and rebuilt; raw sources read from `data/` unchanged.

**144/152 criteria clear across 13 stages**, 22.5 minutes.

| Stage | | Criteria | | Time |
| --- | --- | --- | --- | --- |
| 1 | Design spec | 31/31 | clear | 1s |
| 2 | Surface proteome | 2/2 | clear | 1s |
| 3 | Target discovery | 26/27 | **TRIPPED** | 21s |
| 4 | Target pairing | 10/15 | **TRIPPED** | 9s |
| 4a | Architecture routing | 11/12 | **TRIPPED** | 279s |
| 5 | Binder discovery | 7/7 | clear | 208s |
| 6 | Construct assembly | 8/9 | **TRIPPED** | 0s |
| 9 | Safety gate | 7/7 | clear | 4s |
| 10 | Developability | 6/6 | clear | 1s |
| 11 | Final ranking | 6/6 | clear | 4s |
| 12 | Candidate package | 9/9 | clear | 8s |
| API | HTTP surface | 10/10 | clear | 11s |
| MULTI | Multi-indication | 11/11 | clear | 805s |

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

### Stage 3 — Target discovery (26/27)

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
- clear `G1` — 901 targets carry an unmeasured malignant-to-stromal ratio and 0 of them were rejected; the gate fired on no absent measurement
- clear `G2` — known targets surviving the gate: CEACAM6=193.4, MUC1=34.8, CLDN18=43.3, MSLN=19.7, CEACAM5=exempt
- clear `G3` — the gate rejects 1665 targets including LRRC15=yes, 11 MHC class II (['HLA-DOA', 'HLA-DOB', 'HLA-DPA1']), 8 immunoglobulin (['IGHA1', 'IGHA2', 'IGHD'])
- clear `G4` — 1665 rejected of 2565 with a measured ratio; the gate neither passes everything nor empties the population
- clear `G5` — 1665 rejected against 1665 with a measured ratio at or below 1.0; none carries an absence note
- clear `G6` — risk recomputed independently for all 1665 rejected targets matches the stored value on 1665; the gate changed no risk
- **TRIPPED** `R14` — no criticality tier places two staining levels on opposite sides of the 0.15 ceiling, so the arm gates on presence only — tier1 0.288/0.379/0.460(all above)  tier2 0.173/0.227/0.276(all above)  tier3 0.086/0.114/0.138(all below)
- clear `T1` — the attribution reproduces the reported risk for all 3,400 targets that carry one, to within 1e-12 before rounding
- clear `T2` — every reported organ attains the maximum it is credited with; 314 target(s) reach it on more than one organ and say so
- clear `T3` — every one of 65,077 attributed organ rows recomputes its own score from the measurement it names, and names the arm that won
- clear `T4` — 30,454 organ rows carry NOT_MEASURED on the protein arm, and every organ's staining presence agrees with the atlas entry read independently, so absence is never scored as stained-and-clean
- clear `T5` — both arms decide verdicts that carry a non-zero score: BASELINE=2,769, STAINING=856 winning organs, out of BASELINE=4,893, STAINING=856 counting the zero-scored ties every absent protein produces
- clear `T6` — all 3,400 targets scoring two or more organs report a margin that agrees with their own organ list
- clear `T7` — the attribution payload carries no confidence, evidence-class or measured-weight field; the two scores stay apart
- clear `T8` — no gene symbol appears as a literal anywhere in the attribution code path (79 string literals checked against 3,454 symbols)

### Stage 4 — Target pairing (10/15)

- clear `P1` — combined risk vs min of members, rho=0.7910 (limit 0.95)
- clear `P2` — 6,426 of 19,900 (32.29%) beat the better member by more than 0.05 (limit 1%)
- clear `P3` — 75 blocked targets rescued by some pair: ACHE, ACSL5, ADAM9, ADGRF4, ADGRG6, AMIGO2, AMN, BTNL8
- **TRIPPED** `P4` — f_AB vs f_A x f_B over 19,110 measured pairs, rho=0.9916 (limit 0.98)
- clear `P5` — 0 pairs marked cleared on the optimistic arm
- clear `P6` — 0 recommended pairs are unmeasured
- clear `P7` — 0 cleared pairs contain a ubiquitous immune protein (in pool: none) []
- **TRIPPED** `P8` — 172 of 294 cleared pairs (58.5%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
- clear `P10` — 0 targets recommended dual despite clearing alone
- clear `P11` — outcome spread {'NO_DESIGN': 162, 'DUAL': 30, 'SINGLE': 3, 'ADAPTOR': 5}, largest 81.0% (limit 95%)
- clear `P12` — 0 of 30 dual recommendations change at 2 counts (0.0%, limit 50%)
- **TRIPPED** `P13` — most common partner takes 70.0% of dual recommendations (PRSS21)
- **TRIPPED** `P17` — 22 of 30 dual recommendations (73.3%) had exactly one admissible eligible partner, so no selection was made; limit 50%, above which the majority of recommendations name a partner the stage did not choose between
- clear `P14` — top pair is NPSR1+PTPRN2, top two singles are TMC5+ITGB6
- **TRIPPED** `P15` — pool halved to 100: 17 of 18 shared dual targets change partner (94.4%, limit 50%)

### Stage 4a — Architecture routing (11/12)

- clear `A1` — architecture is order-independent; 0 differ on a reversed pool
- clear `A2` — 0 targets routed ADAPTOR that CONVENTIONAL would have admitted
- clear `A3` — 0 routed risks differ from the Stage 3 risk
- clear `A4` — 0 CONVENTIONAL targets sit above the persistent ceiling
- clear `A5` — NPSR1 (risk 0.0277) routes CONVENTIONAL
- **TRIPPED** `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
- clear `A7` — 0 targets resolve to NO_ARCHITECTURE with no reason
- clear `A8` — with no declared terminable ceiling, 0 targets still route ADAPTOR
- clear `A9` — adaptor admissions across the ceiling sweep: 0.15->0, 0.2->0, 0.25->0, 0.3->5, 0.35->9, 0.4->26, 0.5->70, 0.6->113, 0.7->155
- clear `A10` — declared ceiling 0.35 matches the spec value 0.35
- clear `A11` — 5 adaptor constructs; 0 emit a sequence despite an unsupplied binder
- clear `A12` — 5 of 5 adaptor constructs carry a structure-derived binder the origin check can see; epitope immunogenicity stays NOT_CONNECTED on all of them, so the species gap and the immunogenicity gap remain separate

### Stage 5 — Binder discovery (7/7)

- clear `B1` — 77 targets have entries but no antibody among them (an inert stage would have none); 6 of 121 with entries echo their count
- clear `B3` — known answers hold: 5 sequence targets, 2 pinned structure targets, 2 documented negatives, and CEACAM5 at 1 structure entries of 6 (expected 1)
- clear `B8` — 5 of 5 known targets return a binder on some route
- clear `B5` — 0 candidates carry an affinity value (the source does not have one)
- clear `B7` — pool order carried through unchanged
- clear `B13` — 200 records out of 200 decisions in; 0 genes dropped, 0 added
- clear `B11` — 0 candidates claim a resolved isoform (neither route can determine one)

### Stage 6 — Construct assembly (8/9)

- clear `K0` — the decision set is routed (persistent 0.15, terminable 0.35) and yields 5 construct(s) for the criteria below to read
- clear `K1` — 5 constructs translate back to their own sequence
- **TRIPPED** `K2` — no dual carries a binder on both arms, so the two-arm join is not exercised anywhere in this decision set (5 of 200 rows assembled, 5 of them by the anti-tag route, whose binder is verified above); the criterion has nothing to pin on and reports that rather than clearing on an empty set
- clear `K3` — domain boundaries partition all 5 constructs exactly
- clear `K4` — every part of every construct names its source (10 parts in the first construct, 50 across all 5)
- clear `K5` — part costs sum to the printed total for all 5 constructs
- clear `K6` — all 5 buildable constructs carry the mandatory safety switch
- clear `K7` — every owed construct was built and none was built without the binder its architecture needs; 24 targets have a binder but no recommendation (§5.1)
- clear `K8` — 200 rows and 200 distinct genes against the 200 the Stage 4 manifest records

### Stage 9 — Safety gate (7/7)

- clear `S1` — registry returns trials for both pins (MSLN=169, CLDN18=3)
- clear `S2` — pinned names classify as expected (Amatuximab=chimeric, Zolbetuximab=chimeric)
- clear `S3` — every target with a binder carries a Stage 3 risk
- clear `S4` — no target escapes the ceiling applied to it: 3 admitted against the persistent 0.15, 5 against the terminable 0.35, each on a route declaring the exposure stoppable
- clear `S5` — epitope immunogenicity is NOT_CONNECTED on every row
- clear `S6` — 200 rows and 200 distinct genes against the 200 the Stage 4 manifest records
- clear `S7` — every target with stopped trials is flagged or blocked

### Stage 10 — Developability (6/6)

- clear `D1` — poly-K pI 11.463, poly-E pI 2.623 — the charge model orders the two known answers correctly
- clear `D2` — NST yields 1 sequon, NPT yields 0 — the proline exclusion holds
- clear `D3` — a 3-cysteine control reports parity 'odd', and never 'unpaired: 0'
- clear `D4` — every scored binder has a pI in 1..14 and 0..5 flags
- clear `D5` — 107 rows against 107 binders carrying a sequence
- clear `D6` — no liability is summed into a single score; flags are counted and listed

### Stage 11 — Final ranking (6/6)

- clear `N1` — a dominated point is excluded from the front
- clear `N2` — both non-dominated points are on the front
- clear `N3` — attrition accounts for 200 of 200
- clear `N4` — no weighted or summed score across objectives is emitted
- clear `N5` — status RANKED matches the survivor count 5
- clear `N6` — 200 rows against the 200 the Stage 4 manifest records

### Stage 12 — Candidate package (9/9)

- clear `Q1` — 5 package(s), one per surviving candidate, in the ranking's order: FER1L6, GPR35, TMEM92, TNFSF9, BTNL8
- clear `Q2` — every one of 5 candidates that reached the end is packaged, and none of the 195 that did not
- clear `Q3` — all 9 sections present on every package, and each carries what its stage produced
- clear `Q4` — the packaged DNA translates to the packaged sequence and the packaged domains partition it, for all 5
- clear `Q5` — every packaged attribution reconstructs its own risk to within 1e-12 and matches Stage 3
- clear `Q6` — 17 declared gap(s) probed and all still open; 1 recomputed from the run and 9 stated as judgements
- clear `Q7` — no conservative design exists in this pool and the section says so with the counts behind it, rather than standing blank
- clear `Q8` — 9 connected sources each name a release, and the hash chain is unbroken from Stage 3 to Stage 11
- clear `Q9` — no package emits a section or placeholder for Stage 7 or Stage 8; both are recorded in the gaps section instead (2 absent-stage entries)

### Stage API — HTTP surface (10/10)

- clear `A1` — project created (201), target_antigen None and discovery mode B
- clear `A2` — a view before any run answers 409 RUN_NOT_COMPLETE with instructions, not an empty list
- clear `A3` — a run returns 202 with job 69215abe99ab rather than blocking
- clear `A4` — job finished complete after stages ['sources', 'pairing', 'ranking']
- clear `A5` — 200 BUILDABLE: 5 buildable = 5 complete + 0 awaiting a binder; 0 over budget, 6 reasons
- clear `A6` — end state RANKED, attrition accounts for 195 + 5 of 200; 5 reached = 5 complete + 0 awaiting
- clear `A7` — top target CEACAM5 with a 6-component breakdown
- clear `A8` — pairs carry the span percentile beside the raw fraction (0.006321856890514115 at percentile 0.037)
- clear `A10` — an unknown project answers 404 NOT_FOUND and one that exists without a finished run answers 409 RUN_NOT_COMPLETE: a client can tell a bad id from a run in progress
- clear `A9` — evidence trail for MSLN spans 7 stages: stage3, stage4, stage5, stage6, stage9, stage10, stage11

### Stage MULTI — Multi-indication (11/11)

- clear `M1` — 52 indication-tagged artifacts; 0 changed, 0 disappeared after running both
- clear `M2` — shared sources carry no per-indication copy (0 found)
- clear `M3` — indication-specific module constants remaining: none
- clear `M4` — an atlas-less indication returns NOT_USABLE with no ranking
- clear `M5` — the refusal names malignant_vs_stroma as the missing discriminator, not just a lost weight
- clear `M6` — Mode A on CD19 returns NOT_ASSESSED with 3 reasons
- clear `M7` — CD19 ranks 1303 of 3400 -- outside the top 20, so the verdict is not self-agreement
- clear `M8` — Mode A and Mode B report the same evidence for CD19 (risk 0.5353, composite 0.227)
- clear `M9` — reference unchanged: top3 ['CEACAM5', 'TMC5', 'MUCL3'], pool 200, hash a91c696f2e1318f7, outcomes {'NO_DESIGN': 162, 'DUAL': 30, 'SINGLE': 3, 'ADAPTOR': 5}
- clear `M10` — a degraded indication names its missing source: ["dependency lineage 'NoSuchLineage': ValueError: need at least one array to concatenate"]
- clear `M11` — breast: 1548 rejected by the stromal gate, 0 of them on an absent measurement; known targets CEACAM6=138.7, MUC1=49.3, ERBB2=10.8, TACSTD2=12.5, MSLN=26.5, CEACAM5=exempt

## What the platform returns for this indication

```
    GET /constructs -> BUILDABLE
      - Every surviving design routes to an adaptor architecture. That is two manufactured biologics, not one: the receptor and, separately, the tagged adaptor antibody that gives it its specificity. The second carries its own CMC package and its own regulatory path, and the payload budget the adaptor route saves is paid there instead.
      - The anti-tag binder is a murine anti-GCN4 single-chain Fv retrieved from PDB 1P4B entities 1+2 at revision 1.4. PDB does not record the identifier 52SR4 anywhere in that entry. The identification rests on an exact match between the deposited CDRs and those quoted for 52SR4 in the Calibr/Scripps patent family, together with the shared Zahnd 2004 primary citation. That is an inference drawn by this pipeline, not a fact taken from the source, and a reader is entitled to disagree with it.
      - The retrieved binder is murine. The clinical construct in this tag system is humanized and its sequence is not established, so what is built here is the crystallised murine scFv, not the clinical one. Non-human sequence content is an explicit Stage 9 immunogenicity question and that arm is empty: epitope-level immunogenicity reports NOT_CONNECTED on every row because no epitope source is connected, and the origin check reads INN name stems, which a structure-derived binder does not carry. Nothing in this pipeline has assessed the immunogenicity of this binder.
      - The binder is emitted as deposited, including its crystallisation artifacts, because trimming them is a design decision this pipeline does not take silently. Each construct therefore carries MADYADA at residues 22-28, expression leader carried on the light-chain entity; and ASGADHHHHHH at residues 270-280, purification tag carried on the heavy-chain entity. As emitted these are not manufacturable: the first is a second leader sitting inside the mature protein, the second a His tag between the binder and the hinge. Removing them is a wet-lab step that has not been taken here.
      - No conservative backup exists in this pool. A conservative design is the conventional single-antigen receptor with a clinically-precedented binder, and no such design is buildable here: 3 single-antigen target(s) were recommended (MSLNL, NPSR1, ZPLD1) and none of them assembles, for want of a binder; no dual design assembles at all, because every dual recommendation names a partner that retrieves no binder. This is reported rather than filled by labelling something that does not qualify.
      - 5 advanced design(s) are available, all of them adaptor receptors, which is the architecture row the spec lists for serious normal-tissue expression.
    GET /result     -> RANKED
      blocked on normal tissue risk      - 192     8 remain
      no design recommended              -   0     8 remain
      no binder retrieved                -   3     5 remain
      no construct assembled             -   0     5 remain
      construct over budget              -   0     5 remain
```

## Tripped

- Stage 3 `R14` — no criticality tier places two staining levels on opposite sides of the 0.15 ceiling, so the arm gates on presence only — tier1 0.288/0.379/0.460(all above)  tier2 0.173/0.227/0.276(all above)  tier3 0.086/0.114/0.138(all below)
  - **Not on the accepted list.** Open decision, not a regression.
  - The staining arm vetoes on presence rather than amount: every grade blocks in tiers 1 and 2, none reaches the ceiling in tier 3. Under a conservative tolerance that may be the right design. It is a decision about tolerance, priced in reports/staining-veto-decision.md, and it is not on the accepted list because nobody has taken it yet. The non-zero exit is this decision being open, not a broken build.
- Stage 4 `P4` — f_AB vs f_A x f_B over 19,110 measured pairs, rho=0.9916 (limit 0.98)
  - Accepted: Coverage is span-confounded: f_AB tracks genomic span (+0.68) more than expression (+0.20). Reported beside a span-matched percentile and removed from partner selection.
- Stage 4 `P8` — 172 of 294 cleared pairs (58.5%) stop clearing if the unmeasured antigen saturates its organ (limit 10%)
  - Accepted: 48.6% of cleared pairs stop clearing if an unmeasured antigen saturates its organ. This is the cost of treating missing as a third state instead of imputing it.
- Stage 4 `P13` — most common partner takes 70.0% of dual recommendations (PRSS21)
  - Accepted: Accepted on the evidence in specs/p13-partner-concentration.md. P13 measures partner concentration as a share of dual recommendations, in a set where concentration is structural: 290 of 19,900 pairs are admissible, and 21 of the 21 targets that took the hub had exactly one admissible eligible partner. The statistic therefore cannot distinguish a rule preference from scarce supply, and no change to the selection objective moves it -- forcing every target that had a choice to choose otherwise leaves the share at 70.0%. Accepted rather than fixed, with P17 added to measure the property that is actually true.
- Stage 4 `P17` — 22 of 30 dual recommendations (73.3%) had exactly one admissible eligible partner, so no selection was made; limit 50%, above which the majority of recommendations name a partner the stage did not choose between
  - **Not on the accepted list. This is new.**
- Stage 4 `P15` — pool halved to 100: 17 of 18 shared dual targets change partner (94.4%, limit 50%)
  - Accepted: Partner choice is unstable under pool halving (71.4%). The pairing stage is complete-with-limitations by decision.
- Stage 4a `A6` — MSLN (risk 0.6366, lung) routes NO_ARCHITECTURE
  - Accepted: A positive pin written before the run. It expected MSLN to route to an adaptor because it matches that row's condition in words: serious normal-tissue expression. It does not, because its measured risk 0.6366 is nearly twice the declared terminable ceiling of 0.35. Admitting it needs a ceiling near 0.65, which also admits about 120 others - a clinical policy decision, not a code change. The ceiling stays where the spec pinned it and A9 reports the whole sweep so the trade is visible.
- Stage 6 `K2` — no dual carries a binder on both arms, so the two-arm join is not exercised anywhere in this decision set (5 of 200 rows assembled, 5 of them by the anti-tag route, whose binder is verified above); the criterion has nothing to pin on and reports that rather than clearing on an empty set
  - **Not on the accepted list. This is new.**
