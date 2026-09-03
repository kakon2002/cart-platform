# BTNL8 — candidate package

`Q6UX41` · rank 5 of 5 · ADVANCED · dominated · PACKAGED

This package carries what the pipeline produced. Eight of the reference document's twelve deliverables have something to carry; what the other four are missing is named in **What this package cannot tell you**, at the end, rather than left out.

---

## 1 — Ranking

| objective | value |
| --- | --- |
| attractiveness | 0.4059 |
| safety_margin | -0.1747 |
| binder_count | 0 |
| cleanliness | 0 |

> No weighted total across objectives is emitted. Candidates are compared on a Pareto front, so a design better on one objective and worse on another is not silently averaged into a rank.

## 2 — Construct

|  |  |
| --- | --- |
| architecture | adaptor, anti-tag receptor, antigen on the adaptor |
| verdict | BUILDABLE |
| length | 2868 bp of a 3500 bp payload budget (632 bp spare) |
| residues | 955 |
| safety switch | present |
| binding domain | anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded) |

### Domain map

| domain | residues | aa | bp | provenance | source |
| --- | --- | --- | --- | --- | --- |
| CD8A leader | 21 | 0-21 | 0-63 | proteome | P01732 1-21 |
| anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded) | 259 | 21-280 | 63-840 | structure | 1P4B_1+2 |
| CD8A hinge | 45 | 280-325 | 840-975 | proteome | P01732 138-182 |
| CD8A transmembrane | 21 | 325-346 | 975-1038 | proteome | P01732 183-203 |
| 4-1BB cytoplasmic | 42 | 346-388 | 1038-1164 | proteome | Q07011 214-255 |
| CD3zeta cytoplasmic | 113 | 388-501 | 1164-1503 | proteome | P20963 52-164 |
| T2A skip peptide | 18 | 501-519 | 1503-1557 | synthetic | synthetic, named literal |
| FKBP12 | 107 | 519-626 | 1557-1878 | proteome | P62942 2-108 |
| SGGGS linker | 5 | 626-631 | 1878-1893 | synthetic | synthetic, named literal |
| caspase-9 without CARD | 324 | 631-955 | 1893-2865 | proteome | P55211 93-416 |

> the anti-tag binder is retrieved from 1P4B_1+2, deposited revision 1.4, as deposited and unedited

The amino-acid sequence and the nucleotide map are in `BTNL8.json` beside this file. The DNA is a map under one fixed codon per residue, so the boundaries above are exact. It is not a codon-optimised ordering sequence.

## 3 — Target evidence

|  |  |
| --- | --- |
| composite | 0.4059 |
| measured weight | 1.0 |
| evidence class | PROTEIN_CONFIRMED |
| confidence | 1.0 |
| normal-tissue risk | 0.3247 (gi_tract) |
| risk basis | staining and transcript |
| risk is a lower bound | False |
| tumour-side verdict | TUMOUR_DOMINANT |

### Where the risk came from

Risk 0.3246700695302471 on gi_tract, ahead of the next organ by 0.09751778656247609, across 20 organs that scored.

| organ | weighted | score | tier | arm | staining | transcript |
| --- | --- | --- | --- | --- | --- | --- |
| gi_tract | 0.3247 | 0.5411 | 2 | BASELINE | small intestine High 0.4601 | Colon_Transverse_Mixed_Cell 41.0 TPM 0.5411 |
| marrow_and_blood | 0.2272 | 0.3786 | 2 | STAINING | bone marrow Medium 0.3786 | Whole_Blood 11.5 TPM 0.3657 |
| lung | 0.1957 | 0.1957 | 1 | BASELINE | bronchus Not detected 0.0000 | Lung 2.9 TPM 0.1957 |
| reproductive | 0.138 | 0.4601 | 3 | STAINING | testis High 0.4601 | Testis 6.9 TPM 0.2994 |
| endocrine | 0.0805 | 0.1342 | 2 | BASELINE | adrenal gland Not detected 0.0000 | Thyroid 1.5 TPM 0.1342 |
| liver | 0.0603 | 0.0603 | 1 | BASELINE | liver Not detected 0.0000 | Liver_Portal_Tract 0.5 TPM 0.0603 |
| kidney | 0.0354 | 0.0354 | 1 | BASELINE | kidney Not detected 0.0000 | Kidney_Cortex 0.3 TPM 0.0354 |
| adipose | 0.0321 | 0.107 | 3 | BASELINE | NOT_MEASURED | Adipose_Visceral_Omentum 1.1 TPM 0.1070 |
| pancreas | 0.0267 | 0.0446 | 2 | BASELINE | pancreas Not detected 0.0000 | Pancreas_Mixed_Cell 0.4 TPM 0.0446 |
| heart | 0.0176 | 0.0176 | 1 | BASELINE | heart muscle Not detected 0.0000 | Heart_Atrial_Appendage 0.1 TPM 0.0176 |
| brain | 0.0124 | 0.0124 | 1 | BASELINE | caudate Not detected 0.0000 | Brain_Substantia_nigra 0.1 TPM 0.0124 |
| breast | 0.0114 | 0.0379 | 3 | BASELINE | breast Not detected 0.0000 | Breast_Mammary_Tissue 0.3 TPM 0.0379 |
| vascular | 0.0081 | 0.0081 | 1 | BASELINE | NOT_MEASURED | Artery_Coronary 0.1 TPM 0.0081 |
| muscle | 0.0077 | 0.0129 | 2 | BASELINE | skeletal muscle Not detected 0.0000 | Muscle_Skeletal 0.1 TPM 0.0129 |
| bladder | 0.0063 | 0.0104 | 2 | BASELINE | urinary bladder Not detected 0.0000 | Bladder 0.1 TPM 0.0104 |
| salivary | 0.0047 | 0.0158 | 3 | BASELINE | salivary gland Not detected 0.0000 | Minor_Salivary_Gland 0.1 TPM 0.0158 |
| nerve | 0.0033 | 0.0055 | 2 | BASELINE | NOT_MEASURED | Nerve_Tibial 0.0 TPM 0.0055 |
| skin | 0.0028 | 0.0093 | 3 | BASELINE | skin 1 Not detected 0.0000 | Skin_Not_Sun_Exposed_Suprapubic 0.1 TPM 0.0093 |
| connective | 0.0 | 0.0 | 3 | STAINING | soft tissue 2 Not detected 0.0000 | NOT_MEASURED |
| mucosa | 0.0 | 0.0 | 2 | STAINING | nasopharynx Not detected 0.0000 | NOT_MEASURED |

Every score above recomputes from the measurement beside it, and the largest weighted value is the risk. `NOT_MEASURED` is a third state: it is not a zero and not a clean result.

## 4 — Binders

Stage 5 verdict: **NO_BINDER** · 1 structural entries examined

None.

> The two routes are reported apart and never summed. A target with a named therapeutic but no deposited structure is not a target without a binder.

> This receptor binds a tag, not the antigen, so its binding domain is not a Stage 5 record: it is anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded), retrieved from 1P4B_1+2 and named in the construct section. A Stage 5 verdict of NO_BINDER for this target means no antigen-specific binder was retrieved, which is a different statement from the receptor having no binder.

## 5 — Safety

|  |  |
| --- | --- |
| verdict | FLAGGED |
| risk against the applied ceiling | 0.3247 against 0.35 |
| peak organ | gi_tract |
| binder origin | non-human |
| source organism | Mus musculus |
| epitope immunogenicity | NOT_CONNECTED |
| trials naming this symbol | 1, 0 stopped |

- admitted against the terminable ceiling 0.35, not the persistent 0.15: activation requires a separately dosed adaptor, so the exposure is stoppable
- the receptor carries a structure-derived binder, Mus musculus as deposited in 1P4B_1+2, treated as non-human because no humanised sequence is established for it. This is read from the deposition, not from a name stem, which a structure-derived binder does not carry
- epitope-level immunogenicity remains NOT_CONNECTED for it: no epitope source is connected, so the species is known and the immunogenicity is not

## 6 — Developability

0 binder sequence(s) scored.

None.

> Nothing was scored for this candidate. Stage 10 assesses Stage 5 sequence-route binders, and this design carries none: its binding domain is anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded), retrieved from 1P4B_1+2, which the stage does not read.

> So no developability figure in this platform describes the binder this construct actually carries. That is an absence, not a clean result.

## 7 — Experimental validation plan

### Before any bench work

- The binder is retrieved from 1P4B_1+2 as deposited. Its crystallisation artifacts are still present and must be removed before synthesis; the construct as emitted is not the molecule to order.
- The binder is Mus musculus and no humanised sequence is established for it. No immunogenicity assessment exists in this pipeline: the epitope-level arm is NOT_CONNECTED. Immunogenicity is an open question, not a low risk.
- This is an adaptor design, so it is two products. The tagged adaptor antibody is a second biologic with its own CMC and regulatory path, and nothing below can be run without it.

### In vitro

| step | purpose | material | readout | measures | acceptance |
| --- | --- | --- | --- | --- | --- |
| expression and surface presentation | confirm the receptor reaches the T-cell surface intact | primary human T cells, donor number TO_BE_SET | flow cytometry for the receptor ectodomain | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| antigen binding | confirm the receptor engages its intended tag | the tagged adaptor bearing the targeting antibody | binding titration | apparent affinity, TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| antigen-dependent killing | confirm killing requires both the adaptor and BTNL8 | BTNL8-positive and BTNL8-negative lines; for the adaptor design, a no-adaptor arm is the negative control | cytotoxicity and cytokine release | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| on-target off-tumour check against the declared risk organ | the risk model puts this target's peak normal-tissue signal in gi_tract at 0.3247. That is a prediction from expression data and has not been tested on cells here. | primary cells or organoids from gi_tract | cytotoxicity against the normal-tissue model | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| immunogenicity, currently unassessed | the binder is non-human and this pipeline has no epitope source connected, so nothing upstream has assessed it | donor PBMC panel, donor number TO_BE_SET | T-cell proliferation or an equivalent assay | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| safety switch function | confirm the switch clears the product when triggered | the transduced product from step 1 | viability after dimeriser exposure | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |

### In vivo

| step | purpose | model | readout | measures | acceptance |
| --- | --- | --- | --- | --- | --- |
| tumour control | whether the design controls a BTNL8-positive tumour | xenograft, line and strain TO_BE_SET · untransduced T cells; receptor without adaptor; receptor with adaptor | tumour burden over time | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| adaptor dose dependence | the adaptor is what makes exposure terminable, which is the basis on which this design was admitted against the terminable ceiling rather than the persistent one | as above · adaptor dose titration including withdrawal | tumour burden and receptor engagement over time | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| tolerability | no normal-tissue toxicity has been measured anywhere in this pipeline; the risk figure is derived from expression data | TO_BE_SET; a model expressing the human antigen is required for this to mean anything · as above | weight, histopathology of the declared risk organ | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |

> This is a plan, not a result. Every quantity a step would produce is marked TO_BE_MEASURED, and every acceptance threshold is marked TO_BE_SET_BEFORE_THE_RUN, because setting one after seeing the measurement is how a criterion stops being a test.

> No step below has been run and no outcome is predicted here.

## 8 — Provenance

| source | release | role |
| --- | --- | --- |
| UniProt | 2026_02 | reviewed human proteome |
| Human Protein Atlas | v23 | normal-tissue and pathology atlas |
| GTEx | v10 | bulk normal-tissue baseline |
| GENCODE | 47 | gene span annotation |
| SAbDab | 2.0.10 | antibody structures and named therapeutics |
| ClinicalTrials.gov | v2 | trial registry |
| GDC TCGA | TCGA-PAAD | tumour cohort |
| DepMap | 24Q4 / Pancreas | dependency lineage |
| Single-cell tumour atlas | GSE202051 | single-cell tumour atlas |

Configuration hash chain, each covering the stage before it:

| stage | hash |
| --- | --- |
| stage3 | `a91c696f2e1318f7` |
| stage4 | `5d097e05887e5b28` |
| stage5 | `6418657ed85a4dc6` |
| stage6 | `859cfd24c21ddd81` |
| stage9 | `7b703a7cc4a2d0d6` |
| stage10 | `518e1ef6953ee5fb` |
| stage11 | `7f2036bc6df78d1d` |

---

## What this package cannot tell you

27 elements the reference document asks for are not produced, across 8 deliverables. 17 are checked mechanically by the verifier, 1 recomputed from this run, and 9 are judgements that say so.

### Deliverable 1 — Top 3-5 CAR-T constructs (PARTIAL)

**the six named comparison views: maximum-efficacy, maximum-safety, best balanced, most manufacturable, lowest-cost, best universal-adaptor**

the Pareto front is computed and served, but its points are not labelled against these six. Two of them cannot be computed at all: there is no cost objective and no manufacturability score.

*Blocked by:* 8 and 10

**persistence, escape and uncertainty as ranking objectives**

four objectives are compared: tumour attractiveness, safety margin, binder count and binder cleanliness. Persistence and escape need Stage 8. Uncertainty is held as the separate confidence score and never combined into the ranking, by rule.

*Blocked by:* 8

### Deliverable 4 — Complete sequence and domain map (PARTIAL)

**predicted topology**

no topology is predicted for the assembled receptor.

*Blocked by:* 7

**expression risk**

nothing estimates whether the construct will express.

*Blocked by:* 7

**signalling-strength estimate**

the costimulatory and activation domains are assembled by accession and residue range; nothing estimates what they will signal.

*Blocked by:* 8

**recommended manufacturing format**

the construct reports its size against the payload budget and no format recommendation follows from it.

*Blocked by:* 10

### Deliverable 5 — Target and binder evidence report (PARTIAL)

**de novo binder generation, the whole of the document's section 5.2**

Stage 5 implements retrieval only: a structure route over deposited complexes and a sequence route over named therapeutics. Nothing generates a binder, models a complex, optimises an interface or germlines a framework.

*Blocked by:* none; this is unbuilt rather than blocked

**the binder counts: 20-100 initial, 10-20 computationally validated, 3-5 preferred per target**

these describe the de novo pipeline above. Retrieval returns what the literature holds for a target, which is not a quantity this platform chooses.

**predicted affinity**

the field exists on every candidate and is the constant NOT_CONNECTED: no affinity source is connected, and that is measured rather than assumed.

**epitope location**

what is recorded is the antigen chain and name of a deposited complex, which locates the antigen and not the epitope.

*Blocked by:* 7

**cross-reactivity risk**

no screen against paralogs, family members, normal-tissue proteins, alternative isoforms or polymorphic variants exists.

*Blocked by:* none; a paralog screen needs no stage that is missing

**human-likeness**

read from an INN name stem, which is a naming convention and not a sequence measurement, and which a structure-derived binder does not carry at all.

**uncertainty estimate per binder**

no candidate carries one.

**recommended status in the document's vocabulary**

the platform emits PROTEIN_CONFIRMED, RNA_SUPPORTED or DATA_INSUFFICIENT as an evidence class, and SINGLE, DUAL, ADAPTOR, NO_DESIGN or UNRESOLVED as an outcome. The document asks for high-confidence single, conditional, dual candidate, safety-gated or rejected. These map loosely and are not the same partition.

### Deliverable 6 — Safety-risk matrix (PARTIAL)

**a safety score**

the gate emits a verdict and a measured risk with its peak organ. There is no score, and combining the risk with the confidence to make one is the move this platform does not make.

**cytokine-release and neurotoxicity risk**

expected activation intensity, costimulatory contribution, cytokine profile and expansion kinetics are all unmeasured.

*Blocked by:* 8

**genomic and construct safety: recombination-prone regions, cryptic splice sites, unwanted open reading frames, sequence repeats**

none of it is computed.

*Blocked by:* none

*Note:* This is the cheapest closable gap on the whole list. All four are sequence analysis over the DNA map Stage 6 already emits, with its nucleotide sequence, domain boundaries and per-part provenance in hand. No new data source, no model, no external call, and no stage that does not exist.

**editing-related risks for allogeneic products**

no gene-editing package is assembled, so there is nothing to assess.

*Blocked by:* Stage 6's optional editing module, not built

**epitope-level immunogenicity**

NOT_CONNECTED on every row, because no epitope source is connected. The species of a binder is known; its immunogenicity is not.

**an uncertainty level on the safety verdict**

evidence confidence is measured and is deliberately kept apart from normal-tissue risk. The two are never combined, so the verdict carries no uncertainty term by rule rather than by omission.

### Deliverable 7 — Structural report (ABSENT)

**the entire stage: ten structural models and eight key scores, from extracellular-domain prediction and binder-antigen complexes through epitope accessibility, membrane distance, scFv stability, VH/VL orientation, hinge flexibility, domain interference, oligomerisation and aggregation risk, to synapse geometry, misfolding and tonic-signalling risk**

there is no module, no dataclass, no field and no stub. The project handoff described Stages 7 and 8 as schema only; there is no schema.

*Blocked by:* 7

*Note:* Buildable. It needs structure prediction over sequences Stage 6 already emits, and no data source that is missing. Two of its scores have sequence-level proxies in Stage 10's aggregation-prone regions, which is not the structural claim.

### Deliverable 8 — Functional predictions (ABSENT)

**the entire stage: activation threshold, cytotoxic potential, cytokine-release profile, proliferation, persistence, exhaustion, serial killing, activation-induced cell death, tonic signalling, antigen-density sensitivity, resistance to immunosuppressive conditions, performance under repeated antigen exposure**

there is no module, no dataclass, no field and no stub.

*Blocked by:* 8

*Note:* Not buildable from what is connected. The reference document names partner-generated experimental data among the required training inputs, and no such data is connected. This is the one gap on the list that a decision alone cannot close.

### Deliverable 9 — Manufacturability assessment (PARTIAL)

**expression efficiency, surface-expression probability, tonic signalling, transduction compatibility, product complexity, expected manufacturing yield, release-testing complexity, cost-of-goods, scalability**

nine of the document's thirteen evaluation items have no computation and no connected source. What exists is sequence developability over the binder, plus construct length against the payload budget.

*Blocked by:* 10

**the six named outputs: manufacturability score, vector recommendation, autologous versus allogeneic suitability, critical process risks, simplified backup architecture, recommended analytical assays**

none is produced.

*Blocked by:* 10

*Note:* The manufacturability score is not merely unbuilt. Stage 10 refuses to sum its liability flags into a single number by standing decision, because a flag that fires on every input carries no information and a sum would hide that.

**any developability figure describing the binder these designs carry**

Stage 10 assesses Stage 5 sequence-route binders. None of the 5 candidates in this package carries one, so none is scored: FER1L6, GPR35, TMEM92, TNFSF9, BTNL8. The stage runs and reports its rows, and none of those rows describes a design that ships.

*Blocked by:* none; Stage 10 would have to read the structure route

*Note:* Measured from this run rather than declared, so it disappears the moment one shipping design carries a sequence-route binder.

### Deliverable 12 — Full evidence and decision audit trail (PARTIAL)

**publications linked to each recommendation**

no publication is linked anywhere. Dataset releases and the configuration-hash chain are carried in this package's provenance block; the literature behind a recommendation is not.

**an evidence graph with an entity model**

the trail is per gene and per stage. The document asks for a versioned graph over cancers, antigens, isoforms, organs, epitopes, binders, architectures, trials and toxicity events; nothing here is an entity model.

