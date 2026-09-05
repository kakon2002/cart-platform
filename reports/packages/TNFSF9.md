# TNFSF9 — candidate package

**CAR-PDAC-004** · `P41273` · decision **BACKUP** · PASSED_ALL_GATES · INNOVATIVE_DESIGN · dominated · position 4 of 5 · PACKAGED

This package carries what the pipeline produced. Eight of the reference document's twelve deliverables have something to carry; what the other four are missing is named in **What this package cannot tell you**, at the end, rather than left out.

---

## 1 — Ranking

| field | value |
| --- | --- |
| candidate_id | CAR-PDAC-004 |
| gate_status | PASSED_ALL_GATES |
| decision | BACKUP |
| on Pareto front | no |
| position | 4 of 5 |
| position basis | stage4 composite order, inherited; not a ranking Stage 11 computed |


| objective | value |
| --- | --- |
| attractiveness | 0.4114 |
| safety_margin | -0.1391 |
| binder_count | 0 |
| cleanliness | 0 |

> No weighted total across objectives is emitted. Candidates are compared on a Pareto front, so a design better on one objective and worse on another is not silently averaged into a rank.

> Position is a display index carried from Stage 4's composite order. It is not a ranking Stage 11 computed, and nothing decisional reads it: the decision column reads front membership, which is the comparison this stage actually performs.

## 2 — Scorecard

| component | weight | state | value | source |
| --- | --- | --- | --- | --- |
| tumour_coverage | 0.18 | MEASURED | 0.382 | stage3.patient_prevalence |
| malignant_specificity | 0.16 | MEASURED | 0.1531 | stage3.malignant_vs_stroma |
| normal_tissue_safety | 0.16 | MEASURED | 0.174 | stage9 residual margin below the applied ceiling 0.35 |
| binder_quality | 0.12 | UNKNOWN | — | this is an adaptor design: the receptor binds a tag, no anti-tag binder is retrieved, and there is no binder to score |
| manufacturability | 0.1 | MEASURED | 0.1806 | stage6 headroom against the 3500 bp payload budget |
| developability | 0.1 | UNKNOWN | — | Stage 10 counts sequence liabilities and refuses to sum them into a score, because one flag fires on every binder in the pool and a sum would hide it. That decision stands; this specification does not overturn it, so the component has no value to carry |
| structural_feasibility | 0.08 | UNKNOWN | — | Stage 7 does not exist. It is buildable and unbuilt |
| functional_prediction | 0.06 | UNKNOWN | — | Stage 8 does not exist and is not buildable from what is connected: the required training inputs are not available |
| pairing_robustness | 0.04 | NOT_APPLICABLE | — | an ADAPTOR design names no partner, so pairing robustness is a question that does not arise. This is not missing evidence |


|  |  |
| --- | --- |
| weight version | wm-scoring-1 |
| applicable weight | 0.96 |
| measured weight | 0.6 |
| scored fraction | 0.625 |
| floor | 0.5 |
| evidence confidence | 0.85 |
| prediction uncertainty | UNKNOWN |
| confidence adjustment | 0.85 |
| overall score | 0.197136 |


**UNKNOWN on this candidate:** binder_quality, developability, structural_feasibility, functional_prediction. Each is named above with the reason it is missing. None is imputed.

**NOT_APPLICABLE on this candidate:** pairing_robustness. This is a question that does not arise for this design, not a gap in the evidence.

> evidence confidence applied at exponent 1.0; the prediction-uncertainty penalty is NOT applied because no uncertainty is measured. It is not treated as zero uncertainty, which would flatter an unmeasured candidate.

> Scored on 0.6000 of 0.9600 applicable weight (0.6250). The remaining components are named, not imputed, and the score is normalised over what was measured.

## 3 — Construct

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

The amino-acid sequence and the nucleotide map are in `TNFSF9.json` beside this file. The DNA is a map under one fixed codon per residue, so the boundaries above are exact. It is not a codon-optimised ordering sequence.

## 4 — Target evidence

|  |  |
| --- | --- |
| composite | 0.4114 |
| measured weight | 1.0 |
| evidence class | RNA_SUPPORTED |
| confidence | 0.85 |
| normal-tissue risk | 0.2891 (lung) |
| risk basis | transcript only |
| risk is a lower bound | True |
| tumour-side verdict | TUMOUR_DOMINANT |

### Where the risk came from

Risk 0.2890959284073049 on lung, ahead of the next organ by 4.327185419095159e-06, across 18 organs that scored.

| organ | weighted | score | tier | arm | staining | transcript |
| --- | --- | --- | --- | --- | --- | --- |
| lung | 0.2891 | 0.2891 | 1 | BASELINE | NOT_MEASURED | Lung 6.4 TPM 0.2891 |
| brain | 0.2891 | 0.2891 | 1 | BASELINE | NOT_MEASURED | Brain_Cerebellar_Hemisphere 6.4 TPM 0.2891 |
| gi_tract | 0.2155 | 0.3592 | 2 | BASELINE | NOT_MEASURED | Esophagus_Mucosa 11.0 TPM 0.3592 |
| liver | 0.159 | 0.159 | 1 | BASELINE | NOT_MEASURED | Liver_Portal_Tract 2.0 TPM 0.1590 |
| marrow_and_blood | 0.1148 | 0.1913 | 2 | BASELINE | NOT_MEASURED | Spleen 2.8 TPM 0.1913 |
| reproductive | 0.1129 | 0.3763 | 3 | BASELINE | NOT_MEASURED | Vagina 12.5 TPM 0.3763 |
| vascular | 0.1107 | 0.1107 | 1 | BASELINE | NOT_MEASURED | Artery_Coronary 1.1 TPM 0.1107 |
| nerve | 0.085 | 0.1416 | 2 | BASELINE | NOT_MEASURED | Nerve_Tibial 1.7 TPM 0.1416 |
| kidney | 0.0776 | 0.0776 | 1 | BASELINE | NOT_MEASURED | Kidney_Cortex 0.7 TPM 0.0776 |
| adipose | 0.0539 | 0.1797 | 3 | BASELINE | NOT_MEASURED | Adipose_Subcutaneous 2.5 TPM 0.1797 |
| endocrine | 0.053 | 0.0884 | 2 | BASELINE | NOT_MEASURED | Adrenal_Gland 0.8 TPM 0.0884 |
| skin | 0.047 | 0.1565 | 3 | BASELINE | NOT_MEASURED | Skin_Sun_Exposed_Lower_leg 1.9 TPM 0.1565 |
| heart | 0.0447 | 0.0447 | 1 | BASELINE | NOT_MEASURED | Heart_Atrial_Appendage 0.4 TPM 0.0447 |
| bladder | 0.0398 | 0.0664 | 2 | BASELINE | NOT_MEASURED | Bladder 0.6 TPM 0.0664 |
| breast | 0.0317 | 0.1057 | 3 | BASELINE | NOT_MEASURED | Breast_Mammary_Tissue 1.1 TPM 0.1057 |
| muscle | 0.0317 | 0.0528 | 2 | BASELINE | NOT_MEASURED | Muscle_Skeletal 0.4 TPM 0.0528 |
| salivary | 0.0237 | 0.0791 | 3 | BASELINE | NOT_MEASURED | Minor_Salivary_Gland 0.7 TPM 0.0791 |
| pancreas | 0.0095 | 0.0159 | 2 | BASELINE | NOT_MEASURED | Pancreas_Mixed_Cell 0.1 TPM 0.0159 |

Every score above recomputes from the measurement beside it, and the largest weighted value is the risk. `NOT_MEASURED` is a third state: it is not a zero and not a clean result.

## 5 — Binders

Stage 5 verdict: **NO_BINDER** · 9 structural entries examined

None.

> The two routes are reported apart and never summed. A target with a named therapeutic but no deposited structure is not a target without a binder.

> This receptor binds a tag, not the antigen, so its binding domain is not a Stage 5 record: it is anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded), retrieved from 1P4B_1+2 and named in the construct section. A Stage 5 verdict of NO_BINDER for this target means no antigen-specific binder was retrieved, which is a different statement from the receptor having no binder.

## 6 — Safety

|  |  |
| --- | --- |
| verdict | FLAGGED |
| risk against the applied ceiling | 0.2891 against 0.35 |
| peak organ | lung |
| binder origin | non-human |
| source organism | Mus musculus |
| epitope immunogenicity | NOT_CONNECTED |
| trials naming this symbol | 20, 5 stopped |

- admitted against the terminable ceiling 0.35, not the persistent 0.15: activation requires a separately dosed adaptor, so the exposure is stoppable
- the receptor carries a structure-derived binder, Mus musculus as deposited in 1P4B_1+2, treated as non-human because no humanised sequence is established for it. This is read from the deposition, not from a name stem, which a structure-derived binder does not carry
- epitope-level immunogenicity remains NOT_CONNECTED for it: no epitope source is connected, so the species is known and the immunogenicity is not
- 5 trial(s) mentioning this symbol were terminated, withdrawn or suspended: NCT01861938, NCT04937166, NCT05219578

## 7 — Developability

0 binder sequence(s) scored.

None.

> Nothing was scored for this candidate. Stage 10 assesses Stage 5 sequence-route binders, and this design carries none: its binding domain is anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded), retrieved from 1P4B_1+2, which the stage does not read.

> So no developability figure in this platform describes the binder this construct actually carries. That is an absence, not a clean result.

## 8 — Experimental validation plan

### Before any bench work

- The binder is retrieved from 1P4B_1+2 as deposited. Its crystallisation artifacts are still present and must be removed before synthesis; the construct as emitted is not the molecule to order.
- The binder is Mus musculus and no humanised sequence is established for it. No immunogenicity assessment exists in this pipeline: the epitope-level arm is NOT_CONNECTED. Immunogenicity is an open question, not a low risk.
- This is an adaptor design, so it is two products. The tagged adaptor antibody is a second biologic with its own CMC and regulatory path, and nothing below can be run without it.

### In vitro

| step | purpose | material | readout | measures | acceptance |
| --- | --- | --- | --- | --- | --- |
| expression and surface presentation | confirm the receptor reaches the T-cell surface intact | primary human T cells, donor number TO_BE_SET | flow cytometry for the receptor ectodomain | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| antigen binding | confirm the receptor engages its intended tag | the tagged adaptor bearing the targeting antibody | binding titration | apparent affinity, TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| antigen-dependent killing | confirm killing requires both the adaptor and TNFSF9 | TNFSF9-positive and TNFSF9-negative lines; for the adaptor design, a no-adaptor arm is the negative control | cytotoxicity and cytokine release | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| on-target off-tumour check against the declared risk organ | the risk model puts this target's peak normal-tissue signal in lung at 0.2891. That is a prediction from expression data and has not been tested on cells here. | primary cells or organoids from lung | cytotoxicity against the normal-tissue model | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| immunogenicity, currently unassessed | the binder is non-human and this pipeline has no epitope source connected, so nothing upstream has assessed it | donor PBMC panel, donor number TO_BE_SET | T-cell proliferation or an equivalent assay | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| safety switch function | confirm the switch clears the product when triggered | the transduced product from step 1 | viability after dimeriser exposure | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |

### In vivo

| step | purpose | model | readout | measures | acceptance |
| --- | --- | --- | --- | --- | --- |
| tumour control | whether the design controls a TNFSF9-positive tumour | xenograft, line and strain TO_BE_SET · untransduced T cells; receptor without adaptor; receptor with adaptor | tumour burden over time | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| adaptor dose dependence | the adaptor is what makes exposure terminable, which is the basis on which this design was admitted against the terminable ceiling rather than the persistent one | as above · adaptor dose titration including withdrawal | tumour burden and receptor engagement over time | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |
| tolerability | no normal-tissue toxicity has been measured anywhere in this pipeline; the risk figure is derived from expression data | TO_BE_SET; a model expressing the human antigen is required for this to mean anything · as above | weight, histopathology of the declared risk organ | TO_BE_MEASURED | TO_BE_SET_BEFORE_THE_RUN |

> This is a plan, not a result. Every quantity a step would produce is marked TO_BE_MEASURED, and every acceptance threshold is marked TO_BE_SET_BEFORE_THE_RUN, because setting one after seeing the measurement is how a criterion stops being a test.

> No step below has been run and no outcome is predicted here.

## 9 — Provenance

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
| stage11 | `2591188b434f185d` |

---

## What this package cannot tell you

26 elements the reference document asks for are not produced, across 8 deliverables. 16 are checked mechanically by the verifier, 1 recomputed from this run, and 9 are judgements that say so.

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

