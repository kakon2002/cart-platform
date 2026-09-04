"""Stage 12 — the final candidate package, assembled and never recomputed."""

from __future__ import annotations

from dataclasses import dataclass

from car_pipeline.data import domains
from car_pipeline.stages import stage6, validation

PACKAGED = "PACKAGED"
NO_CANDIDATE_REACHES_THE_END = "NO_CANDIDATE_REACHES_THE_END"

ABSENT = "ABSENT"
PARTIAL = "PARTIAL"

PROBE_MODULE = "module"
PROBE_FIELD = "field"
PROBE_KEY = "key"

JUDGEMENT = "no probe: this claim is a judgement, not a presence"


@dataclass(frozen=True)
class Gap:
    """One element the reference document asks for and the platform does not produce."""

    deliverable: int
    title: str
    state: str
    element: str
    reason: str
    blocking_stage: str | None = None
    probe: tuple[str, ...] | None = None
    note: str = ""

    def as_payload(self) -> dict:
        """The gap, with its probe stated so a reader can check it too."""
        return {
            "deliverable": self.deliverable,
            "title": self.title,
            "state": self.state,
            "element": self.element,
            "reason": self.reason,
            "blocking_stage": self.blocking_stage,
            "probe": list(self.probe) if self.probe else None,
            "probe_note": None if self.probe else JUDGEMENT,
            "note": self.note or None,
        }


GAPS: tuple[Gap, ...] = (
    Gap(1, "Top 3-5 CAR-T constructs", PARTIAL,
        "the six named comparison views: maximum-efficacy, maximum-safety, "
        "best balanced, most manufacturable, lowest-cost, best universal-adaptor",
        "the Pareto front is computed and served, but its points are not "
        "labelled against these six. Two of them cannot be computed at all: "
        "there is no cost objective and no manufacturability score.",
        blocking_stage="8 and 10", probe=(PROBE_KEY, "ranking", "comparison_views")),
    Gap(1, "Top 3-5 CAR-T constructs", PARTIAL,
        "persistence, escape and uncertainty as ranking objectives",
        "four objectives are compared: tumour attractiveness, safety margin, "
        "binder count and binder cleanliness. Persistence and escape need "
        "Stage 8. Uncertainty is held as the separate confidence score and "
        "never combined into the ranking, by rule.",
        blocking_stage="8"),

    Gap(4, "Complete sequence and domain map", PARTIAL,
        "predicted topology",
        "no topology is predicted for the assembled receptor.",
        blocking_stage="7", probe=(PROBE_KEY, "construct", "predicted_topology")),
    Gap(4, "Complete sequence and domain map", PARTIAL,
        "expression risk",
        "nothing estimates whether the construct will express.",
        blocking_stage="7", probe=(PROBE_KEY, "construct", "expression_risk")),
    Gap(4, "Complete sequence and domain map", PARTIAL,
        "signalling-strength estimate",
        "the costimulatory and activation domains are assembled by accession "
        "and residue range; nothing estimates what they will signal.",
        blocking_stage="8", probe=(PROBE_KEY, "construct", "signalling_strength")),
    Gap(4, "Complete sequence and domain map", PARTIAL,
        "recommended manufacturing format",
        "the construct reports its size against the payload budget and no "
        "format recommendation follows from it.",
        blocking_stage="10", probe=(PROBE_KEY, "construct", "manufacturing_format")),

    Gap(5, "Target and binder evidence report", PARTIAL,
        "de novo binder generation, the whole of the document's section 5.2",
        "Stage 5 implements retrieval only: a structure route over deposited "
        "complexes and a sequence route over named therapeutics. Nothing "
        "generates a binder, models a complex, optimises an interface or "
        "germlines a framework.",
        blocking_stage="none; this is unbuilt rather than blocked",
        probe=(PROBE_KEY, "binders", "generated")),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "the binder counts: 20-100 initial, 10-20 computationally validated, "
        "3-5 preferred per target",
        "these describe the de novo pipeline above. Retrieval returns what the "
        "literature holds for a target, which is not a quantity this platform "
        "chooses."),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "predicted affinity",
        "the field exists on every candidate and is the constant "
        "NOT_CONNECTED: no affinity source is connected, and that is measured "
        "rather than assumed."),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "epitope location",
        "what is recorded is the antigen chain and name of a deposited "
        "complex, which locates the antigen and not the epitope.",
        blocking_stage="7", probe=(PROBE_KEY, "binders", "epitope")),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "cross-reactivity risk",
        "no screen against paralogs, family members, normal-tissue proteins, "
        "alternative isoforms or polymorphic variants exists.",
        blocking_stage="none; a paralog screen needs no stage that is missing",
        probe=(PROBE_KEY, "binders", "cross_reactivity")),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "human-likeness",
        "read from an INN name stem, which is a naming convention and not a "
        "sequence measurement, and which a structure-derived binder does not "
        "carry at all."),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "uncertainty estimate per binder",
        "no candidate carries one.",
        probe=(PROBE_KEY, "binders", "uncertainty")),
    Gap(5, "Target and binder evidence report", PARTIAL,
        "recommended status in the document's vocabulary",
        "the platform emits PROTEIN_CONFIRMED, RNA_SUPPORTED or "
        "DATA_INSUFFICIENT as an evidence class, and SINGLE, DUAL, ADAPTOR, "
        "NO_DESIGN or UNRESOLVED as an outcome. The document asks for "
        "high-confidence single, conditional, dual candidate, safety-gated or "
        "rejected. These map loosely and are not the same partition."),

    Gap(6, "Safety-risk matrix", PARTIAL,
        "a safety score",
        "the gate emits a verdict and a measured risk with its peak organ. "
        "There is no score, and combining the risk with the confidence to make "
        "one is the move this platform does not make.",
        probe=(PROBE_FIELD, "car_pipeline.stages.stage9.SafetyRecord",
               "safety_score")),
    Gap(6, "Safety-risk matrix", PARTIAL,
        "cytokine-release and neurotoxicity risk",
        "expected activation intensity, costimulatory contribution, cytokine "
        "profile and expansion kinetics are all unmeasured.",
        blocking_stage="8",
        probe=(PROBE_FIELD, "car_pipeline.stages.stage9.SafetyRecord",
               "cytokine_risk")),
    Gap(6, "Safety-risk matrix", PARTIAL,
        "editing-related risks for allogeneic products",
        "no gene-editing package is assembled, so there is nothing to assess.",
        blocking_stage="Stage 6's optional editing module, not built"),
    Gap(6, "Safety-risk matrix", PARTIAL,
        "epitope-level immunogenicity",
        "NOT_CONNECTED on every row, because no epitope source is connected. "
        "The species of a binder is known; its immunogenicity is not."),
    Gap(6, "Safety-risk matrix", PARTIAL,
        "an uncertainty level on the safety verdict",
        "evidence confidence is measured and is deliberately kept apart from "
        "normal-tissue risk. The two are never combined, so the verdict "
        "carries no uncertainty term by rule rather than by omission."),

    Gap(9, "Manufacturability assessment", PARTIAL,
        "expression efficiency, surface-expression probability, tonic "
        "signalling, transduction compatibility, product complexity, expected "
        "manufacturing yield, release-testing complexity, cost-of-goods, "
        "scalability",
        "nine of the document's thirteen evaluation items have no computation "
        "and no connected source. What exists is sequence developability over "
        "the binder, plus construct length against the payload budget.",
        blocking_stage="10", probe=(PROBE_KEY, "developability", "process")),
    Gap(9, "Manufacturability assessment", PARTIAL,
        "the six named outputs: manufacturability score, vector "
        "recommendation, autologous versus allogeneic suitability, critical "
        "process risks, simplified backup architecture, recommended "
        "analytical assays",
        "none is produced.",
        blocking_stage="10", probe=(PROBE_KEY, "developability", "assessment"),
        note="The manufacturability score is not merely unbuilt. Stage 10 "
             "refuses to sum its liability flags into a single number by "
             "standing decision, because a flag that fires on every input "
             "carries no information and a sum would hide that."),

    Gap(7, "Structural report", ABSENT,
        "the entire stage: ten structural models and eight key scores, from "
        "extracellular-domain prediction and binder-antigen complexes through "
        "epitope accessibility, membrane distance, scFv stability, VH/VL "
        "orientation, hinge flexibility, domain interference, oligomerisation "
        "and aggregation risk, to synapse geometry, misfolding and "
        "tonic-signalling risk",
        "there is no module, no dataclass, no field and no stub. The project "
        "handoff described Stages 7 and 8 as schema only; there is no schema.",
        blocking_stage="7", probe=(PROBE_MODULE, "stage7"),
        note="Buildable. It needs structure prediction over sequences Stage 6 "
             "already emits, and no data source that is missing. Two of its "
             "scores have sequence-level proxies in Stage 10's "
             "aggregation-prone regions, which is not the structural claim."),

    Gap(8, "Functional predictions", ABSENT,
        "the entire stage: activation threshold, cytotoxic potential, "
        "cytokine-release profile, proliferation, persistence, exhaustion, "
        "serial killing, activation-induced cell death, tonic signalling, "
        "antigen-density sensitivity, resistance to immunosuppressive "
        "conditions, performance under repeated antigen exposure",
        "there is no module, no dataclass, no field and no stub.",
        blocking_stage="8", probe=(PROBE_MODULE, "stage8"),
        note="Not buildable from what is connected. The reference document "
             "names partner-generated experimental data among the required "
             "training inputs, and no such data is connected. This is the one "
             "gap on the list that a decision alone cannot close."),

    Gap(12, "Full evidence and decision audit trail", PARTIAL,
        "publications linked to each recommendation",
        "no publication is linked anywhere. Dataset releases and the "
        "configuration-hash chain are carried in this package's provenance "
        "block; the literature behind a recommendation is not.",
        probe=(PROBE_KEY, "provenance", "publications")),
    Gap(12, "Full evidence and decision audit trail", PARTIAL,
        "an evidence graph with an entity model",
        "the trail is per gene and per stage. The document asks for a "
        "versioned graph over cancers, antigens, isoforms, organs, epitopes, "
        "binders, architectures, trials and toxicity events; nothing here is "
        "an entity model."),
)


def measured_gaps(packages: list[dict]) -> list[Gap]:
    """Gaps that are properties of this run, recomputed rather than declared."""
    if not packages:
        return []
    unscored = [p["gene"] for p in packages
                if p["developability"]["binders_scored"] == 0]
    if len(unscored) != len(packages):
        return []
    return [Gap(
        9, "Manufacturability assessment", PARTIAL,
        "any developability figure describing the binder these designs carry",
        f"Stage 10 assesses Stage 5 sequence-route binders. None of the "
        f"{len(packages)} candidates in this package carries one, so none is "
        f"scored: {', '.join(unscored)}. The stage runs and reports its rows, "
        f"and none of those rows describes a design that ships.",
        blocking_stage="none; Stage 10 would have to read the structure route",
        note="Measured from this run rather than declared, so it disappears "
             "the moment one shipping design carries a sequence-route binder.",
    )]


def gap_payload(packages: list[dict] | None = None) -> dict:
    """The gaps section: what is missing, per deliverable, with its probe."""
    gaps = GAPS + tuple(measured_gaps(packages or []))
    by_deliverable: dict[int, dict] = {}
    for gap in gaps:
        entry = by_deliverable.setdefault(gap.deliverable, {
            "deliverable": gap.deliverable,
            "title": gap.title,
            "state": gap.state,
            "missing": [],
        })
        entry["missing"].append(gap.as_payload())
    ordered = [by_deliverable[k] for k in sorted(by_deliverable)]
    probed = sum(1 for g in gaps if g.probe)
    measured = sum(1 for g in gaps if g not in GAPS)
    return {
        "deliverables_with_gaps": len(ordered),
        "elements_missing": len(gaps),
        "elements_probed": probed,
        "elements_measured": measured,
        "elements_asserted": len(gaps) - probed - measured,
        "by_deliverable": ordered,
        "reasons": [
            "Every element here is something the reference document asks for "
            "and this platform does not produce. It is named rather than "
            "omitted, because a package that carried eight sections silently "
            "would read as twelve.",
            f"{probed} of {len(gaps)} entries carry a probe the verifier "
            f"executes and {measured} are recomputed from this run; the "
            f"remaining {len(gaps) - probed - measured} are judgements and say "
            "so. A gap that stops being true trips criterion Q6.",
        ],
    }


def _segment_payload(segment) -> dict:
    """One domain, with the boundary and the provenance it was built from."""
    return {
        "name": segment.name,
        "provenance": segment.provenance,
        "accession": segment.accession or None,
        "feature": segment.feature or None,
        "source_residues": (
            f"{segment.start_residue}-{segment.end_residue}"
            if segment.start_residue else None),
        "aa_start": segment.aa_start,
        "aa_end": segment.aa_end,
        "bp_start": segment.bp_start,
        "bp_end": segment.bp_end,
        "residues": segment.residues,
    }


def _construct_payload(construct) -> dict:
    """Deliverable 4, carried whole from the construct stage."""
    return {
        "verdict": construct.verdict,
        "architecture": construct.architecture,
        "outcome": construct.outcome,
        "partner": construct.partner,
        "binder_name": construct.binder_name or None,
        "partner_binder_name": construct.partner_binder_name or None,
        "amino_acid_sequence": construct.amino_acid_sequence,
        "residues": len(construct.amino_acid_sequence),
        "dna": construct.dna,
        "total_bp": construct.total_bp,
        "budget_bp": stage6.BUDGET_BP,
        "headroom_bp": construct.headroom_bp,
        "safety_switch": construct.has_switch,
        "domains": [_segment_payload(s) for s in construct.segments],
        "reason": construct.reason or None,
    }


def _target_payload(ranked, attribution) -> dict:
    """Deliverables 5 and 11 on the target side, including the per-organ map."""
    return {
        "composite": ranked.composite,
        "composite_supported": ranked.composite_supported,
        "measured_weight": ranked.measured_weight,
        "evidence_class": ranked.evidence_class,
        "tier_rank": ranked.tier_rank,
        "confidence": ranked.confidence,
        "below_evidence_floor": ranked.below_evidence_floor,
        "sources_disagree": ranked.sources_disagree,
        "tumour_side_verdict": ranked.tumour_side_verdict,
        "protein_arm_measured": ranked.protein_arm_measured,
        "risk": ranked.risk,
        "risk_organ": ranked.risk_organ,
        "risk_basis": ranked.risk_basis,
        "risk_is_lower_bound": ranked.risk_is_lower_bound,
        "cleared": ranked.cleared,
        "components": {k: ranked.component_value(k) for k in ranked.components},
        "risk_attribution": attribution.as_payload(),
    }


def _candidate_payload(candidate) -> dict:
    """One retrieved binder, on the route that retrieved it."""
    return {
        "route": candidate.route,
        "identifier": candidate.identifier,
        "name": candidate.name or None,
        "format": candidate.fmt or None,
        "clinical_stage": candidate.clinical_stage or None,
        "status": candidate.status or None,
        "heavy_residues": len(candidate.heavy_sequence) or None,
        "light_residues": len(candidate.light_sequence) or None,
        "antigen_chain": candidate.antigen_chain or None,
        "antigen_name": candidate.antigen_name or None,
        "method": candidate.method or None,
        "affinity": candidate.affinity,
        "isoform": candidate.isoform,
    }


def _binding_domain(construct):
    """The structure-derived segment carrying the binder, if there is one."""
    return next((s for s in construct.segments
                 if s.provenance == domains.STRUCTURE), None)


def _adaptor_binder_note(construct) -> str | None:
    """Where an adaptor receptor's binder came from, since Stage 5 has none."""
    carried = _binding_domain(construct)
    if carried is None:
        return None
    return (
        "This receptor binds a tag, not the antigen, so its binding domain is "
        f"not a Stage 5 record: it is {carried.name}, retrieved from "
        f"{carried.accession} and named in the construct section. A Stage 5 "
        "verdict of NO_BINDER for this target means no antigen-specific "
        "binder was retrieved, which is a different statement from the "
        "receptor having no binder."
    )


def _binder_payload(record, construct) -> dict:
    """Deliverable 5 on the binder side, by route and never summed."""
    note = _adaptor_binder_note(construct)
    if record is None:
        reasons = ["No Stage 5 record exists for this gene."]
        if note:
            reasons.append(note)
        return {
            "verdict": None,
            "entries": 0,
            "structure_route": [],
            "sequence_route": [],
            "reasons": reasons,
        }
    reasons = [
        "The two routes are reported apart and never summed. A target with a "
        "named therapeutic but no deposited structure is not a target without "
        "a binder.",
    ]
    if note:
        reasons.append(note)
    return {
        "verdict": record.verdict,
        "entries": len(record.entries),
        "entries_without_antibody": record.entries_without_antibody,
        "entries_excluded_as_model": record.entries_excluded_as_model,
        "structure_route": [_candidate_payload(c) for c in record.structure],
        "sequence_route": [_candidate_payload(c) for c in record.sequence],
        "reasons": reasons,
    }


def _safety_payload(record) -> dict:
    """Deliverable 6, carried whole from the gate."""
    if record is None:
        return {"verdict": None, "reasons": ["No Stage 9 record for this gene."]}
    return {
        "verdict": record.verdict,
        "risk": record.risk,
        "risk_organ": record.risk_organ,
        "ceiling_applied": record.ceiling,
        "binder_name": record.binder_name or None,
        "binder_origin": record.binder_origin,
        "binder_origins": record.binder_origins,
        "binder_source_organism": record.binder_source_organism or None,
        "binder_structure_accession": record.binder_structure_accession or None,
        "epitope_immunogenicity": record.epitope_immunogenicity,
        "trials_total": record.trials_total,
        "trials_stopped": record.trials_stopped,
        "trials_stopped_ids": record.trials_stopped_ids,
        "trials_truncated": record.trials_truncated,
        "construct_safety": record.construct_safety,
        "reasons": record.reasons,
    }


def _developability_payload(rows, construct) -> dict:
    """Deliverable 9, flags listed and never summed into a score."""
    if not rows:
        carried = _binding_domain(construct)
        why = (
            f"its binding domain is {carried.name}, retrieved from "
            f"{carried.accession}, which the stage does not read"
            if carried is not None else
            "this construct carries no binder residues at all: its binding "
            "domain declares a size and no sequence")
        return {
            "binders_scored": 0,
            "rows": [],
            "reasons": [
                "Nothing was scored for this candidate. Stage 10 assesses "
                f"Stage 5 sequence-route binders, and this design carries "
                f"none: {why}.",
                "So no developability figure in this platform describes the "
                "binder this construct actually carries. That is an absence, "
                "not a clean result.",
            ],
        }
    return {
        "binders_scored": len(rows),
        "rows": [
            {
                "binder": r.binder,
                "residues": r.residues,
                "isoelectric_point": r.isoelectric_point,
                "net_charge": r.net_charge,
                "cysteines": r.cysteines,
                "cysteine_parity": r.cysteine_parity,
                "glycosylation_sequons": r.sequons,
                "aggregation_prone_starts": r.apr_starts,
                "gravy": r.gravy,
                "flags": [{"kind": k, "detail": d} for k, d in r.flags],
                "flag_count": r.flag_count,
            }
            for r in rows
        ],
        "reasons": [
            "Flags are counted and listed, never summed into a developability "
            "score. One of them fires on every binder in this pool, which a "
            "sum would hide.",
        ],
    }


def _ranking_payload(entry, position, total) -> dict:
    """Deliverable 1, carried from the ranking with no re-ordering."""
    return {
        "position": position,
        "of": total,
        "on_pareto_front": entry.on_front,
        "objectives": {
            "attractiveness": entry.attractiveness,
            "safety_margin": entry.safety_margin,
            "binder_count": entry.binder_count,
            "cleanliness": entry.cleanliness,
        },
        "binder_supplied": entry.binder_supplied,
        "reasons": [
            "No weighted total across objectives is emitted. Candidates are "
            "compared on a Pareto front, so a design better on one objective "
            "and worse on another is not silently averaged into a rank.",
        ],
    }


def build(run: dict) -> tuple[list[dict], str]:
    """One package per candidate that reached the end, in the ranking's order."""
    final = run.get("final") or []
    survivors = [r for r in final if r.survived]
    if not survivors:
        return [], NO_CANDIDATE_REACHES_THE_END

    by_construct = {c.gene: c for c in run["constructs"]}
    by_gate = {g.gene: g for g in run["gated"]}
    by_ranked = {r.gene: r for r in run["ranked"] if r.gene}
    dev_by_gene: dict[str, list] = {}
    for row in run["developability"]:
        dev_by_gene.setdefault(row.gene, []).append(row)
    inputs = run["risk_inputs"]

    packages = []
    for position, entry in enumerate(survivors, 1):
        gene = entry.gene
        construct = by_construct.get(gene)
        ranked = by_ranked.get(gene)
        safety = by_gate.get(gene)
        attribution = inputs.attribute(entry.accession, gene)
        target = {
            "risk": getattr(ranked, "risk", None),
            "risk_organ": getattr(ranked, "risk_organ", None),
            "evidence_class": getattr(ranked, "evidence_class", None),
        }
        packages.append({
            "gene": gene,
            "accession": entry.accession,
            "ranking": _ranking_payload(entry, position, len(survivors)),
            "design_class": validation.design_class(construct),
            "construct": _construct_payload(construct),
            "target_evidence": _target_payload(ranked, attribution),
            "binders": _binder_payload(run["binders"].get(gene), construct),
            "safety": _safety_payload(safety),
            "developability": _developability_payload(
                dev_by_gene.get(gene, []), construct),
            "validation_plan": validation.plan(construct, target, safety),
            "provenance": run["provenance"],
        })
    return packages, PACKAGED


SECTIONS = (
    "ranking", "design_class", "construct", "target_evidence", "binders",
    "safety", "developability", "validation_plan", "provenance",
)
