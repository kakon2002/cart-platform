"""A validation plan for one candidate construct: what to measure, never what it will show."""

from __future__ import annotations

from car_pipeline.data.domains import STRUCTURE

TO_BE_MEASURED = "TO_BE_MEASURED"
TO_BE_SET = "TO_BE_SET_BEFORE_THE_RUN"


def _binder(construct) -> dict:
    """The receptor's binding domain, and where its residues came from."""
    for segment in construct.segments:
        if segment.provenance in (STRUCTURE, "stage5"):
            return {
                "name": segment.name,
                "provenance": segment.provenance,
                "accession": segment.accession or None,
                "residues": f"{segment.aa_start + 1}-{segment.aa_end}",
            }
    return {"name": construct.binder_name or None, "provenance": None,
            "accession": None, "residues": None}


def _safety_module(construct) -> dict:
    """The switch this construct carries, if it carries one."""
    parts = [s for s in construct.segments
             if "caspase" in s.name.lower() or "FKBP" in s.name]
    if not parts:
        return {"present": False, "parts": [], "residues": None}
    return {
        "present": True,
        "parts": [s.name for s in parts],
        "residues": f"{parts[0].aa_start + 1}-{parts[-1].aa_end}",
    }


def _blockers(construct, safety) -> list[str]:
    """What must be resolved before any bench work, stated as blockers."""
    out = []
    binder = _binder(construct)
    if binder["provenance"] == STRUCTURE:
        out.append(
            f"The binder is retrieved from {binder['accession']} as deposited. "
            "Its crystallisation artifacts are still present and must be "
            "removed before synthesis; the construct as emitted is not the "
            "molecule to order."
        )
        organism = getattr(safety, "binder_source_organism", "") if safety else ""
        out.append(
            f"The binder is {organism or 'non-human'} and no humanised sequence "
            "is established for it. No immunogenicity assessment exists in this "
            "pipeline: the epitope-level arm is NOT_CONNECTED. Immunogenicity "
            "is an open question, not a low risk."
        )
    if construct.outcome == "ADAPTOR":
        out.append(
            "This is an adaptor design, so it is two products. The tagged "
            "adaptor antibody is a second biologic with its own CMC and "
            "regulatory path, and nothing below can be run without it."
        )
    return out


def _in_vitro(construct, target, binder, safety) -> list[dict]:
    """The minimal in-vitro sequence, each step naming its readout."""
    gene = construct.gene
    steps = [
        {
            "step": "expression and surface presentation",
            "purpose": f"confirm the receptor reaches the T-cell surface intact",
            "material": "primary human T cells, donor number TO_BE_SET",
            "readout": "flow cytometry for the receptor ectodomain",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
        {
            "step": "antigen binding",
            "purpose": f"confirm the receptor engages its intended tag",
            "material": (
                "the tagged adaptor bearing the targeting antibody"
                if construct.outcome == "ADAPTOR"
                else f"recombinant {gene} ectodomain"),
            "readout": "binding titration",
            "measures": "apparent affinity, " + TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
        {
            "step": "antigen-dependent killing",
            "purpose": (
                f"confirm killing requires both the adaptor and {gene}"
                if construct.outcome == "ADAPTOR"
                else f"confirm killing requires {gene}"),
            "material": (
                f"{gene}-positive and {gene}-negative lines; for the adaptor "
                "design, a no-adaptor arm is the negative control"),
            "readout": "cytotoxicity and cytokine release",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
        {
            "step": "on-target off-tumour check against the declared risk organ",
            "purpose": (
                f"the risk model puts this target's peak normal-tissue signal "
                f"in {target.get('risk_organ') or 'an organ it could not name'}"
                + (f" at {target['risk']}" if target.get("risk") is not None else "")
                + ". That is a prediction from expression data and has not been "
                "tested on cells here."),
            "material": (
                f"primary cells or organoids from "
                f"{target.get('risk_organ') or 'the peak organ'}"),
            "readout": "cytotoxicity against the normal-tissue model",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
    ]
    if safety and getattr(safety, "binder_origin", "") == "non-human":
        steps.append({
            "step": "immunogenicity, currently unassessed",
            "purpose": (
                "the binder is non-human and this pipeline has no epitope "
                "source connected, so nothing upstream has assessed it"),
            "material": "donor PBMC panel, donor number TO_BE_SET",
            "readout": "T-cell proliferation or an equivalent assay",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        })
    if _safety_module(construct)["present"]:
        steps.append({
            "step": "safety switch function",
            "purpose": "confirm the switch clears the product when triggered",
            "material": "the transduced product from step 1",
            "readout": "viability after dimeriser exposure",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        })
    return steps


def _in_vivo(construct, target) -> list[dict]:
    """The minimal in-vivo sequence, contingent on the in-vitro result."""
    gene = construct.gene
    return [
        {
            "step": "tumour control",
            "purpose": f"whether the design controls a {gene}-positive tumour",
            "model": "xenograft, line and strain TO_BE_SET",
            "arms": (
                "untransduced T cells; receptor without adaptor; receptor with "
                "adaptor" if construct.outcome == "ADAPTOR"
                else "untransduced T cells; receptor"),
            "readout": "tumour burden over time",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
        {
            "step": "adaptor dose dependence" if construct.outcome == "ADAPTOR"
                    else "dose dependence",
            "purpose": (
                "the adaptor is what makes exposure terminable, which is the "
                "basis on which this design was admitted against the "
                "terminable ceiling rather than the persistent one"
                if construct.outcome == "ADAPTOR"
                else "establish the dose-response relationship"),
            "model": "as above",
            "arms": "adaptor dose titration including withdrawal"
                    if construct.outcome == "ADAPTOR" else "dose titration",
            "readout": "tumour burden and receptor engagement over time",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
        {
            "step": "tolerability",
            "purpose": (
                "no normal-tissue toxicity has been measured anywhere in this "
                "pipeline; the risk figure is derived from expression data"),
            "model": "TO_BE_SET; a model expressing the human antigen is "
                     "required for this to mean anything",
            "arms": "as above",
            "readout": "weight, histopathology of the declared risk organ",
            "measures": TO_BE_MEASURED,
            "acceptance": TO_BE_SET,
        },
    ]


def plan(construct, target: dict, safety=None) -> dict:
    """The validation plan for one construct. A template, never a prediction."""
    binder = _binder(construct)
    return {
        "status": "PLAN",
        "gene": construct.gene,
        "accession": construct.accession,
        "architecture": construct.architecture,
        "binder": binder,
        "safety_module": _safety_module(construct),
        "construct_bp": construct.total_bp,
        "target": {
            "gene": construct.gene,
            "risk": target.get("risk"),
            "risk_organ": target.get("risk_organ"),
            "evidence_class": target.get("evidence_class"),
        },
        "before_any_bench_work": _blockers(construct, safety),
        "in_vitro": _in_vitro(construct, target, binder, safety),
        "in_vivo": _in_vivo(construct, target),
        "reasons": [
            "This is a plan, not a result. Every quantity a step would produce "
            "is marked TO_BE_MEASURED, and every acceptance threshold is marked "
            "TO_BE_SET_BEFORE_THE_RUN, because setting one after seeing the "
            "measurement is how a criterion stops being a test.",
            "No step below has been run and no outcome is predicted here.",
        ],
    }


CONSERVATIVE = "CONSERVATIVE_BACKUP"
ADVANCED = "ADVANCED"


def design_class(construct) -> str | None:
    """Which of the two the client asked for this design is, if either."""
    if construct.verdict != "BUILDABLE":
        return None
    clinical = any(s.provenance == "stage5" for s in construct.segments)
    if construct.outcome == "SINGLE" and clinical:
        return CONSERVATIVE
    if construct.outcome in ("DUAL", "ADAPTOR"):
        return ADVANCED
    return None


def design_class_summary(constructs) -> dict:
    """What the pool can and cannot supply, named rather than approximated."""
    labelled = [(c, design_class(c)) for c in constructs]
    conservative = [c.gene for c, k in labelled if k == CONSERVATIVE]
    advanced = [c.gene for c, k in labelled if k == ADVANCED]
    reasons = []
    if not conservative:
        singles = [c for c in constructs if c.outcome == "SINGLE"]
        singles_built = [c for c in singles if c.amino_acid_sequence]
        duals_built = [c for c in constructs
                       if c.outcome == "DUAL" and c.amino_acid_sequence]
        over = [c for c in duals_built if c.verdict != "BUILDABLE"]
        names = ", ".join(sorted(c.gene for c in singles))
        reasons.append(
            "No conservative backup exists in this pool. A conservative design "
            "is the conventional single-antigen receptor with a "
            "clinically-precedented binder, and no such design is buildable "
            "here: "
            + ("no target clears the ceiling alone" if not singles else
               f"{len(singles)} single-antigen target(s) were recommended "
               f"({names}) and none of them assembles, for want of a binder"
               if not singles_built else
               f"{len(singles)} single-antigen target(s) were recommended "
               f"({names}) and the {len(singles_built)} that assemble carry no "
               "binder retrieved as a named therapeutic")
            + "; "
            + (f"{len(duals_built)} dual design(s) assemble, {len(over)} of "
               "them over the payload budget" if duals_built else
               "no dual design assembles at all, because every dual "
               "recommendation names a partner that retrieves no binder")
            + ". This is reported rather than filled by labelling something "
            "that does not qualify."
        )
    if advanced:
        reasons.append(
            f"{len(advanced)} advanced design(s) are available, all of them "
            "adaptor receptors, which is the architecture row the spec lists "
            "for serious normal-tissue expression."
        )
    return {
        "conservative_backup": conservative or None,
        "advanced": advanced or None,
        "reasons": reasons,
    }
