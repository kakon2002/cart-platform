"""Writes one candidate package per surviving design, from the run's own output."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.stages import stage3, stage12

OUT = Path(__file__).resolve().parent / "reports" / "packages"

WRAP = 96


def _p(lines: list[str], text: str = "") -> None:
    """One paragraph."""
    lines.append(text)


def _table(lines: list[str], header: list[str], rows: list[list[str]]) -> None:
    """A markdown table, or a line saying the table is empty."""
    if not rows:
        _p(lines, "None.")
        _p(lines)
        return
    _p(lines, "| " + " | ".join(header) + " |")
    _p(lines, "| " + " | ".join("---" for _ in header) + " |")
    for row in rows:
        _p(lines, "| " + " | ".join(str(c) for c in row) + " |")
    _p(lines)


def _score_cell(package: dict) -> str:
    """The overall score, or why there is none. Never a blank cell."""
    card = package.get("scorecard")
    if not card:
        return "not scored"
    return "null" if card["overall"] is None else f"{card['overall']:.4f}"


def _fraction_cell(package: dict) -> str:
    """How much of the applicable frame the score rests on."""
    card = package.get("scorecard")
    return "\u2014" if not card else f"{card['scored_fraction']:.3f}"


def render(package: dict, gaps: dict, provenance: dict, status: str) -> str:
    """One candidate's package as a document."""
    gene = package["gene"]
    ranking = package["ranking"]
    construct = package["construct"]
    target = package["target_evidence"]
    lines: list[str] = []

    _p(lines, f"# {gene} — candidate package")
    _p(lines)
    _p(lines, f"**{ranking['candidate_id']}** · `{package['accession']}` · "
              f"decision **{ranking['decision']}** · {ranking['gate_status']} · "
              f"{package['design_class']} · "
              f"{'on the Pareto front' if ranking['on_pareto_front'] else 'dominated'}"
              f" · position {ranking['position']} of {ranking['of']} · {status}")
    _p(lines)
    _p(lines, "This package carries what the pipeline produced. Eight of the "
              "reference document's twelve deliverables have something to "
              "carry; what the other four are missing is named in **What this "
              "package cannot tell you**, at the end, rather than left out.")
    _p(lines)
    _p(lines, "---")
    _p(lines)

    _p(lines, "## 1 — Ranking")
    _p(lines)
    _table(lines, ["field", "value"], [
        ["candidate_id", ranking["candidate_id"]],
        ["gate_status", ranking["gate_status"]],
        ["decision", ranking["decision"]],
        ["on Pareto front", "yes" if ranking["on_pareto_front"] else "no"],
        ["position", f"{ranking['position']} of {ranking['of']}"],
        ["position basis", ranking["position_basis"]],
    ])
    _p(lines)
    _table(lines, ["objective", "value"],
           [[k, v] for k, v in ranking["objectives"].items()])
    for reason in ranking["reasons"]:
        _p(lines, f"> {reason}")
        _p(lines)

    card = package.get("scorecard")
    if card:
        _p(lines, "## 2 \u2014 Scorecard")
        _p(lines)
        _table(lines, ["component", "weight", "state", "value", "source"],
               [[c["component"], c["weight"], c["state"],
                 "\u2014" if c["value"] is None else round(c["value"], 4),
                 c["source"]] for c in card["components"]])
        _p(lines)
        _table(lines, ["", ""], [
            ["weight version", card["weight_version"]],
            ["applicable weight", card["applicable_weight"]],
            ["measured weight", card["measured_weight"]],
            ["scored fraction", card["scored_fraction"]],
            ["floor", card["minimum_scored_fraction"]],
            ["evidence confidence", card["evidence_confidence"]],
            ["prediction uncertainty",
             "UNKNOWN" if card["prediction_uncertainty"] is None
             else card["prediction_uncertainty"]],
            ["confidence adjustment", card["confidence_adjustment"]],
            ["overall score",
             "not emitted" if card["overall"] is None else card["overall"]],
        ])
        _p(lines)
        if card["unknown_components"]:
            _p(lines, "**UNKNOWN on this candidate:** "
                      + ", ".join(card["unknown_components"])
                      + ". Each is named above with the reason it is missing. "
                        "None is imputed.")
            _p(lines)
        if card["not_applicable_components"]:
            _p(lines, "**NOT_APPLICABLE on this candidate:** "
                      + ", ".join(card["not_applicable_components"])
                      + ". This is a question that does not arise for this "
                        "design, not a gap in the evidence.")
            _p(lines)
        for reason in card["reasons"]:
            _p(lines, f"> {reason}")
            _p(lines)

    _p(lines, "## 3 — Construct")
    _p(lines)
    _table(lines, ["", ""], [
        ["architecture", construct["architecture"]],
        ["verdict", construct["verdict"]],
        ["length", f"{construct['total_bp']} bp of a "
                   f"{construct['budget_bp']} bp payload budget "
                   f"({construct['headroom_bp']} bp spare)"],
        ["residues", construct["residues"]],
        ["safety switch", "present" if construct["safety_switch"] else "absent"],
        ["binding domain", construct["binder_name"] or "none"],
    ])
    _p(lines, "### Domain map")
    _p(lines)
    _table(lines, ["domain", "residues", "aa", "bp", "provenance", "source"],
           [[d["name"], d["residues"], f"{d['aa_start']}-{d['aa_end']}",
             f"{d['bp_start']}-{d['bp_end']}", d["provenance"],
             (f"{d['accession']} {d['source_residues']}"
              if d["accession"] and d["source_residues"]
              else d["accession"] or "synthetic, named literal")]
            for d in construct["domains"]])
    if construct["reason"]:
        _p(lines, f"> {construct['reason']}")
        _p(lines)
    _p(lines, "The amino-acid sequence and the nucleotide map are in "
              f"`{gene}.json` beside this file. The DNA is a map under one "
              "fixed codon per residue, so the boundaries above are exact. It "
              "is not a codon-optimised ordering sequence.")
    _p(lines)

    _p(lines, "## 4 — Target evidence")
    _p(lines)
    _table(lines, ["", ""], [
        ["composite", target["composite"]],
        ["measured weight", target["measured_weight"]],
        ["evidence class", target["evidence_class"]],
        ["confidence", target["confidence"]],
        ["normal-tissue risk", f"{target['risk']} ({target['risk_organ']})"],
        ["risk basis", target["risk_basis"]],
        ["risk is a lower bound", target["risk_is_lower_bound"]],
        ["tumour-side verdict", target["tumour_side_verdict"]],
    ])
    _p(lines, "### Where the risk came from")
    _p(lines)
    attribution = target["risk_attribution"]
    _p(lines, f"Risk {attribution['risk']} on "
              f"{', '.join(attribution['winning_organs'])}, ahead of the next "
              f"organ by {attribution['margin']}, across "
              f"{attribution['organs_scored']} organs that scored.")
    _p(lines)
    _table(lines, ["organ", "weighted", "score", "tier", "arm", "staining",
                   "transcript"],
           [[o["organ"], round(o["weighted"], 4), round(o["score"], 4),
             o["tier"], o["arm"],
             (stage3.NOT_MEASURED if o["staining"] == stage3.NOT_MEASURED
              else f"{o['staining']['label']} {o['staining']['level_name']} "
                   f"{o['staining']['score']:.4f}"),
             (stage3.NOT_MEASURED if o["baseline"] == stage3.NOT_MEASURED
              else f"{o['baseline']['label']} {o['baseline']['tpm']:.1f} TPM "
                   f"{o['baseline']['score']:.4f}")]
            for o in attribution["organs"]])
    _p(lines, "Every score above recomputes from the measurement beside it, "
              "and the largest weighted value is the risk. `NOT_MEASURED` is a "
              "third state: it is not a zero and not a clean result.")
    _p(lines)

    _p(lines, "## 5 — Binders")
    _p(lines)
    binders = package["binders"]
    _p(lines, f"Stage 5 verdict: **{binders['verdict'] or 'no record'}** · "
              f"{binders['entries']} structural entries examined")
    _p(lines)
    _table(lines, ["route", "identifier", "name", "format", "stage", "affinity"],
           [[c["route"], c["identifier"], c["name"] or "-", c["format"] or "-",
             c["clinical_stage"] or "-", c["affinity"]]
            for c in binders["structure_route"] + binders["sequence_route"]])
    for reason in binders["reasons"]:
        _p(lines, f"> {reason}")
        _p(lines)

    _p(lines, "## 6 — Safety")
    _p(lines)
    safety = package["safety"]
    _table(lines, ["", ""], [
        ["verdict", safety["verdict"]],
        ["risk against the applied ceiling",
         f"{safety['risk']} against {safety['ceiling_applied']}"],
        ["peak organ", safety["risk_organ"]],
        ["binder origin", safety["binder_origin"]],
        ["source organism", safety["binder_source_organism"] or "not established"],
        ["epitope immunogenicity", safety["epitope_immunogenicity"]],
        ["trials naming this symbol",
         f"{safety['trials_total']}, {safety['trials_stopped']} stopped"],
    ])
    for reason in safety["reasons"]:
        _p(lines, f"- {reason}")
    _p(lines)

    _p(lines, "## 7 — Developability")
    _p(lines)
    dev = package["developability"]
    _p(lines, f"{dev['binders_scored']} binder sequence(s) scored.")
    _p(lines)
    _table(lines, ["binder", "residues", "pI", "net charge", "Cys", "sequons",
                   "APR", "flags"],
           [[r["binder"], r["residues"], r["isoelectric_point"],
             r["net_charge"], f"{r['cysteines']} ({r['cysteine_parity']})",
             len(r["glycosylation_sequons"]),
             len(r["aggregation_prone_starts"]), r["flag_count"]]
            for r in dev["rows"]])
    for reason in dev["reasons"]:
        _p(lines, f"> {reason}")
        _p(lines)

    _p(lines, "## 8 — Experimental validation plan")
    _p(lines)
    plan = package["validation_plan"]
    _p(lines, "### Before any bench work")
    _p(lines)
    for blocker in plan["before_any_bench_work"]:
        _p(lines, f"- {blocker}")
    _p(lines)
    for label, steps, setting in (("In vitro", plan["in_vitro"], "material"),
                                  ("In vivo", plan["in_vivo"], "model")):
        _p(lines, f"### {label}")
        _p(lines)
        _table(lines, ["step", "purpose", setting, "readout", "measures",
                       "acceptance"],
               [[s["step"], s["purpose"],
                 " · ".join(x for x in (s.get(setting), s.get("arms")) if x),
                 s["readout"], s["measures"], s["acceptance"]]
                for s in steps])
    for reason in plan["reasons"]:
        _p(lines, f"> {reason}")
        _p(lines)

    _p(lines, "## 9 — Provenance")
    _p(lines)
    _table(lines, ["source", "release", "role"],
           [[s["source"], s["release"], s["role"]]
            for s in provenance["sources"]])
    _p(lines, "Configuration hash chain, each covering the stage before it:")
    _p(lines)
    _table(lines, ["stage", "hash"],
           [[k, f"`{v}`"] for k, v in provenance["hashes"].items()])

    _p(lines, "---")
    _p(lines)
    _p(lines, "## What this package cannot tell you")
    _p(lines)
    _p(lines, f"{gaps['elements_missing']} elements the reference document "
              f"asks for are not produced, across "
              f"{gaps['deliverables_with_gaps']} deliverables. "
              f"{gaps['elements_probed']} are checked mechanically by the "
              f"verifier, {gaps['elements_measured']} recomputed from this "
              f"run, and {gaps['elements_asserted']} are judgements that say "
              "so.")
    _p(lines)
    for entry in gaps["by_deliverable"]:
        _p(lines, f"### Deliverable {entry['deliverable']} — {entry['title']} "
                  f"({entry['state']})")
        _p(lines)
        for element in entry["missing"]:
            _p(lines, f"**{element['element']}**")
            _p(lines)
            _p(lines, element["reason"])
            _p(lines)
            if element["blocking_stage"]:
                _p(lines, f"*Blocked by:* {element['blocking_stage']}")
                _p(lines)
            if element["note"]:
                _p(lines, f"*Note:* {element['note']}")
                _p(lines)
    return "\n".join(lines) + "\n"


def main() -> int:
    """Write one package per surviving candidate."""
    from car_pipeline.api import pipeline

    print("running the pipeline", flush=True)
    run = pipeline.run("Pancreatic Ductal Adenocarcinoma",
                       progress=lambda s, n="": None)
    packages = run["packages"]
    status = run["package_status"]
    gaps = stage12.gap_payload(packages)
    provenance = run["provenance"]

    if not packages:
        print(f"  {status}: no candidate reached the end, so no package was "
              "assembled. That is the result, not an empty directory.")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    for package in packages:
        gene = package["gene"]
        document = OUT / f"{gene}.md"
        document.write_text(render(package, gaps, provenance, status),
                            encoding="utf-8")
        payload = OUT / f"{gene}.json"
        payload.write_text(json.dumps(package, indent=2), encoding="utf-8")
        print(f"  {gene:10s} {len(document.read_text(encoding='utf-8')):6d} "
              f"chars  {document.relative_to(OUT.parent.parent).as_posix()}")

    index = OUT / "README.md"
    lines = ["# Candidate packages", "",
             f"{len(packages)} candidate(s), status {status}.", "",
             "| candidate | id | decision | score | fraction | class | construct | safety |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for p in packages:
        lines.append(
            f"| [{p['gene']}]({p['gene']}.md) | {p['ranking']['candidate_id']} "
            f"| {p['ranking']['decision']} "
            f"| {_score_cell(p)} | {_fraction_cell(p)} | {p['design_class']} | "
            f"{p['construct']['total_bp']} bp | {p['safety']['verdict']} |")
    lines += ["", f"Each package names what it cannot tell you: "
                  f"{gaps['elements_missing']} elements across "
                  f"{gaps['deliverables_with_gaps']} deliverables.", ""]
    index.write_text("\n".join(lines), encoding="utf-8")
    print(f"  index      {index.relative_to(OUT.parent.parent).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
