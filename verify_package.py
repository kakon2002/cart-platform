"""Runs the candidate package and tests it against the criteria fixed in the spec."""

from __future__ import annotations

import importlib.util
import os
import sys

os.environ.setdefault("CART_NO_MATRIX_FETCH", "1")

from car_pipeline.stages import stage3, stage6, stage11, stage12, validation


def probe_open(kind: str, target: str, name: str, packages: list[dict]) -> bool:
    """Whether the gap this probe describes is still open."""
    if kind == stage12.PROBE_MODULE:
        return importlib.util.find_spec(f"car_pipeline.stages.{target}") is None
    if kind == stage12.PROBE_FIELD:
        module_path, _, attribute = target.rpartition(".")
        try:
            module = importlib.import_module(module_path)
        except ImportError:
            return True
        holder = getattr(module, attribute, None)
        if holder is None:
            return True
        return name not in getattr(holder, "__annotations__", {})
    if kind == stage12.PROBE_KEY:
        return all(name not in (p.get(target) or {}) for p in packages)
    raise ValueError(f"unknown probe kind {kind!r}")


def main() -> int:
    """Run the candidate-package criteria."""
    print("running the pipeline", flush=True)
    from car_pipeline.api import pipeline

    run = pipeline.run("Pancreatic Ductal Adenocarcinoma",
                       progress=lambda s, n="": None)
    packages = run["packages"]
    status = run["package_status"]
    survivors = [r for r in run["final"] if r.survived]
    gaps = stage12.gap_payload(packages)
    provenance = run["provenance"]
    print(f"  {len(packages)} package(s), status {status}, "
          f"{len(survivors)} survivor(s)")

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    checked: list[str] = []
    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    packaged = [p["gene"] for p in packages]
    expected = [r.gene for r in survivors]
    criterion(
        "Q1", not packages or packaged != expected,
        f"{len(packages)} package(s), one per surviving candidate, in the "
        f"ranking's order: {', '.join(packaged)}"
        if packages and packaged == expected else
        "no package was assembled, so every criterion below would read an "
        "empty set" if not packages else
        f"packaged {packaged} against survivors {expected}")

    survived = {r.gene for r in survivors}
    failed = {r.gene for r in run["final"] if not r.survived}
    q2_bad = ([f"{g}: packaged but did not reach the end"
               for g in packaged if g in failed]
              + [f"{g}: reached the end and was not packaged"
                 for g in survived if g not in set(packaged)])
    criterion("Q2", bool(q2_bad),
              f"every one of {len(survived)} candidates that reached the end is "
              f"packaged, and none of the {len(failed)} that did not"
              if not q2_bad else "; ".join(q2_bad[:3]))

    q3_bad = []
    for p in packages:
        for section in stage12.SECTIONS:
            if section not in p:
                q3_bad.append(f"{p['gene']}: no {section} section")
        if not p["construct"].get("amino_acid_sequence"):
            q3_bad.append(f"{p['gene']}: construct section carries no sequence")
        if p["target_evidence"].get("risk") is None:
            q3_bad.append(f"{p['gene']}: target evidence carries no risk")
        if not p["validation_plan"].get("in_vitro"):
            q3_bad.append(f"{p['gene']}: validation plan carries no steps")
        if not p["safety"].get("verdict"):
            q3_bad.append(f"{p['gene']}: safety section carries no verdict")
    criterion("Q3", bool(q3_bad) or not packages,
              f"all {len(stage12.SECTIONS)} sections present on every package, "
              f"and each carries what its stage produced"
              if packages and not q3_bad else
              "; ".join(q3_bad[:3]) if q3_bad else "no package to inspect")

    q4_bad = []
    for p in packages:
        c = p["construct"]
        if stage6.translate(c["dna"]) != c["amino_acid_sequence"]:
            q4_bad.append(f"{p['gene']}: packaged DNA does not translate back")
            continue
        pos = 0
        for d in c["domains"]:
            if d["aa_start"] != pos or d["aa_end"] > len(c["amino_acid_sequence"]):
                q4_bad.append(f"{p['gene']}: domain boundaries do not partition")
                break
            pos = d["aa_end"]
        else:
            if pos != len(c["amino_acid_sequence"]):
                q4_bad.append(f"{p['gene']}: domains stop short of the sequence")
    criterion("Q4", bool(q4_bad) or not packages,
              f"the packaged DNA translates to the packaged sequence and the "
              f"packaged domains partition it, for all {len(packages)}"
              if packages and not q4_bad else
              "; ".join(q4_bad[:3]) if q4_bad else "no construct to check")

    by_ranked = {r.gene: r for r in run["ranked"] if r.gene}
    q5_bad = []
    for p in packages:
        attribution = p["target_evidence"]["risk_attribution"]
        organs = attribution["organs"]
        if not organs:
            q5_bad.append(f"{p['gene']}: no organ attributed")
            continue
        top = max(o["weighted"] for o in organs)
        if abs(top - attribution["risk"]) > 1e-12:
            q5_bad.append(f"{p['gene']}: attribution {attribution['risk']} "
                          f"against max {top}")
        source = by_ranked[p["gene"]].risk
        if round(attribution["risk"], 4) != source:
            q5_bad.append(f"{p['gene']}: packaged risk rounds to "
                          f"{round(attribution['risk'], 4)}, Stage 3 says {source}")
        if p["target_evidence"]["risk"] != source:
            q5_bad.append(f"{p['gene']}: carried risk differs from Stage 3")
    criterion("Q5", bool(q5_bad) or not packages,
              f"every packaged attribution reconstructs its own risk to within "
              f"1e-12 and matches Stage 3"
              if packages and not q5_bad else
              "; ".join(q5_bad[:3]) if q5_bad else "no attribution to check")

    executed, q6_bad = 0, []
    for entry in gaps["by_deliverable"]:
        for element in entry["missing"]:
            probe = element["probe"]
            if not probe:
                continue
            executed += 1
            kind, target = probe[0], probe[1]
            name = probe[2] if len(probe) > 2 else ""
            if not probe_open(kind, target, name, packages):
                q6_bad.append(
                    f"deliverable {element['deliverable']}: "
                    f"{element['element'][:48]} is no longer missing "
                    f"({kind} {target}{'.' + name if name else ''} exists)")
    criterion(
        "Q6", bool(q6_bad) or executed == 0,
        f"{executed} declared gap(s) probed and all still open; "
        f"{gaps['elements_measured']} recomputed from the run and "
        f"{gaps['elements_asserted']} stated as judgements"
        if executed and not q6_bad else
        "; ".join(q6_bad[:3]) if q6_bad else
        "no gap carried a probe, so the whole table is unverified assertion")

    classes = validation.design_class_summary(run["constructs"])
    backup = classes["conservative_design"]
    refusal = [r for r in classes["reasons"] if "conservative" in r.lower()]
    has_counts = any(any(ch.isdigit() for ch in r) for r in refusal)
    labelled = [c.gene for c in run["constructs"]
                if validation.design_class(c) == validation.CONSERVATIVE]
    q7_bad = []
    if backup is None and not refusal:
        q7_bad.append("no conservative backup and no refusal stating why")
    if backup is None and refusal and not has_counts:
        q7_bad.append("the refusal carries no counts behind it")
    if backup is not None and sorted(backup) != sorted(labelled):
        q7_bad.append(f"reports {backup} as conservative against {labelled}")
    criterion(
        "Q7", bool(q7_bad),
        "no conservative design exists in this pool and the section says so "
        "with the counts behind it, rather than standing blank"
        if backup is None and not q7_bad else
        f"conservative backup: {backup}" if not q7_bad else
        "; ".join(q7_bad))

    pins = provenance["sources"]
    missing_pin = [s["source"] for s in pins if not s["release"]]
    chain = provenance["hashes"]
    wanted = ("stage3", "stage4", "stage5", "stage6", "stage9", "stage10",
              "stage11")
    absent = [k for k in wanted if not chain.get(k)]
    criterion(
        "Q8", bool(missing_pin) or bool(absent),
        f"{len(pins)} connected sources each name a release, and the hash "
        f"chain is unbroken from Stage 3 to Stage 11"
        if not (missing_pin or absent) else
        f"sources without a release {missing_pin}; chain gaps {absent}")

    forbidden = ("structural_report", "structure_report", "stage7",
                 "functional_predictions", "stage8", "functional_response")
    q9_bad = []
    for p in packages:
        for key in p:
            if key in forbidden:
                q9_bad.append(f"{p['gene']}: emits a {key} section")
        for section in ("construct", "safety", "target_evidence"):
            for key in p[section]:
                if key in forbidden:
                    q9_bad.append(f"{p['gene']}.{section}: emits {key}")
    absent_gaps = [e for entry in gaps["by_deliverable"] for e in entry["missing"]
                   if e["state"] == stage12.ABSENT]
    criterion(
        "Q9", bool(q9_bad) or len(absent_gaps) < 2,
        f"no package emits a section or placeholder for Stage 7 or Stage 8; "
        f"both are recorded in the gaps section instead ({len(absent_gaps)} "
        f"absent-stage entries)"
        if not q9_bad and len(absent_gaps) >= 2 else
        "; ".join(q9_bad[:3]) if q9_bad else
        f"only {len(absent_gaps)} absent stage(s) recorded, expected 2")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHAT A READER RECEIVES")
    print("=" * 72)
    print(f"    {len(packages)} candidate package(s), status {status}")
    for p in packages:
        r = p["ranking"]
        print(f"      {p['gene']:10s} rank {r['position']}/{r['of']}  "
              f"front={'yes' if r['on_pareto_front'] else 'no ':3s}  "
              f"{p['design_class']:9s}  {p['construct']['total_bp']:5d} bp  "
              f"{p['safety']['verdict']}")

    print()
    print("  WHAT THE PACKAGE CANNOT CARRY")
    print("  " + "-" * 68)
    print(f"    {gaps['elements_missing']} element(s) missing across "
          f"{gaps['deliverables_with_gaps']} deliverable(s)")
    print(f"    {gaps['elements_probed']} probed, "
          f"{gaps['elements_measured']} measured from the run, "
          f"{gaps['elements_asserted']} judgements")
    for entry in gaps["by_deliverable"]:
        print(f"      deliverable {entry['deliverable']:2d} {entry['state']:8s} "
              f"{len(entry['missing'])} element(s)  {entry['title']}")

    print()
    print("  THE TWO ABSENT STAGES")
    print("  " + "-" * 68)
    for element in absent_gaps:
        print(f"    Stage {element['blocking_stage']}: {element['note']}")

    print()
    print("  PROVENANCE")
    print("  " + "-" * 68)
    for source in pins:
        print(f"    {source['source']:26s} {str(source['release']):22s} "
              f"{source['role']}")
    print(f"    hash chain  " + " -> ".join(chain[k] for k in wanted))
    return 0


if __name__ == "__main__":
    sys.exit(main())
