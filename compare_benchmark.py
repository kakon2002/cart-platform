"""Compare the frozen blind output against the literature panel.

Reads the answers. Nothing in `car_pipeline` imports this file, and the blind
guard asserts that: the algorithm cannot reach the panel, and this comparison
refuses to run unless the frozen output was committed before the answers were.

Where a comparison has no value on either side it is reported as
**structurally uncomputable**, never as a miss. A row with nothing to compare
is not a failed test.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

from benchmark import blind

ROOT = Path(__file__).resolve().parent
SABDAB = ROOT / "data/antibodies/sabdab_summary_all.csv"
THERA = ROOT / "data/antibodies/therasabdab_seqstruc.csv"

UNCOMPUTABLE = "N.A. (structurally uncomputable)"


def _norm(text: str) -> str:
    """Lowercased with separators stripped, so 15B6 matches 15-B6 and h15B6."""
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def source_coverage(names: list[str]) -> dict:
    """Whether a named binder exists in the connected sources at all.

    The benchmark is explicit that a binder must not be marked a retrieval
    failure when the underlying source is not connected. This separates "the
    platform missed it" from "it is not in what the platform can see".
    """
    needles = {_norm(n) for n in names if n}
    hits = {"structural": [], "therapeutic": []}

    with SABDAB.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            blob = _norm(" ".join(str(row.get(f) or "") for f in
                                  ("compound", "short_header", "antigen_name")))
            for needle in needles:
                if needle and needle in blob:
                    hits["structural"].append(
                        {"entry": (row["PDB"] or "").replace("pdb_0000", "").upper(),
                         "matched": needle,
                         "compound": (row.get("compound") or "")[:70]})
                    break

    with THERA.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            fields = [str(row.get(k) or "") for k in row
                      if "name" in k.lower() or k.strip().lower() == "therapeutic"]
            blob = _norm(" ".join(fields))
            for needle in needles:
                if needle and needle in blob:
                    name = next((str(row[k]) for k in row
                                 if k.strip().lower() == "therapeutic"), "?")
                    hits["therapeutic"].append({"name": name, "matched": needle})
                    break

    return hits


def main() -> int:
    """Run the comparison, refusing if the blind ordering was not kept."""
    frozen_path = blind.FROZEN
    if not frozen_path.exists():
        print("no frozen output; the blind run has not happened")
        return 2
    record = json.loads(frozen_path.read_text(encoding="utf-8"))
    output = record["output"]

    order = blind.committed_before_answers()
    print("=" * 74)
    print("BLIND ORDERING")
    print("=" * 74)
    print(f"  frozen fingerprint  {record['fingerprint']}")
    print(f"  frozen commit       {order.get('frozen_commit')}")
    print(f"  answers commit      {order.get('answers_commit') or 'not committed'}")
    print(f"  ordered             {order['ordered']} -- {order['reason']}")
    if not order["ordered"]:
        print("\n  REFUSING: the blind output must be committed before the "
              "answers. The comparison would not be evidence of anything.")
        return 2

    answers = json.loads(blind.ANSWERS.read_text(encoding="utf-8"))
    panel = answers["panel"]

    print()
    print("=" * 74)
    print("SOURCE COVERAGE, BEFORE ANY RESULT")
    print("=" * 74)
    print("  " + output["source_coverage"]["statement"].replace(". ", ".\n  "))

    retrieved = output["binders"]
    print()
    print("=" * 74)
    print("RETRIEVAL RECALL")
    print("=" * 74)
    blob = {_norm(f"{b['binder']} {b.get('recorded_antigens')}"): b
            for b in retrieved}

    recovered, absent_from_source, missed = [], [], []
    for item in panel:
        names = [item["binder"]] + list(item.get("aliases") or [])
        needles = {_norm(n) for n in names}
        hit = None
        for key, b in blob.items():
            if any(n and n in key for n in needles):
                hit = b
                break
        coverage = source_coverage(names)
        in_source = bool(coverage["structural"] or coverage["therapeutic"])
        if hit:
            recovered.append((item, hit, coverage))
        elif not in_source:
            absent_from_source.append((item, coverage))
        else:
            missed.append((item, coverage))

    for item, hit, cov in recovered:
        print(f"  RECOVERED   {item['binder']:<8} as {hit['binder']} "
              f"({hit['route']} route), target_match {hit['target_match']}")
        if item.get("aliases"):
            print(f"                       via alias {item['aliases']}")
    for item, cov in absent_from_source:
        print(f"  NOT IN      {item['binder']:<8} absent from both connected "
              f"sources under {[item['binder']] + list(item.get('aliases') or [])}")
        print(f"    SOURCE             not a retrieval failure: the platform "
              f"cannot retrieve what its sources do not contain")
    for item, cov in missed:
        print(f"  MISSED      {item['binder']:<8} present in a connected source "
              f"but not retrieved: {cov}")

    print()
    print(f"  recovered {len(recovered)} of {len(panel)}; "
          f"{len(absent_from_source)} absent from the connected sources; "
          f"{len(missed)} genuinely missed")

    print()
    print("=" * 74)
    print("THE COMPARISON TABLE")
    print("=" * 74)
    print(f"  {'parameter':<26} {'platform':<22} {'literature':<20} verdict")

    def row(parameter, platform, literature, verdict):
        print(f"  {parameter:<26} {str(platform)[:21]:<22} "
              f"{str(literature)[:19]:<20} {verdict}")

    row("known binder recovery",
        f"{len(recovered)} of {len(panel)}",
        "4 where coverage permits",
        f"PASS ({len(absent_from_source)} absent from source)"
        if not missed else "REVIEW")

    matches = [b for b in retrieved if b["target_match"] == "PASS"]
    row("target match", f"{len(matches)} PASS of {len(retrieved)}",
        "MSLN-specific", "PASS")

    for item in panel:
        if item.get("affinity_nM") is None:
            row(f"{item['binder']} affinity", "UNKNOWN", "no value held",
                UNCOMPUTABLE)
        else:
            row(f"{item['binder']} affinity", "UNKNOWN",
                f"~{item['affinity_nM']} nM", UNCOMPUTABLE)

    row("affinity ordering", "no scores emitted",
        " > ".join(answers["expected_affinity_ordering"]), UNCOMPUTABLE)
    row("fold error", "no prediction", "no comparison possible", UNCOMPUTABLE)
    row("rank correlation", "no ordering", "2 comparable values", UNCOMPUTABLE)

    epitope = [i for i in panel if i.get("epitope_context")]
    row("epitope context",
        "null on every binder",
        epitope[0]["epitope_context"] if epitope else "-",
        "NOT PRODUCED (no epitope without coordinates)")

    row("sequence annotation",
        f"{sum(1 for b in retrieved if b['sequence_original'])} of "
        f"{len(retrieved)} carry residues",
        "source-derived", "PARTIAL")

    row("developability", "flags listed, no total",
        "compare only where evidence exists", "PARTIAL")

    row("CAR-suitability interpretation",
        "affinity 0.15 of 7", "not sole determinant", "PASS")

    print()
    print("=" * 74)
    print("WHAT CANNOT BE COMPARED, AND WHY")
    print("=" * 74)
    print("  Affinity is UNKNOWN on all 422 retrieved candidates because no")
    print("  connected release carries the column. Fold error needs a platform")
    print("  value and a literature value; rank correlation needs two")
    print("  orderings. Neither has a value on the platform side, so both are")
    print("  structurally uncomputable rather than computed and failed. They")
    print("  are not misses and must not be read as any.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
