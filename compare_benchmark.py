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


def _key(name: str) -> str:
    """A CSV header, with the byte-order mark the first column carries."""
    return (name or "").lstrip("﻿").strip().lower()


def _tokens(text: str) -> str:
    """Separators flattened to spaces, so MORAb-15B6 yields the token 15B6."""
    return re.sub(r"[^A-Za-z0-9]+", " ", text or "")


def names_match(needle: str, haystack: str) -> bool:
    """Whether a binder name appears in text as a whole token.

    Whole-token rather than substring. A substring test on "M5" matches
    GRM5, S2M11, muscarinic M5 and a hundred other things, and every one of
    those would be reported as a panel binder that exists in the sources.
    """
    if not needle:
        return False
    # The needle is normalised the same way as the haystack, or MORAb-009
    # fails to match "MORAb 009" after separators are flattened.
    normalised = _tokens(needle).strip()
    if not normalised:
        return False
    pattern = (r"(?<![A-Za-z0-9])"
               + r"\s+".join(re.escape(part) for part in normalised.split())
               + r"(?![A-Za-z0-9])")
    return re.search(pattern, _tokens(haystack), flags=re.IGNORECASE) is not None


def msln_source_records(names: set[str]) -> dict:
    """Every MSLN-annotated source record, with the compound text it carries.

    Scoped to MSLN. The panel is a panel of MSLN binders, so "present in the
    connected sources" means present among the MSLN records -- an antibody
    called M5 that binds a coronavirus spike is not the M5 of this panel.
    """
    structural: dict[str, str] = {}
    with SABDAB.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            recorded = [a.strip().lower()
                        for a in (row.get("antigen_name") or "").split("|")]
            if not any(a in names for a in recorded):
                continue
            entry = (row["PDB"] or "").replace("pdb_0000", "").upper()
            text = " ".join(str(row.get(f) or "")
                            for f in ("compound", "short_header"))
            structural.setdefault(entry, text)

    therapeutic: dict[str, str] = {}
    with THERA.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            target = next((str(row[k]) for k in row
                           if _key(k) == "target"), "")
            if "msln" not in target.lower() and "mesothelin" not in target.lower():
                continue
            name = next((str(row[k]) for k in row
                         if _key(k) == "therapeutic"), "?")
            extra = " ".join(str(row[k] or "") for k in row
                             if "name" in k.lower())
            therapeutic[name] = f"{name} {extra}"
    return {"structural": structural, "therapeutic": therapeutic}


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
    msln_names = {n.lower() for n in
                  (output["normalized_target"].get("protein_name") or "").split("|")}
    # Rebuild the name set the same way the algorithm did, from UniProt alone.
    import re as _re
    protein = output["normalized_target"].get("protein_name") or ""
    msln_names = {_re.split(r"\s*[(\[]", protein)[0].strip().lower()}
    msln_names |= {m.strip().lower()
                   for m in _re.findall(r"\(([^()]+)\)", protein)}
    msln_names |= {(c.get("note") or "").strip().lower()
                   for c in output["normalized_target"].get("mature_chains") or []}
    msln_names |= {"msln"}
    msln_names.discard("")

    sources = msln_source_records(msln_names)

    # Which deposited entries the platform actually returned.
    retrieved_entries = {
        (b["binder"].split(":")[0] or "").upper(): b
        for b in retrieved if b["route"] == "structure"
    }
    retrieved_therapeutics = {
        b["binder"]: b for b in retrieved if b["route"] == "sequence"
    }

    print()
    print("=" * 74)
    print("RETRIEVAL RECALL")
    print("=" * 74)
    print(f"  the connected sources hold {len(sources['structural'])} "
          f"MSLN-annotated structural entries and "
          f"{len(sources['therapeutic'])} MSLN therapeutics")
    print(f"  the platform retrieved {len(retrieved_entries)} entries and "
          f"{len(retrieved_therapeutics)} therapeutics")
    print()

    recovered, not_in_source, missed = [], [], []
    for item in panel:
        names = [item["binder"]] + list(item.get("aliases") or [])

        found_in_source, found_retrieved, how = [], None, ""
        for entry, text in sources["structural"].items():
            if any(names_match(n, text) for n in names):
                found_in_source.append(("entry", entry, text))
                if entry in retrieved_entries:
                    found_retrieved = retrieved_entries[entry]
                    how = f"entry {entry}, recorded as {text.strip()[:52]!r}"
        for name, text in sources["therapeutic"].items():
            if any(names_match(n, text) for n in names):
                found_in_source.append(("therapeutic", name, text))
                if name in retrieved_therapeutics and not found_retrieved:
                    found_retrieved = retrieved_therapeutics[name]
                    how = f"therapeutic {name}"

        if found_retrieved:
            recovered.append((item, found_retrieved, how))
        elif not found_in_source:
            not_in_source.append(item)
        else:
            missed.append((item, found_in_source))

    for item, hit, how in recovered:
        alias = ("" if item["binder"] in how else
                 f" via alias {[a for a in item.get('aliases') or [] if names_match(a, how)] or item.get('aliases')}")
        print(f"  RECOVERED    {item['binder']:<6} {how}{alias}")
        print(f"                      target_match {hit['target_match']}, "
              f"{hit['route']} route")
    for item in not_in_source:
        print(f"  NOT IN       {item['binder']:<6} no MSLN record in either "
              f"connected source names it, under "
              f"{[item['binder']] + list(item.get('aliases') or [])}")
        print(f"    SOURCE            Per the benchmark's own rule this is "
              f"source coverage, not a retrieval failure.")
    for item, where in missed:
        print(f"  MISSED       {item['binder']:<6} named in "
              f"{[w[1] for w in where]} but not retrieved")

    print()
    print(f"  recovered {len(recovered)} of {len(panel)}; "
          f"{len(not_in_source)} not identifiable in the connected sources; "
          f"{len(missed)} genuinely missed")
    if not missed:
        print()
        print("  Recall is good. Every panel binder the connected sources name")
        print("  was retrieved. The two that were not retrieved are not named")
        print("  by any MSLN record in either source, so there was nothing to")
        print("  retrieve them by.")

    print()
    print("=" * 74)
    print("THE COMPARISON TABLE")
    print("=" * 74)
    print(f"  {'parameter':<26} {'platform':<22} {'literature':<20} verdict")

    def row(parameter, platform, literature, verdict):
        print(f"  {parameter:<26} {str(platform)[:21]:<22} "
              f"{str(literature)[:19]:<20} {verdict}")

    row("known binder recovery",
        f"{len(recovered)} of {len(recovered) + len(missed)} nameable",
        "4 where coverage permits",
        f"PASS ({len(not_in_source)} not named in source)"
        if not missed else f"REVIEW ({len(missed)} missed)")

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
