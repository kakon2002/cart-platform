"""B3 — does the recorded antigen name the target we asked about?

This asks a deliberately weaker question than the one the platform refused
elsewhere. It does not ask what an antibody *binds*: that needs coordinates,
and no coordinate file is connected. It asks what the depositors *recorded* as
its antigen, which is annotation, and which is enough to catch the failure this
check exists for — an entry found by searching on the target's accession whose
antibody is against a different chain of the same complex.

The two claims are not the same and the output says so. A PASS means the entry's
own annotation names this target. It does not mean the interface has been
inspected.
"""

from __future__ import annotations

import re

PASS = "PASS"
FAIL = "FAIL"
UNKNOWN = "UNKNOWN"

VERDICTS = (PASS, FAIL, UNKNOWN)

# Why a row reached its verdict. The verdict alone is not enough to audit the
# counting rule below: three unrelated situations all read UNKNOWN, and a
# reader needs to see which one admitted a binder to the count. Parsing the
# reason prose for that would make the rule depend on wording.
MATCHED = "matched"
WRONG_ANTIGEN = "wrong_antigen"
NO_ANNOTATION = "no_annotation"
SEQUENCE_ROUTE = "sequence_route"
NO_RECORD = "no_record"

PATHS = (MATCHED, WRONG_ANTIGEN, NO_ANNOTATION, SEQUENCE_ROUTE, NO_RECORD)

# The paths that admit a binder to the objective. One clause: count a binder
# unless its recorded antigen names something else.
COUNTED_PATHS = (MATCHED, NO_ANNOTATION, SEQUENCE_ROUTE, NO_RECORD)

# Values the source uses for "nothing recorded". They are absence, not a name,
# and they must not be matched against or treated as a mismatch.
ABSENT = {"", "na", "n/a", "none", "null", "unknown", "-"}

# Recorded antigens that are laboratory reagents rather than the subject of the
# experiment. They appear beside a real antigen in a pipe-delimited list and
# must not make an otherwise-matching entry fail.
REAGENT_WORDS = (
    "ion", "sulfate", "sulphate", "phosphate", "acetate", "chloride",
    "glycerol", "water", "buffer", "detergent", "peg", "polyethylene glycol",
)


def name_set(record) -> set[str]:
    """Every name UniProt gives this target, derived rather than written down.

    The recommended name, each parenthesised synonym, and the note on every
    mature chain. The chain notes matter: an entry may record the antigen as a
    processed form -- "Mesothelin, cleaved form" -- which is the target under a
    name the recommended name alone does not contain.
    """
    names: set[str] = set()
    protein = getattr(record, "protein_name", "") or ""
    head = re.split(r"\s*[\(\[]", protein)[0].strip()
    if head:
        names.add(head)
    names |= {m.strip() for m in re.findall(r"\(([^()]+)\)", protein) if m.strip()}
    for chain in getattr(record, "chains", None) or []:
        note = (getattr(chain, "note", "") or "").strip()
        if note:
            names.add(note)
    gene = (getattr(record, "gene", "") or "").strip()
    if gene:
        names.add(gene)
    return {n for n in names if n}


def antigen_elements(antigen_name: str | None) -> list[str]:
    """The recorded antigens, one per element.

    The source records more than one antigen as a pipe-delimited list, and the
    elements are different molecules. Matching the whole string would fail an
    entry whose antigen is the target plus a crystallisation additive, and it
    would also let a two-protein complex match on neither while looking like a
    single unrecognised name.
    """
    return [part.strip() for part in (antigen_name or "").split("|")
            if part.strip()]


def _is_absent(element: str) -> bool:
    """Whether this element records nothing rather than something."""
    return element.strip().lower() in ABSENT


def _is_reagent(element: str) -> bool:
    """Whether this element is a laboratory reagent rather than an antigen."""
    low = element.strip().lower()
    return any(word in low for word in REAGENT_WORDS)


def evaluate(antigen_name: str | None, names: set[str],
             antigen_type: str | None = None,
             antigen_species: str | None = None) -> dict:
    """One target-match verdict, with the evidence it rests on.

    Species is carried and never gates. It is recorded as NA on entries that
    are genuinely the target, so requiring it would fail those, and a check
    that fails correct records is worse than no check.
    """
    elements = antigen_elements(antigen_name)
    lowered = {n.strip().lower() for n in names}

    informative = [e for e in elements if not _is_absent(e)]
    matched = [e for e in informative if e.strip().lower() in lowered]
    reagents = [e for e in informative if _is_reagent(e)]
    others = [e for e in informative if e not in matched and e not in reagents]

    payload = {
        "target_match_score": None,
        "wrong_antigen_flag": None,
        "recorded_antigens": elements,
        "matched_antigens": matched,
        "unmatched_antigens": others,
        "reagent_antigens": reagents,
        "antigen_type": antigen_type or None,
        "antigen_species": antigen_species or None,
        "species_is_a_gate": False,
    }

    if not informative:
        return {
            **payload,
            "target_match": UNKNOWN,
            "match_path": NO_ANNOTATION,
            "reason": (
                "The entry records no antigen name, so there is nothing to "
                "match against. This is absence of annotation, not evidence "
                "of a different antigen: it reads UNKNOWN rather than FAIL, "
                "and settling it would need the coordinates that are not "
                "connected."),
        }

    if matched:
        return {
            **payload,
            "target_match": PASS,
            "match_path": MATCHED,
            "target_match_score": 1.0,
            "wrong_antigen_flag": False,
            "reason": (
                f"The entry's own annotation names this target: "
                f"{', '.join(matched)}. That is a claim about what the "
                f"depositors recorded, not about what the antibody was "
                f"observed to contact"
                + (f"; {', '.join(reagents)} recorded alongside is a "
                   "laboratory reagent and is disregarded" if reagents else "")
                + "."),
        }

    return {
        **payload,
        "target_match": FAIL,
        "match_path": WRONG_ANTIGEN,
        "target_match_score": 0.0,
        "wrong_antigen_flag": True,
        "reason": (
            f"The entry records its antigen as {', '.join(others)}, which is "
            f"not this target under any name UniProt gives it. The entry was "
            f"found by searching on the target's accession, so it contains "
            f"the target somewhere; the antibody in it is annotated against "
            f"something else."),
    }


def check(rows: list[dict], record) -> list[dict]:
    """Apply the target-match check to every retrieved binder row.

    Sequence-route binders carry no deposited antigen annotation at all, so
    they read UNKNOWN for want of the evidence this check reads rather than
    because anything about them is doubtful.
    """
    names = name_set(record) if record is not None else set()
    out = []
    for row in rows:
        if record is None:
            verdict = {
                "target_match": UNKNOWN,
                "match_path": NO_RECORD,
                "target_match_score": None,
                "wrong_antigen_flag": None,
                "recorded_antigens": [],
                "matched_antigens": [],
                "unmatched_antigens": [],
                "reagent_antigens": [],
                "antigen_type": None,
                "antigen_species": None,
                "species_is_a_gate": False,
                "reason": ("No UniProt record for this target, so no name set "
                           "could be derived to match against."),
            }
        elif row.get("route") == "sequence":
            verdict = {
                "target_match": UNKNOWN,
                "match_path": SEQUENCE_ROUTE,
                "target_match_score": None,
                "wrong_antigen_flag": None,
                "recorded_antigens": [],
                "matched_antigens": [],
                "unmatched_antigens": [],
                "reagent_antigens": [],
                "antigen_type": None,
                "antigen_species": None,
                "species_is_a_gate": False,
                "reason": (
                    "This binder comes from the named-therapeutic route, which "
                    "records the target it was raised against but no deposited "
                    "antigen annotation. There is no recorded antigen to match, "
                    "so the check that catches a wrong-antigen structural hit "
                    "does not apply and reads UNKNOWN rather than PASS."),
            }
        else:
            verdict = evaluate(
                row.get("antigen_name"), names,
                row.get("antigen_type"), row.get("antigen_species"))
        out.append({**row, **verdict})
    return out


def counts_towards_objective(row: dict) -> bool:
    """Whether one checked binder contributes to the ranking objective.

    One clause: count it unless its recorded antigen names something else.

    Written as a single negative test on purpose. The positive form -- count a
    binder whose annotation names the target -- reads more natural and is a
    different rule: it would drop every sequence-route binder, whose UNKNOWN is
    about a check that does not apply rather than about doubt, and it would drop
    entries with no annotation at all, treating absence as evidence of a
    different antigen. Both are imputations against the binder.
    """
    return row.get("target_match") != FAIL


def objective_counts(rows: list[dict]) -> dict:
    """The three counts and the path breakdown, over checked binder rows.

    The counts are arithmetically closed -- retrieved equals counted plus
    wrong-antigen -- so a reader can see that no binder was dropped between the
    number the search returned and the number the ranking used. The breakdown
    says which path admitted each counted binder, because three different
    situations read UNKNOWN and a count resting mostly on absent annotation is a
    weaker claim than one resting on annotation that matched.
    """
    counted = [r for r in rows if counts_towards_objective(r)]
    paths = {p: 0 for p in PATHS}
    for row in rows:
        path = row.get("match_path")
        if path in paths:
            paths[path] += 1
    return {
        "retrieved": len(rows),
        "counted": len(counted),
        "wrong_antigen": len(rows) - len(counted),
        "paths": paths,
        "counted_by_path": {p: paths[p] for p in COUNTED_PATHS},
    }


NOTES = [
    "The ranking objective counts a binder unless its recorded antigen names "
    "something else. A binder whose entry records no antigen still counts: "
    "absent annotation is not evidence of a different antigen, and refusing to "
    "count it would impute against the binder. The counted set is broken down "
    "by path so a reader can see which of those admitted it.",
    "A target-match verdict is a claim about the annotation a depositor "
    "recorded, not about an interface anyone inspected. Deciding what an "
    "antibody contacts needs coordinates, and no coordinate file is connected.",
    "Recorded antigens are matched one element at a time. The source records "
    "several antigens as a pipe-delimited list, and whole-string matching "
    "would fail an entry whose antigen is the target beside a crystallisation "
    "additive.",
    "An entry with no recorded antigen reads UNKNOWN, never PASS and never "
    "FAIL. Absent annotation is not evidence of a different antigen.",
    "Antigen species is carried and never gates. Entries that are genuinely "
    "the target record it as NA, so requiring it would fail correct records.",
]
