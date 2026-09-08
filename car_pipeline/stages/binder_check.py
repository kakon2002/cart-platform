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
MATCHED_TOLERANT = "matched_tolerant"
WRONG_ANTIGEN = "wrong_antigen"
NO_ANNOTATION = "no_annotation"
SEQUENCE_ROUTE = "sequence_route"
NO_RECORD = "no_record"

PATHS = (MATCHED, MATCHED_TOLERANT, WRONG_ANTIGEN, NO_ANNOTATION,
         SEQUENCE_ROUTE, NO_RECORD)

# The paths that admit a binder to the objective. One clause: count a binder
# unless its recorded antigen names something else.
COUNTED_PATHS = (MATCHED, MATCHED_TOLERANT, NO_ANNOTATION, SEQUENCE_ROUTE,
                 NO_RECORD)

# What the matcher will accept, in the configuration hash. A run under a
# different matching rule counts different binders and must not compare equal
# to one under this rule.
MATCH_BASIS = (
    "exact name, or the target name appearing in order inside a longer label "
    "that adds only form qualifiers and short isoform designators; "
    "comma-joined fusion partners separated first")

# The surplus a tolerant match will tolerate, and nothing else.
#
# A blacklist of dangerous words was tried first and is not sufficient. The
# words that separate two proteins are ordinary content words: CELSR2 is
# recorded as "EGF-like protein 2" and CELSR3 as "Multiple EGF-like domains
# protein 2", so the first is carried inside the second with "multiple" and
# "domains" left over. Nothing marks those two words as dangerous except that
# they name a different protein.
#
# So the rule is inverted. A longer label matches only when everything it adds
# is a qualifier -- a word that describes which form of the target is present,
# not which protein it is. Anything else declines, which withholds credit from
# a real binder rather than inventing evidence for one.
QUALIFIER_WORDS = frozenset({
    "isoform", "isoforms", "variant", "variants", "form", "of", "the",
    "peptide", "peptides", "fragment", "fragments", "epitope",
    "construct", "inactive", "recombinant", "synthetic", "soluble",
    "mature", "precursor", "truncated", "tagged", "full", "length",
    "ecto", "ectodomain", "extracellular",
})

# Isoform designators -- the A2 of "Isoform A2 of Claudin-18", the CRA_a of a
# submitted-name isoform. They are labels rather than words, and they are
# admitted by length: a surplus word of one or two characters cannot be the
# name of a different protein.
#
# A bare number is never admitted, whatever its length. Numbering is what
# separates one family member from another -- APLP-1 from APLP-2, and the
# bestrophins from each other -- and the reference set carries bare stems that
# sit inside every member of such a family.
MAXIMUM_DESIGNATOR_LENGTH = 2


def _admissible_surplus(word: str) -> bool:
    """Whether one leftover word still lets a longer label name the target."""
    if word.isdigit():
        return False
    return word in QUALIFIER_WORDS or len(word) <= MAXIMUM_DESIGNATOR_LENGTH

# A reference name shorter than this, or carrying no letter, is not used for
# tolerant matching. The derived name set contains degenerate entries -- a bare
# "+" appears on 115 targets, pulled out of constructions such as
# Na(+)/K(+)-transporting ATPase, and bare single letters on others. Exact
# equality is accidentally immune to these because no recorded antigen is ever
# the single character "+". Containment is not: an empty or one-character token
# set is contained in nearly everything. The immunity was an accident of the
# old rule and this guard replaces it with an intended one.
MINIMUM_NAME_LENGTH = 4

# A reference name carrying more than one digit-only word is a structured
# identifier, not a protein name -- an enzyme classification code, a catalogue
# number. Containment compares word SETS, which is exactly what lets a
# word-order variant match, and is exactly what makes these unsafe: the codes
# 7.6.2.1 and 6.2.1.7 have identical word sets and name different enzymes.
# Order-insensitivity cannot be had for prose and refused for codes, so the
# codes are excluded instead.
MAXIMUM_NUMERIC_WORDS = 1

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


def name_tokens(text: str) -> frozenset[str]:
    """The lowercased alphanumeric words of a name.

    Nothing is discarded. Dropping words that look generic is what makes a
    tolerant matcher dangerous: drop "receptor" and the hepatocyte growth
    factor beta chain becomes a match for the hepatocyte growth factor
    receptor, which is a different molecule and the antibody against one is not
    evidence about the other.

    Splitting on alphanumeric runs is also what keeps the trailing numbers
    apart. Claudin-1 and Claudin-18 differ only in that number, and both are
    carried in the same pool.
    """
    out: list[str] = []
    current: list[str] = []
    for character in text.lower():
        if character.isalnum():
            current.append(character)
        elif current:
            out.append("".join(current))
            current = []
    if current:
        out.append("".join(current))
    return frozenset(out)


def name_words(text: str) -> list[str]:
    """The words of a name in the order written, keeping repeats.

    The set form loses both, and both carry meaning: "lectin-like 1" and
    "lectin 1" differ only by a repeated word, and "System N" and "N-system"
    differ only by order. Each pair names a different protein.
    """
    out: list[str] = []
    current: list[str] = []
    for character in text.lower():
        if character.isalnum():
            current.append(character)
        elif current:
            out.append("".join(current))
            current = []
    if current:
        out.append("".join(current))
    return out


def ordered_surplus(reference: list[str], candidate: list[str]) -> list[str] | None:
    """The words left over when the reference appears in order, else None.

    The reference name must be recoverable from the candidate by deleting
    words, never by reordering them.
    """
    position = 0
    surplus: list[str] = []
    for word in candidate:
        if position < len(reference) and word == reference[position]:
            position += 1
        else:
            surplus.append(word)
    return surplus if position == len(reference) else None


def usable_names(names: set[str]) -> list[tuple[str, frozenset[str]]]:
    """The reference names admissible for tolerant matching, with their words.

    Two exclusions, both applied here rather than at the point of use so that
    every caller of the tolerant path gets them: names too short or carrying no
    letter, and names whose meaning depends on the order of several numbers.
    """
    usable = []
    for name in names:
        stripped = name.strip()
        if len(stripped) < MINIMUM_NAME_LENGTH:
            continue
        if not any(c.isalpha() for c in stripped):
            continue
        tokens = name_tokens(stripped)
        if not tokens:
            continue
        if sum(1 for t in tokens if t.isdigit()) > MAXIMUM_NUMERIC_WORDS:
            continue
        usable.append((stripped, tokens))
    return usable


def fusion_fragments(element: str) -> list[str]:
    """One recorded antigen split at the comma the source joins partners with.

    An expression construct is recorded as a single antigen naming both the
    fusion partner and the protein -- "Ubiquitin-like protein SMT3,Cadherin-1".
    The protein is present under its exact reference name; it is only the
    partner beside it that hides it. Separating them first means the fusion
    case is settled by exact equality rather than by tolerance.
    """
    parts = [part.strip() for part in element.split(",") if part.strip()]
    return parts if len(parts) > 1 else []


def match_element(element: str, lowered: set[str],
                  usable: list[tuple[str, frozenset[str]]]) -> str | None:
    """How this recorded antigen names the target, or None if it does not.

    Returns MATCHED for exact equality, MATCHED_TOLERANT where the target's
    name is carried inside a longer label, and None where the element names
    something else.
    """
    candidates = [element] + fusion_fragments(element)

    for candidate in candidates:
        if candidate.strip().lower() in lowered:
            return MATCHED

    for candidate in candidates:
        words = name_words(candidate)
        if not words:
            continue
        for name, _tokens in usable:
            surplus = ordered_surplus(name_words(name), words)
            if surplus is None:
                continue
            # Everything the longer label adds must be a qualifier or a short
            # designator. One unrecognised word is enough to decline: the
            # reference name set carries bracket-truncated stems such as
            # "Amine oxidase", derived from "Amine oxidase [copper-containing]
            # 2", and such a stem sits inside every member of its family.
            if not all(_admissible_surplus(w) for w in surplus):
                continue
            return MATCHED_TOLERANT
    return None


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
    usable = usable_names(names)

    informative = [e for e in elements if not _is_absent(e)]
    how = {e: match_element(e, lowered, usable) for e in informative}
    matched = [e for e in informative if how[e] is not None]
    tolerant = [e for e in informative if how[e] == MATCHED_TOLERANT]
    reagents = [e for e in informative if _is_reagent(e)]
    others = [e for e in informative if e not in matched and e not in reagents]

    payload = {
        "target_match_score": None,
        "wrong_antigen_flag": None,
        "recorded_antigens": elements,
        "matched_antigens": matched,
        "tolerantly_matched_antigens": tolerant,
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
        exact = [e for e in matched if how[e] == MATCHED]
        return {
            **payload,
            "target_match": PASS,
            "match_path": MATCHED if exact else MATCHED_TOLERANT,
            "target_match_score": 1.0,
            "wrong_antigen_flag": False,
            "reason": (
                f"The entry's own annotation names this target: "
                f"{', '.join(matched)}"
                + ("" if exact else
                   ", carried inside a longer label rather than recorded under "
                   "the reference name exactly")
                + ". That is a claim about what the "
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
                "tolerantly_matched_antigens": [],
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
                "tolerantly_matched_antigens": [],
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
