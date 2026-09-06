"""B1, B4, B5 and B7 — what is derivable about a target and a binder.

Each module here builds what connected data supports and names what it does
not, rather than emitting a weaker quantity under a stronger name.

**B1** normalises the target from UniProt: identity, synonyms, membrane
attachment, and the mature chain boundaries that carry cleavage and shedding.
It does not carry residues, because the proteome cache holds annotation and no
sequence column.

**B4** separates the original sequence from a cleaned one and never overwrites
provenance. Removal happens only where the source's own metadata justifies it.

**B5** reports what the source states about an antibody -- chain identity,
length, cysteine parity, sequons -- and marks CDR boundaries and germline
similarity UNKNOWN, because no numbering scheme and no germline reference are
connected. A regex approximation reported as "CDRs" would be a weaker
measurement under a stronger name.

**B7** extends the existing developability flags with the motif classes the
benchmark names. It emits no total: Stage 10's standing decision that liability
flags are not summed into a score is not overturned here.
"""

from __future__ import annotations

import math
import re
from collections import Counter

from car_pipeline.stages import stage10

# ---------------------------------------------------------------------------
# B1 -- target normalisation
# ---------------------------------------------------------------------------

def normalise_target(record) -> dict:
    """What UniProt establishes about the target, and what it does not."""
    if record is None:
        return {"status": "UNKNOWN",
                "reason": "no UniProt record for this target"}

    chains = [
        {"start": c.start, "end": c.end, "note": getattr(c, "note", None),
         "chain_id": getattr(c, "chain_id", None),
         "length": (c.end - c.start + 1) if c.start and c.end else None}
        for c in (getattr(record, "chains", None) or [])
    ]
    full = max(chains, key=lambda c: c["length"] or 0) if chains else None
    processed = [c for c in chains if full and c is not full]

    return {
        "status": "NORMALISED",
        "accession": record.accession,
        "gene": record.gene,
        "protein_name": record.protein_name,
        "membrane_class": getattr(record, "membrane_class", None),
        "gpi_anchored": getattr(record, "gpi_anchored", None),
        "transmembrane_count": getattr(record, "transmem_count", None),
        "mature_chains": chains,
        "full_length_chain": full,
        "processed_fragments": processed,
        "shedding_relevant": bool(processed),
        "sequence": None,
        "epitope_regions": None,
        "reasons": [
            (f"{len(processed)} processed fragment(s) are annotated beside the "
             f"full-length chain, so this target is cleaved and a fragment can "
             f"be released. Which fragment a binder engages therefore matters, "
             f"and a binder against a released fragment meets soluble antigen."
             if processed else
             "No processed fragment is annotated, so no cleavage product is "
             "established for this target."),
            "sequence is null: the connected proteome cache holds annotation "
            "with no sequence column, so residue-level work is not available "
            "from it.",
            "epitope_regions is null: no epitope is established without "
            "coordinates or curated epitope evidence, and neither is "
            "connected.",
        ],
    }


# ---------------------------------------------------------------------------
# B4 -- sequence sanitation
# ---------------------------------------------------------------------------

ARTIFACTS = {
    "poly-histidine tag": re.compile(r"H{6,}"),
    "glycine-serine linker": re.compile(r"(?:GGGGS){2,}|(?:GGGS){3,}"),
    "FLAG tag": re.compile(r"DYKDDDDK"),
    "myc tag": re.compile(r"EQKLISEEDL"),
    "Strep-tag II": re.compile(r"WSHPQFEK"),
    "HA tag": re.compile(r"YPYDVPDYA"),
    "TEV protease site": re.compile(r"ENLYFQ[GS]"),
    "thrombin site": re.compile(r"LVPR[GS]S"),
}

# How close to a terminus a motif must sit to be a construct addition rather
# than part of the domain. A tag hanging off an end is a tag; the same motif in
# the middle of a fold is sequence.
TERMINAL_MARGIN = 12


def find_artifacts(sequence: str) -> list[dict]:
    """Every construct artifact in a sequence, with its span and position."""
    out = []
    n = len(sequence or "")
    for name, pattern in ARTIFACTS.items():
        for match in pattern.finditer(sequence or ""):
            start, end = match.start() + 1, match.end()
            terminal = start <= TERMINAL_MARGIN or end >= n - TERMINAL_MARGIN + 1
            out.append({
                "artifact": name,
                "start": start,
                "end": end,
                "residues": match.group(),
                "terminal": terminal,
                "where": ("N-terminal" if start <= TERMINAL_MARGIN
                          else "C-terminal" if terminal else "internal"),
            })
    return sorted(out, key=lambda a: a["start"])


def sanitise(sequence: str | None, domain_declared: bool = False,
             note: str = "") -> dict:
    """Original and cleaned sequences, with every change justified.

    `domain_declared` is the justification: the source states this sequence is
    an antibody variable domain, so a tag or linker hanging off a terminus is a
    construct addition rather than part of it. Without that statement a motif is
    flagged and left in place, because a histidine run in the middle of a
    sequence of unknown extent may be sequence.
    """
    if not sequence:
        return {"sequence_original": None, "sequence_clean": None,
                "modifications": [], "flagged_not_removed": [],
                "reason": "no sequence to sanitise"}

    found = find_artifacts(sequence)
    removable = [a for a in found if a["terminal"] and domain_declared]
    flagged = [a for a in found if a not in removable]

    clean = sequence
    modifications = []
    # Trim from the C-terminal end backwards so earlier spans keep their offsets.
    for artifact in sorted(removable, key=lambda a: a["start"], reverse=True):
        clean = clean[:artifact["start"] - 1] + clean[artifact["end"]:]
        modifications.append({
            **artifact,
            "action": "removed",
            "justification": (
                f"the source declares this sequence an antibody variable "
                f"domain{f' ({note})' if note else ''}, and this "
                f"{artifact['artifact']} sits at its {artifact['where']} end, "
                f"so it is a construct addition rather than part of the domain"),
        })

    return {
        "sequence_original": sequence,
        "sequence_clean": clean if modifications else None,
        "modifications": list(reversed(modifications)),
        "flagged_not_removed": [
            {**a, "action": "flagged",
             "justification": (
                 "no source metadata justifies removing this. An internal "
                 "motif may be sequence rather than a construct addition, and "
                 "removing it on suspicion would destroy provenance."
                 if not a["terminal"] else
                 "the source does not declare what this sequence is, so a "
                 "terminal motif cannot be attributed to a construct.")}
            for a in flagged
        ],
        "reason": (
            f"{len(modifications)} artifact(s) removed and {len(flagged)} "
            f"flagged. The original is carried unchanged in every case."
            if found else
            "No construct artifact detected. The original is carried and no "
            "cleaned sequence is emitted, because there is nothing to clean."),
    }


# ---------------------------------------------------------------------------
# B5 -- antibody annotation
# ---------------------------------------------------------------------------

def annotate(heavy: str | None, light: str | None, fmt: str | None = None) -> dict:
    """What the source establishes about an antibody, and what it does not."""
    def chain(seq, label):
        if not seq:
            return {"present": False, "length": None, "cysteines": None,
                    "cysteine_parity": None, "glycosylation_sequons": None}
        cys = seq.count("C")
        return {
            "present": True,
            "length": len(seq),
            "cysteines": cys,
            "cysteine_parity": "even" if cys % 2 == 0 else "odd",
            "glycosylation_sequons": stage10.glycosylation_sequons(seq),
        }

    return {
        "format": fmt or None,
        "heavy_chain": chain(heavy, "heavy"),
        "light_chain": chain(light, "light"),
        "vh_vl_identified": bool(heavy) and bool(light),
        "cdr_boundaries": None,
        "germline_similarity": None,
        "humanness": None,
        "reasons": [
            "Heavy and light chains are identified because the source records "
            "them separately, not because they were inferred.",
            "cdr_boundaries is null. Assigning complementarity-determining "
            "regions requires a numbering scheme, and no numbering tool is "
            "connected. A motif-anchored guess reported as CDRs would be a "
            "weaker measurement under a stronger name.",
            "germline_similarity and humanness are null. No germline reference "
            "is connected. The naming-convention signal this platform holds "
            "elsewhere is a convention, not a sequence measurement, and it is "
            "not reported here as one.",
        ],
    }


# ---------------------------------------------------------------------------
# B7 -- developability, extended with the motif classes the benchmark names
# ---------------------------------------------------------------------------

MOTIF_LIABILITIES = {
    "deamidation": re.compile(r"N[GSNAH]"),
    "isomerisation": re.compile(r"D[GSDTH]"),
    "oxidation": re.compile(r"[MW]"),
    "hydrolysis": re.compile(r"DP"),
}

LOW_COMPLEXITY_WINDOW = 12
LOW_COMPLEXITY_BITS = 2.6


def _entropy(window: str) -> float:
    """Shannon entropy of one window, in bits."""
    counts = Counter(window)
    n = len(window)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def low_complexity(sequence: str) -> list[int]:
    """Start positions of windows below the entropy floor, one-based."""
    out = []
    for i in range(len(sequence) - LOW_COMPLEXITY_WINDOW + 1):
        if _entropy(sequence[i:i + LOW_COMPLEXITY_WINDOW]) < LOW_COMPLEXITY_BITS:
            out.append(i + 1)
    merged: list[int] = []
    for start in out:
        if not merged or start > merged[-1] + LOW_COMPLEXITY_WINDOW - 1:
            merged.append(start)
    return merged


def developability(sequence: str | None) -> dict:
    """The full liability scorecard. Flags are listed and never summed."""
    if not sequence:
        return {"status": "NOTHING_TO_SCORE",
                "reason": "no sequence to assess"}

    motifs = {
        name: [m.start() + 1 for m in pattern.finditer(sequence)]
        for name, pattern in MOTIF_LIABILITIES.items()
    }
    cys = sequence.count("C")
    return {
        "status": "SCORED",
        "length": len(sequence),
        "isoelectric_point": stage10.isoelectric_point(sequence),
        "net_charge_at_formulation_ph": round(
            stage10.net_charge(sequence, stage10.FORMULATION_PH), 4),
        "gravy": stage10.gravy(sequence),
        "aggregation_prone_regions": stage10.aggregation_prone(sequence),
        "glycosylation_sequons": stage10.glycosylation_sequons(sequence),
        "cysteines": cys,
        "unpaired_cysteine": cys % 2 == 1,
        "deamidation_motifs": motifs["deamidation"],
        "isomerisation_motifs": motifs["isomerisation"],
        "oxidation_sites": motifs["oxidation"],
        "hydrolysis_motifs": motifs["hydrolysis"],
        "low_complexity_windows": low_complexity(sequence),
        "developability_score": None,
        "reasons": [
            "developability_score is null by standing decision, not for want "
            "of inputs. Liability flags are counted and listed and never "
            "summed, because a flag that fires on every binder in a pool "
            "carries no information and a total would hide that.",
            "Motif counts are positions where a liability is possible, not "
            "measurements that it occurs. Whether a site is exposed enough to "
            "deamidate or oxidise depends on the folded structure, and no "
            "coordinates are connected.",
        ],
    }
