"""Stage 9 — genomic and construct safety, read off Stage 6's DNA map."""

from __future__ import annotations

import re
from dataclasses import dataclass

CODON_INVARIANT = "CODON_INVARIANT"
MAP_SPECIFIC = "MAP_SPECIFIC"


PEPTIDE_REPEAT_MIN = 8
NT_REPEAT_MIN = 24
TRACT_MIN = 8
ORF_MIN = 30
HOMOPOLYMER_MIN = 8


LISTED = 6


THRESHOLDS = {
    "peptide_repeat_min_residues": PEPTIDE_REPEAT_MIN,
    "nucleotide_repeat_min_bp": NT_REPEAT_MIN,
    "pyrimidine_tract_min": TRACT_MIN,
    "orf_min_codons": ORF_MIN,
    "homopolymer_min_bp": HOMOPOLYMER_MIN,
}


COMPLEMENT = str.maketrans("ACGT", "TGCA")

DONOR = re.compile(r"GT[AG]AGT")
ACCEPTOR = re.compile(rf"[CT]{{{TRACT_MIN},}}[ACGT]?AG")
HOMOPOLYMER = re.compile(rf"(A{{{HOMOPOLYMER_MIN},}}|C{{{HOMOPOLYMER_MIN},}}"
                         rf"|G{{{HOMOPOLYMER_MIN},}}|T{{{HOMOPOLYMER_MIN},}})")

STOPS = ("TAA", "TAG", "TGA")


@dataclass(frozen=True)
class Finding:
    """One construct-safety observation, and whether the encoding could move it."""

    kind: str
    basis: str
    at: str
    detail: str

    def as_payload(self) -> dict:
        """The finding, with the basis stated rather than implied."""
        return {"kind": self.kind, "basis": self.basis, "at": self.at,
                "detail": self.detail}


def reverse_complement(dna: str) -> str:
    """The other strand, read 5' to 3'."""
    return dna.translate(COMPLEMENT)[::-1]


def repeated_parts(segments) -> list[Finding]:
    """Parts appearing twice in one construct, read from the domain map."""
    seen: dict[tuple, list] = {}
    for segment in segments:
        if not segment.accession:
            continue
        key = (segment.accession, segment.start_residue, segment.end_residue)
        seen.setdefault(key, []).append(segment)
    out = []
    for (accession, start, end), members in sorted(
            seen.items(), key=lambda kv: -len(kv[1])):
        if len(members) < 2:
            continue
        span = f"{start}-{end}" if start else "whole part"
        out.append(Finding(
            "repeated_part", CODON_INVARIANT,
            ", ".join(f"{m.bp_start}-{m.bp_end}" for m in members),
            f"{members[0].name} ({accession} {span}) appears {len(members)} "
            f"times, {members[0].residues} residues each. Identical stretches "
            f"of this length are the dominant recombination hazard in a vector, "
            f"and no codon choice removes them"))
    return out


def repeated_peptides(protein: str) -> list[Finding]:
    """Peptides recurring in the construct, which every encoding preserves."""
    positions: dict[str, list[int]] = {}
    for i in range(len(protein) - PEPTIDE_REPEAT_MIN + 1):
        positions.setdefault(protein[i:i + PEPTIDE_REPEAT_MIN], []).append(i)
    out = []
    for peptide, at in positions.items():
        if len(at) < 2:
            continue
        spaced = _spaced(at, PEPTIDE_REPEAT_MIN)
        if len(spaced) < 2:
            continue
        out.append(Finding(
            "repeated_peptide", CODON_INVARIANT,
            ", ".join(str(p + 1) for p in spaced),
            f"{peptide} occurs {len(spaced)} times without overlapping"))
    return out


def internal_methionines(protein: str) -> list[Finding]:
    """Methionines after the initiator, each a candidate internal start."""
    at = [i + 1 for i, residue in enumerate(protein) if residue == "M" and i]
    if not at:
        return []
    return [Finding(
        "internal_methionine", CODON_INVARIANT,
        ", ".join(str(p) for p in at[:LISTED]),
        f"{len(at)} methionine(s) after the initiator. Each is a possible "
        f"internal initiation site in the coding frame under any encoding")]


def _spaced(at: list[int], size: int) -> list[int]:
    """Occurrences that do not overlap the one kept before them."""
    kept = [at[0]]
    for position in at[1:]:
        if position - kept[-1] >= size:
            kept.append(position)
    return kept


def _kmers(dna: str, size: int) -> dict[str, list[int]]:
    """Every substring of this size and where it starts."""
    out: dict[str, list[int]] = {}
    for i in range(len(dna) - size + 1):
        out.setdefault(dna[i:i + size], []).append(i)
    return out


def direct_repeats(dna: str) -> list[Finding]:
    """Nucleotide stretches occurring more than once on this strand."""
    out = []
    for kmer, at in _kmers(dna, NT_REPEAT_MIN).items():
        if len(at) < 2:
            continue
        spaced = _spaced(at, NT_REPEAT_MIN)
        if len(spaced) < 2:
            continue
        out.append(Finding(
            "direct_repeat", MAP_SPECIFIC,
            ", ".join(str(p + 1) for p in spaced),
            f"{NT_REPEAT_MIN} bp repeated {len(spaced)} times without "
            f"overlapping"))
    return out


def inverted_repeats(dna: str) -> list[Finding]:
    """Stretches whose reverse complement also occurs."""
    table = _kmers(dna, NT_REPEAT_MIN)
    out, seen = [], set()
    for kmer, at in table.items():
        mirror = reverse_complement(kmer)
        if mirror not in table or kmer in seen:
            continue
        elsewhere = [q for q in table[mirror] if q != at[0]]
        if not elsewhere:
            continue
        seen.add(kmer)
        seen.add(mirror)
        out.append(Finding(
            "inverted_repeat", MAP_SPECIFIC,
            f"{at[0] + 1} and {elsewhere[0] + 1}",
            f"{NT_REPEAT_MIN} bp and its reverse complement both occur"))
    return out


def splice_sites(dna: str) -> list[Finding]:
    """Donor and acceptor motifs, which the codon table places arbitrarily."""
    out = []
    donors = [m.start() + 1 for m in DONOR.finditer(dna)]
    if donors:
        out.append(Finding(
            "splice_donor", MAP_SPECIFIC,
            ", ".join(str(p) for p in donors[:LISTED]),
            f"{len(donors)} match(es) to the GT[AG]AGT donor consensus"))
    acceptors = [m.end() - 1 for m in ACCEPTOR.finditer(dna)]
    if acceptors:
        out.append(Finding(
            "splice_acceptor", MAP_SPECIFIC,
            ", ".join(str(p) for p in acceptors[:LISTED]),
            f"{len(acceptors)} AG preceded by a pyrimidine tract of "
            f"{TRACT_MIN} or more"))
    return out


def _orfs(dna: str, frames: tuple[int, ...], kind: str) -> list[Finding]:
    """Reading frames other than the annotated one, ATG to stop."""
    where = (", numbered along the reverse strand"
             if kind == "reverse_orf" else "")
    out = []
    for frame in frames:
        i = frame
        while i < len(dna) - 2:
            if dna[i:i + 3] != "ATG":
                i += 3
                continue
            j, stopped = i + 3, False
            while j < len(dna) - 2:
                if dna[j:j + 3] in STOPS:
                    stopped = True
                    break
                j += 3
            if not stopped:
                break
            codons = (j - i) // 3
            if codons >= ORF_MIN:
                out.append(Finding(
                    kind, MAP_SPECIFIC, f"{i + 1}-{j + 3}",
                    f"{codons} codons in frame {frame + 1}{where}"))
            i = j + 3
    return out


def alternate_orfs(dna: str) -> list[Finding]:
    """Open reading frames in the two frames the construct does not use."""
    return _orfs(dna, (1, 2), "alternate_orf")


def reverse_orfs(dna: str) -> list[Finding]:
    """Open reading frames on the other strand."""
    return _orfs(reverse_complement(dna), (0, 1, 2), "reverse_orf")


def homopolymers(dna: str) -> list[Finding]:
    """Single-base runs long enough to trouble synthesis."""
    out = []
    for match in HOMOPOLYMER.finditer(dna):
        out.append(Finding(
            "homopolymer", MAP_SPECIFIC, f"{match.start() + 1}-{match.end()}",
            f"{len(match.group(0))} x {match.group(0)[0]}"))
    return out


def findings(protein: str, dna: str, segments) -> list[Finding]:
    """Every construct-safety finding, invariant ones first."""
    return (
        repeated_parts(segments)
        + repeated_peptides(protein)
        + internal_methionines(protein)
        + direct_repeats(dna)
        + inverted_repeats(dna)
        + splice_sites(dna)
        + alternate_orfs(dna)
        + reverse_orfs(dna)
        + homopolymers(dna)
    )


def analyse(protein: str, dna: str, segments) -> dict:
    """The arm's report for one construct. It gates nothing and changes nothing."""
    found = findings(protein, dna, segments)
    invariant = [f for f in found if f.basis == CODON_INVARIANT]
    specific = [f for f in found if f.basis == MAP_SPECIFIC]
    counts: dict[str, int] = {}
    for f in found:
        counts[f.kind] = counts.get(f.kind, 0) + 1
    return {
        "analysed": True,
        "thresholds": dict(THRESHOLDS),
        "counts": counts,
        "codon_invariant": [f.as_payload() for f in invariant[:LISTED * 2]],
        "codon_invariant_total": len(invariant),
        "map_specific": [f.as_payload() for f in specific[:LISTED * 2]],
        "map_specific_total": len(specific),
        "reasons": [
            "This is read from Stage 6's DNA map, which is a reverse "
            "translation under one fixed codon per residue and not an ordering "
            "sequence. Findings marked MAP_SPECIFIC are properties of that "
            "encoding: a different codon table moves them, and the sequence "
            "anyone would order has different ones.",
            "Findings marked CODON_INVARIANT survive any encoding, because they "
            "are properties of the protein or of the domain layout. Those are "
            "the ones a reader can act on.",
            "Nothing here gates. No construct is blocked and no verdict "
            "changes; setting a tolerance on these counts would need outcome "
            "data this platform does not have.",
        ],
    }
