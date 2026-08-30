"""Stage 6 — construct assembly."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from car_pipeline.stages.stage4 import ADAPTOR
from car_pipeline.data.domains import (
    SYNTHETIC_PARTS,
    Part,
    anti_tag_binder,
    build_parts,
)

BUILDABLE = "BUILDABLE"
BUDGET_EXCEEDED = "BUDGET_EXCEEDED"
NO_CONSTRUCT = "NO_CONSTRUCT"


BUDGET_BP = 3500
STOP = "TAA"


CODON = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
    "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG",
    "Y": "TAC", "V": "GTG",
}
REVERSE_CODON = {v: k for k, v in CODON.items()}


BUILDABLE_OUTCOMES = ("SINGLE", "DUAL", "ADAPTOR")


def assemblable(sequence: str) -> set[str]:
    """Residues in this sequence that the codon table cannot encode."""
    return set(sequence) - set(CODON)


def translate(dna: str) -> str:
    """Back to residues, stopping at the stop codon."""
    out = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        if codon == STOP:
            break
        out.append(REVERSE_CODON.get(codon, "?"))
    return "".join(out)


@dataclass
class Segment:
    name: str
    provenance: str
    accession: str
    feature: str
    start_residue: int | None
    end_residue: int | None
    aa_start: int
    aa_end: int
    bp_start: int
    bp_end: int

    @property
    def residues(self) -> int:
        """The segment's length in residues."""
        return self.aa_end - self.aa_start


@dataclass
class Construct:
    gene: str
    accession: str
    pool_index: int
    outcome: str
    partner: str | None
    verdict: str
    architecture: str = ""
    binder_name: str = ""
    partner_binder_name: str = ""
    amino_acid_sequence: str = ""
    dna: str = ""
    segments: list[Segment] = field(default_factory=list)
    reason: str = ""

    binder_supplied: bool = True

    declared_bp: int | None = None

    @property
    def total_bp(self) -> int:
        """The construct's length, falling back to the declared size."""
        if not self.dna and self.declared_bp is not None:
            return self.declared_bp
        return len(self.dna)

    @property
    def headroom_bp(self) -> int:
        """What is left of the payload budget."""
        return BUDGET_BP - self.total_bp

    @property
    def has_switch(self) -> bool:
        """Whether a safety switch is among the segments."""
        return any("caspase" in s.name.lower() for s in self.segments)


def _scfv(vh: str, vl: str, linker: Part) -> list[tuple[str, Part]]:
    """The variable regions joined by the linker, in order."""
    return [
        ("VH", Part("binder VH", "stage5", vh)),
        ("linker", linker),
        ("VL", Part("binder VL", "stage5", vl)),
    ]


def _assemble(pieces: list[tuple[str, Part]]) -> tuple[str, str, list[Segment]]:
    """Concatenate, reverse-translate, and record every boundary."""
    aa_parts, segments = [], []
    aa_pos = 0

    complete = all(part.supplied for _l, part in pieces)
    for _label, part in pieces:
        seq = part.sequence
        segments.append(
            Segment(
                name=part.name,
                provenance=part.provenance,
                accession=part.accession,
                feature=part.feature,
                start_residue=part.start,
                end_residue=part.end,
                aa_start=aa_pos,
                aa_end=aa_pos + part.residues,
                bp_start=aa_pos * 3,
                bp_end=(aa_pos + part.residues) * 3,
            )
        )
        aa_parts.append(seq)
        aa_pos += part.residues
    if not complete:
        return "", "", segments
    protein = "".join(aa_parts)
    unknown = assemblable(protein)
    if unknown:
        raise ValueError(
            "cannot encode residues " + ",".join(sorted(unknown))
            + "; a construct whose DNA does not encode its own protein must not "
              "be produced"
        )
    dna = "".join(CODON[residue] for residue in protein) + STOP
    return protein, dna, segments


def build(
    decisions: list[dict],
    binders: dict[str, object],
    parts: dict[str, Part] | None = None,
) -> list[Construct]:
    """One record per pool member, in the order Stage 4 emitted them."""
    parts = parts or build_parts()
    linker = SYNTHETIC_PARTS["linker"]
    skip = SYNTHETIC_PARTS["skip"]
    switch_linker = SYNTHETIC_PARTS["switch_linker"]
    adaptor_binder = anti_tag_binder()

    def best_binder(gene: str):
        """The shortest sequence-route binder: the smallest that fits is the"""
        record = binders.get(gene)
        if record is None:
            return None
        usable = [
            c for c in record.sequence
            if c.heavy_sequence and c.light_sequence
            and not assemblable(c.heavy_sequence + c.light_sequence)
        ]
        if not usable:
            return None
        return min(usable, key=lambda c: len(c.heavy_sequence) + len(c.light_sequence))

    switch = [
        ("skip", skip),
        ("FKBP12", parts["switch_fkbp"]),
        ("switch linker", switch_linker),
        ("caspase", parts["switch_caspase"]),
    ]

    out: list[Construct] = []
    for row in decisions:
        gene = row["gene"]
        target = best_binder(gene)
        partner_gene = row.get("partner")
        partner = best_binder(partner_gene) if partner_gene else None

        if row["outcome"] == ADAPTOR:
            pieces = (
                [("leader", parts["leader"]),
                 ("anti-tag", adaptor_binder),
                 ("hinge", parts["hinge"]), ("TM", parts["transmembrane"]),
                 ("4-1BB", parts["costimulatory"]),
                 ("CD3zeta", parts["activation"])]
                + switch
            )
            protein, dna, segments = _assemble(pieces)
            total = sum(seg.bp_end - seg.bp_start for seg in segments) + len(STOP)
            out.append(Construct(
                gene=gene, accession=row["accession"], pool_index=row["pool_index"],
                outcome=row["outcome"], partner=None,
                verdict=BUILDABLE if total <= BUDGET_BP else BUDGET_EXCEEDED,
                architecture="adaptor, anti-tag receptor, antigen on the adaptor",
                binder_name=adaptor_binder.name,
                amino_acid_sequence=protein, dna=dna, segments=segments,
                binder_supplied=adaptor_binder.supplied, declared_bp=total,
                reason=(
                    f"the anti-tag binder is retrieved from {adaptor_binder.accession}, "
                    f"{adaptor_binder.feature}, as deposited and unedited"
                    if adaptor_binder.supplied else
                    "the anti-tag binder declares a size but no sequence; "
                    "no anti-tag antibody exists in the cached structural "
                    "set and none was invented"
                ),
            ))
            continue

        if target is None:
            out.append(Construct(
                gene=gene, accession=row["accession"], pool_index=row["pool_index"],
                outcome=row["outcome"], partner=partner_gene,
                verdict=NO_CONSTRUCT,
                reason="no binder with a usable variable-region sequence "
                       "(the structure route carries no sequence to assemble)",
            ))
            continue

        if row["outcome"] not in BUILDABLE_OUTCOMES:
            out.append(Construct(
                gene=gene, accession=row["accession"], pool_index=row["pool_index"],
                outcome=row["outcome"], partner=partner_gene,
                verdict=NO_CONSTRUCT,
                binder_name=target.name,
                reason=f"a binder exists, but Stage 4 returned {row['outcome']}; "
                       "this stage does not build for a target the pairing stage "
                       "did not recommend",
            ))
            continue

        dual = row["outcome"] == "DUAL"
        if dual and partner is None:
            out.append(Construct(
                gene=gene, accession=row["accession"], pool_index=row["pool_index"],
                outcome=row["outcome"], partner=partner_gene,
                verdict=NO_CONSTRUCT,
                binder_name=target.name,
                reason=f"dual design, but the partner {partner_gene} has no binder",
            ))
            continue

        if dual:
            pieces = (
                [("leader", parts["leader"])]
                + _scfv(target.heavy_sequence, target.light_sequence, linker)
                + [("hinge", parts["hinge"]), ("TM", parts["transmembrane"]),
                   ("CD3zeta", parts["activation"]), ("skip", skip),
                   ("leader", parts["leader"])]
                + _scfv(partner.heavy_sequence, partner.light_sequence, linker)
                + [("hinge", parts["hinge"]), ("TM", parts["transmembrane"]),
                   ("4-1BB", parts["costimulatory"])]
                + switch
            )
            architecture = "dual, split signal, two receptors"
        else:
            pieces = (
                [("leader", parts["leader"])]
                + _scfv(target.heavy_sequence, target.light_sequence, linker)
                + [("hinge", parts["hinge"]), ("TM", parts["transmembrane"]),
                   ("4-1BB", parts["costimulatory"]),
                   ("CD3zeta", parts["activation"])]
                + switch
            )
            architecture = "single, second generation"

        protein, dna, segments = _assemble(pieces)
        construct = Construct(
            gene=gene, accession=row["accession"], pool_index=row["pool_index"],
            outcome=row["outcome"], partner=partner_gene,
            verdict=BUILDABLE if len(dna) <= BUDGET_BP else BUDGET_EXCEEDED,
            architecture=architecture,
            binder_name=target.name,
            partner_binder_name=partner.name if (dual and partner) else "",
            amino_acid_sequence=protein, dna=dna, segments=segments,
        )
        if construct.verdict == BUDGET_EXCEEDED:
            construct.reason = (
                f"{construct.total_bp} bp against a {BUDGET_BP} bp budget, "
                f"over by {-construct.headroom_bp}"
            )
        out.append(construct)
    return out


def configuration_hash(stage5_hash: str, genes: list[str]) -> str:
    """Covers Stage 5's configuration, not Stage 4's."""
    payload = {
        "stage5": stage5_hash,
        "genes": genes,
        "budget_bp": BUDGET_BP,
        "codon": CODON,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
