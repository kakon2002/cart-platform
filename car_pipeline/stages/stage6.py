"""Stage 6 — construct assembly.

Implements `specs/stage6-construct-assembly.md`. Mechanical by design: every
constraint is fixed upstream and nothing here re-derives one.

The DNA is a **map, not an ordering sequence**. It is reverse-translated under one
fixed codon per amino acid so that domain boundaries are exact and the round trip
is checkable. It is not codon-optimised and must not be read as though it were.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from car_pipeline.stages.stage4 import ADAPTOR
from car_pipeline.data.domains import (
    SYNTHETIC_PARTS,
    Part,
    build_parts,
)

BUILDABLE = "BUILDABLE"
BUDGET_EXCEEDED = "BUDGET_EXCEEDED"
NO_CONSTRUCT = "NO_CONSTRUCT"

#: Stage 1's payload budget, in bases. Carried, not recomputed.
BUDGET_BP = 3500
STOP = "TAA"

#: One codon per amino acid, fixed so the map is deterministic and the round trip
#: is a real check. Human-frequent choices, but frequency is not the point:
#: reproducibility is.
CODON = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
    "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG",
    "Y": "TAC", "V": "GTG",
}
REVERSE_CODON = {v: k for k, v in CODON.items()}


#: Outcomes this stage will build for. A construct assembled for a target the
#: pairing stage rejected would be a design presented for something upstream says
#: is not designable, and a reader cannot be expected to carry that caveat.
#: ADAPTOR joins these: it is a real routed architecture that
#: assembles, even though its binder sequence is not supplied.
BUILDABLE_OUTCOMES = ("SINGLE", "DUAL", "ADAPTOR")


def assemblable(sequence: str) -> set[str]:
    """Residues in this sequence that the codon table cannot encode.

    Returned rather than silently mapped. One retrieved therapeutic carries
    lowercase residues, and `CODON.get(residue, "NNN")` turned those into an
    ambiguous codon that translated back as a mismatch — a construct whose DNA
    did not encode its own protein, produced without an error.
    """
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
    aa_end: int          # half-open
    bp_start: int
    bp_end: int          # half-open

    @property
    def residues(self) -> int:
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
    #: False when a part in this construct declares a size but no sequence --
    #: the adaptor receptor's anti-tag binder is the only such part. The layout
    #: and the length are real; the residues are not supplied and are not
    #: invented, so `amino_acid_sequence` and `dna` are left empty rather than
    #: filled with something that would read as a design.
    binder_supplied: bool = True
    #: The construct's length when there is no DNA to measure it from.
    declared_bp: int | None = None

    @property
    def total_bp(self) -> int:
        if not self.dna and self.declared_bp is not None:
            return self.declared_bp
        return len(self.dna)

    @property
    def headroom_bp(self) -> int:
        return BUDGET_BP - self.total_bp

    @property
    def has_switch(self) -> bool:
        return any("caspase" in s.name.lower() for s in self.segments)


def _scfv(vh: str, vl: str, linker: Part) -> list[tuple[str, Part]]:
    return [
        ("VH", Part("binder VH", "stage5", vh)),
        ("linker", linker),
        ("VL", Part("binder VL", "stage5", vl)),
    ]


def _assemble(pieces: list[tuple[str, Part]]) -> tuple[str, str, list[Segment]]:
    """Concatenate, reverse-translate, and record every boundary."""
    aa_parts, segments = [], []
    aa_pos = 0
    # A part may declare a length without residues. Coordinates still advance by
    # its declared size so the map and the budget are right; the protein is
    # withheld below rather than padded with filler.
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
        # Length and layout are real; the residues are not supplied. Returning
        # an empty protein rather than a padded one means nothing downstream can
        # mistake this for a sequence that was designed.
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
    adaptor_binder = SYNTHETIC_PARTS["adaptor_binder"]

    def best_binder(gene: str):
        """The shortest sequence-route binder: the smallest that fits is the
        most favourable reading of the budget, and it is labelled as a bound."""
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

        # An adaptor receptor binds the tag, not the antigen, so it needs no
        # target binder and is built before the binder guard below. Its
        # antigen specificity lives in a separately manufactured adaptor that
        # is not in this vector -- which is the whole reason it fits.
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
                binder_supplied=False, declared_bp=total,
                reason="the anti-tag binder declares a size but no sequence; "
                       "no anti-tag antibody exists in the cached structural "
                       "set and none was invented",
            ))
            continue

        if target is None:
            out.append(Construct(
                gene=gene, accession=row["accession"], pool_index=row["pool_index"],
                outcome=row["outcome"], partner=partner_gene,
                verdict=NO_CONSTRUCT,
                # Assembly needs a sequence, so this reads the sequence route
                # only. A target with a structure-route binder and no sequence
                # lands here, and saying "no binder" would contradict Stage 5.
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
            # Split signal: activation on one receptor, costimulation on the
            # other, so neither antigen alone gives a complete signal.
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
            # Only for a dual. Stage 4 fills `partner` on the SINGLE branch too,
            # so labelling a single-arm construct from it names a binder the
            # sequence does not contain.
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
    """Covers Stage 5's configuration, not Stage 4's.

    Passing the upstream-of-upstream hash here would leave a Stage 5 change
    invisible to this stage's identity, which is the whole reason the hash is
    carried.
    """
    payload = {
        "stage5": stage5_hash,
        "genes": genes,
        "budget_bp": BUDGET_BP,
        "codon": CODON,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
