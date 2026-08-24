"""Stage 10 — developability.

Implements `specs/stage10-developability.md`. Arithmetic over sequences that
already exist; no source is added.

Scores the binder variable regions from Stage 5 rather than constructs, because
Stage 6 produces zero buildable constructs and scoring those would score nothing.

**Nothing here predicts manufacturing failure.** These are sequence liabilities —
properties that correlate with difficulty — reported as counts and values. They
are never summed: one number would let a strong liability be averaged away by
four weak absences.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

NOTHING_TO_SCORE = "NOTHING_TO_SCORE"

#: Side-chain pKa values, plus the termini. One fixed table so the isoelectric
#: point is reproducible; the exact values differ between published sets and the
#: choice is stated rather than left implicit.
PKA_SIDE = {"C": 8.5, "D": 3.9, "E": 4.1, "H": 6.0, "K": 10.5, "R": 12.5, "Y": 10.1}
PKA_N_TERM = 9.7
PKA_C_TERM = 2.3
POSITIVE = {"K", "R", "H"}

#: Kyte–Doolittle hydropathy.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: Fixed before the run. See the specification.
FORMULATION_PH = 7.4
PI_WINDOW = 1.0
CHARGE_FLOOR = 1.0
APR_WINDOW = 7
APR_THRESHOLD = 1.0


def net_charge(sequence: str, ph: float) -> float:
    """Net charge at a given pH, from the fixed pKa table."""
    charge = 1.0 / (1.0 + 10 ** (ph - PKA_N_TERM))
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERM - ph))
    for residue in sequence:
        pka = PKA_SIDE.get(residue)
        if pka is None:
            continue
        if residue in POSITIVE:
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        else:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(sequence: str) -> float:
    """pH of zero net charge, by bisection.

    Bounded to 1..14 by construction. A sequence with no ionisable side chain
    still has termini, so the point always exists inside those bounds.
    """
    low, high = 1.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        if net_charge(sequence, mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 3)


def glycosylation_sequons(sequence: str) -> list[int]:
    """Positions of `N-X-S/T` where X is not proline, one-based.

    The proline exclusion is the part of this rule most easily dropped, and
    dropping it silently inflates every count. It is pinned by a criterion.
    """
    out = []
    for i in range(len(sequence) - 2):
        n, x, t = sequence[i], sequence[i + 1], sequence[i + 2]
        if n == "N" and x != "P" and t in ("S", "T"):
            out.append(i + 1)
    return out


def aggregation_prone(sequence: str) -> list[int]:
    """Start positions of hydrophobic windows, one-based, merged if overlapping."""
    hits: list[int] = []
    for i in range(len(sequence) - APR_WINDOW + 1):
        window = sequence[i:i + APR_WINDOW]
        values = [HYDROPATHY.get(r) for r in window]
        if any(v is None for v in values):
            continue
        if sum(values) / APR_WINDOW >= APR_THRESHOLD:
            if not hits or i + 1 > hits[-1] + APR_WINDOW - 1:
                hits.append(i + 1)
    return hits


def gravy(sequence: str) -> float | None:
    values = [HYDROPATHY.get(r) for r in sequence]
    usable = [v for v in values if v is not None]
    if not usable:
        return None
    return round(sum(usable) / len(usable), 4)


@dataclass
class Developability:
    gene: str
    binder: str
    residues: int
    isoelectric_point: float
    net_charge: float
    cysteines: int
    cysteine_parity: str
    sequons: list[int] = field(default_factory=list)
    apr_starts: list[int] = field(default_factory=list)
    gravy: float | None = None
    flags: list[tuple[str, str]] = field(default_factory=list)

    @property
    def flag_count(self) -> int:
        return len(self.flags)

    @property
    def kinds(self) -> list[str]:
        return [kind for kind, _ in self.flags]


def score(sequence: str, gene: str, binder: str) -> Developability:
    pi = isoelectric_point(sequence)
    charge = round(net_charge(sequence, FORMULATION_PH), 4)
    cys = sequence.count("C")
    sequons = glycosylation_sequons(sequence)
    aprs = aggregation_prone(sequence)

    # (kind, detail). The kind is a stable key so the report can group without
    # parsing prose back out of a sentence.
    flags = []
    if abs(pi - FORMULATION_PH) <= PI_WINDOW:
        flags.append(("pI near formulation pH",
                      f"pI {pi} within {PI_WINDOW} of {FORMULATION_PH}"))
    if abs(charge) < CHARGE_FLOOR:
        flags.append(("low net charge",
                      f"net charge {charge} below {CHARGE_FLOOR}"))
    if cys % 2 == 1:
        flags.append(("odd cysteine count",
                      f"{cys} cysteines, so at least one is unpaired"))
    if sequons:
        flags.append(("N-glycosylation sequon",
                      f"{len(sequons)} sequon(s) at {sequons[:4]}"))
    if aprs:
        flags.append(("aggregation-prone region",
                      f"{len(aprs)} region(s) starting {aprs[:4]}"))

    return Developability(
        gene=gene,
        binder=binder,
        residues=len(sequence),
        isoelectric_point=pi,
        net_charge=charge,
        cysteines=cys,
        # Never "unpaired: 0". Pairing cannot be read from sequence, and an even
        # count is only the absence of a guarantee, not evidence of pairing.
        cysteine_parity="odd" if cys % 2 else "even",
        sequons=sequons,
        apr_starts=aprs,
        gravy=gravy(sequence),
        flags=flags,
    )


def assess(binders: dict) -> tuple[list[Developability], str]:
    """Score every binder carrying a sequence, in pool order.

    Returns the rows and a status. An empty result is `NOTHING_TO_SCORE`, never
    an empty table: a table with no rows reads as "nothing had liabilities",
    which is the opposite of "nothing was examined".
    """
    rows: list[Developability] = []
    for gene, record in binders.items():
        for candidate in record.sequence:
            vh, vl = candidate.heavy_sequence, candidate.light_sequence
            if not vh:
                continue
            rows.append(score(vh + vl, gene, candidate.name))
    status = NOTHING_TO_SCORE if not rows else "SCORED"
    return rows, status


def configuration_hash(stage5_hash: str, genes: list[str]) -> str:
    payload = {
        "stage5": stage5_hash,
        "genes": genes,
        "pka": PKA_SIDE,
        "ph": FORMULATION_PH,
        "thresholds": [PI_WINDOW, CHARGE_FLOOR, APR_WINDOW, APR_THRESHOLD],
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
