"""B9 — deposited structural evidence, gated on the target-match check.

The benchmark's own rule: use deposited structural evidence only where target
match is verified, otherwise return UNKNOWN. That is the right gate, and it is
the reason this module reads a verdict before it reads a structure. An entry
found by searching on a target's accession contains that target somewhere; the
antibody in it may be annotated against another chain entirely, and reporting
its geometry as evidence about the target would be the error the check exists
to catch, one stage later.

What this returns is **retrieved, not predicted**: the entry, its experimental
method, its resolution and the chains involved. What it does not return is
interface geometry -- which residues contact which -- because that needs
coordinates and no coordinate file is connected anywhere in this platform.

Resolution is carried per entry and never averaged, and no cutoff is imposed.
The connected MSLN entries span 1.52 to 4.31 angstrom, which is wide enough
that a single figure would flatten a real difference; choosing a threshold after
seeing that distribution would be a bound fitted to an observation.
"""

from __future__ import annotations

from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.stages import binder_check

VERIFIED = "VERIFIED_EXPERIMENTAL"
UNKNOWN = "UNKNOWN"

# What the evidence rests on, said rather than implied.
BASIS = "deposited experimental structure, retrieved; not predicted"

_INDEX: dict[str, list] | None = None


def _index() -> dict[str, list]:
    """The deposited-structure index, loaded once."""
    global _INDEX
    if _INDEX is None:
        _INDEX = AntibodySource().structures()
    return _INDEX


def _instance(identifier: str | None):
    """The deposited instance an identifier names, if it can be found.

    Stage 5 forms the identifier as entry:heavy+light. Matching on the entry
    alone is not enough -- an entry carries several antibody instances and they
    differ in antigen, method is shared but chains are not -- so the chains are
    matched too.
    """
    if not identifier or ":" not in identifier:
        return None
    entry, _, chains = identifier.partition(":")
    for candidate in _index().get(entry.strip().lower(), []):
        if f"{candidate.heavy_chain}{candidate.light_chain}" == chains:
            return candidate
    return None


def evidence(row: dict) -> dict:
    """Structural evidence for one retrieved binder, or a named absence."""
    verdict = row.get("target_match")

    if verdict != binder_check.PASS:
        return {
            "structural_status": UNKNOWN,
            "entry": None,
            "method": None,
            "resolution": None,
            "antigen_chain": None,
            "antigen_type": None,
            "heavy_chain": None,
            "light_chain": None,
            "interface_geometry": None,
            "basis": None,
            "reason": (
                "Target match did not pass, so no structural evidence is "
                "reported. A deposited entry found by searching on this "
                "target's accession contains the target, but the antibody in "
                "it is annotated against something else or against nothing "
                "recorded; its geometry is not evidence about this target."
                if verdict == binder_check.FAIL else
                "Target match is UNKNOWN, so no structural evidence is "
                "reported. Absent or inapplicable antigen annotation is not a "
                "reason to treat a structure as verified."),
        }

    instance = _instance(row.get("identifier"))
    if instance is None:
        return {
            "structural_status": UNKNOWN,
            "entry": None,
            "method": row.get("method") or None,
            "resolution": None,
            "antigen_chain": row.get("antigen_chain") or None,
            "antigen_type": None,
            "heavy_chain": None,
            "light_chain": None,
            "interface_geometry": None,
            "basis": None,
            "reason": (
                "Target match passed, but the deposited instance this binder "
                "names could not be located in the connected release, so its "
                "method and resolution cannot be reported."),
        }

    resolution = (instance.resolution or "").strip()
    return {
        "structural_status": VERIFIED,
        "entry": instance.pdb.upper(),
        "method": instance.method or None,
        "resolution": float(resolution) if _numeric(resolution) else None,
        "resolution_recorded": resolution or None,
        "antigen_chain": instance.antigen_chain or None,
        "antigen_type": instance.antigen_type or None,
        "heavy_chain": instance.heavy_chain or None,
        "light_chain": instance.light_chain or None,
        "interface_geometry": None,
        "basis": BASIS,
        "reason": (
            f"A deposited {instance.method or 'experimental'} structure whose "
            f"recorded antigen names this target"
            + (f", solved at {resolution} angstrom" if resolution else "")
            + ". Interface geometry is null: which residues contact which "
              "needs coordinates, and no coordinate file is connected. What is "
              "reported is that a verified complex exists and what kind it is, "
              "never what its interface looks like."),
    }


def _numeric(value: str) -> bool:
    """Whether a recorded resolution is a number rather than a placeholder."""
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def attach(rows: list[dict]) -> list[dict]:
    """Every row with its structural evidence, gated on its own verdict."""
    return [{**row, "structural_evidence": evidence(row)} for row in rows]


NOTES = [
    "Structural evidence is gated on the target-match verdict. A binder whose "
    "recorded antigen names something else carries none, rather than carrying "
    "geometry that would be read as evidence about this target.",
    "What is reported is retrieved, not predicted: the entry, its method, its "
    "resolution and its chains. Interface geometry is null everywhere, because "
    "no coordinate file is connected.",
    "Resolution is carried per entry and never averaged, and no cutoff is "
    "applied. A structure at 1.5 angstrom and one at 4.3 support different "
    "claims, and a threshold chosen after seeing the spread would be fitted to "
    "it.",
]
