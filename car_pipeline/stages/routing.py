"""Stage 4a — the risk profile selects the architecture, not the other way round."""

from __future__ import annotations

from dataclasses import dataclass


CONVENTIONAL = "CONVENTIONAL"
AND_GATE = "AND_GATE"
ADAPTOR = "ADAPTOR"


AND_NOT = "AND_NOT"
SWITCHABLE = "SWITCHABLE"

NO_ARCHITECTURE = "NO_ARCHITECTURE"
NOT_IMPLEMENTED = "NOT_IMPLEMENTED"
NOT_CONFIGURED = "NOT_CONFIGURED"

BUILT = (CONVENTIONAL, AND_GATE, ADAPTOR)

UNBUILT_REASON = {
    AND_NOT: "no exclusion-antigen source: pairing selects for tumour "
             "co-expression and an inhibitory CAR needs the opposite relation",
    SWITCHABLE: "the FKBP12 in the build is wild-type, so a rapamycin ON-switch "
                "and the mandatory rimiducid suicide switch answer to the same "
                "drug; the fix is a point mutation, which has no provenance "
                "class",
}


@dataclass(frozen=True)
class Route:
    """One target's architecture, and the ceiling that decided it."""

    gene: str
    architecture: str

    risk: float | None
    risk_organ: str | None

    ceiling: float | None

    exposure: str | None
    partner: str | None
    reason: str

    @property
    def admitted(self) -> bool:
        """Whether this route reached an architecture that can be built."""
        return self.architecture in BUILT


@dataclass(frozen=True)
class Tolerances:
    """The two declared ceilings."""

    persistent: float
    terminable: float | None = None


def route(
    gene: str,
    risk: float | None,
    risk_organ: str | None,
    tolerances: Tolerances,
    pair_risk: float | None = None,
    partner: str | None = None,
    exclusion_marker: str | None = None,
) -> Route:
    """The first architecture that admits this target, simplest first."""
    def made(architecture, ceiling, exposure, reason, prt=None):
        """Build a Route carrying this gene's risk and the reason for the verdict."""
        return Route(gene=gene, architecture=architecture, risk=risk,
                     risk_organ=risk_organ, ceiling=ceiling, exposure=exposure,
                     partner=prt, reason=reason)

    if risk is None:
        return made(NO_ARCHITECTURE, None, None,
                    "normal tissue risk not measured")

    if risk <= tolerances.persistent:
        return made(CONVENTIONAL, tolerances.persistent, "persistent",
                    f"risk {risk:.4f} within the persistent ceiling")

    if pair_risk is not None and pair_risk <= tolerances.persistent:
        return made(AND_GATE, tolerances.persistent, "persistent",
                    f"alone {risk:.4f} exceeds the persistent ceiling; paired "
                    f"{pair_risk:.4f} clears it", partner)

    if tolerances.terminable is None:
        adaptor_state = (NOT_CONFIGURED,
                         "no terminable ceiling declared for this project")
    elif risk <= tolerances.terminable:
        return made(ADAPTOR, tolerances.terminable, "terminable",
                    f"risk {risk:.4f} exceeds the persistent ceiling but is "
                    f"within the terminable one; exposure is dose-limited")
    else:
        adaptor_state = None

    if exclusion_marker:
        return made(AND_NOT, None, None,
                    f"{NOT_IMPLEMENTED}: {UNBUILT_REASON[AND_NOT]}")

    if adaptor_state:
        return made(adaptor_state[0], None, None, adaptor_state[1])

    return made(NO_ARCHITECTURE, None, None,
                f"risk {risk:.4f} exceeds every declared ceiling"
                + (f"; best pair {pair_risk:.4f} also exceeds the persistent "
                   "ceiling" if pair_risk is not None else "; no measured pair"))


def sweep(routes_input, tolerances: Tolerances, values) -> dict[float, int]:
    """How many targets the adaptor row admits across a range of ceilings."""
    out: dict[float, int] = {}
    for value in values:
        probe = Tolerances(persistent=tolerances.persistent, terminable=value)
        out[value] = sum(
            1 for item in routes_input
            if route(*item[:3], probe, *item[3:]).architecture == ADAPTOR
        )
    return out


def configuration_payload(tolerances: Tolerances) -> dict:
    """What routing contributes to the Stage 4 configuration hash."""
    return {
        "routing_version": 1,
        "persistent_ceiling": tolerances.persistent,
        "terminable_ceiling": tolerances.terminable,
        "built": list(BUILT),
    }
