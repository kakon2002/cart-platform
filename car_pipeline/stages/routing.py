"""Stage 4a — the risk profile selects the architecture, not the other way round.

Before this module the pipeline applied one risk ceiling at Stage 3, before any
architecture was known. That is backwards: the architectures that exist to make
a risky target tolerable were only ever offered to targets that had already
cleared without them. 199 of 200 pool members died at that gate.

Here the profile routes first and the ceiling follows from the architecture.

**Two ceilings, never blended.** The existing ceiling governs a *persistent*
exposure — a self-amplifying T cell that cannot be withdrawn. An adaptor design
does not make the antigen safer; the adaptor still binds it. It makes the
exposure *terminable*, because activation needs a separately dosed protein with
a finite half-life. Magnitude and reversibility are different axes and get
different numbers.

**Routing never reduces a risk number.** The target's risk is carried through
unchanged; what routing changes is which ceiling it is compared against.
Substituting the receptor's risk for the target's would make the gate vacuous,
because every adaptor receptor binds the same tag and looks equally harmless.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Architectures this stage can route to. Ordered by increasing product
#: complexity, which is also the order they are tried in: one receptor and one
#: product before two receptors, and two receptors before two manufactured
#: products. A stated preference, not a tuned one.
CONVENTIONAL = "CONVENTIONAL"
AND_GATE = "AND_GATE"
ADAPTOR = "ADAPTOR"

#: Rows the spec names that this stage does not build. These are reported by
#: name with a reason. A target that would have routed to an unbuilt row is a
#: different finding from one that routes nowhere, and collapsing the two would
#: understate exactly what the missing architectures are worth.
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
    #: The Stage 3 risk, carried through unchanged. Never rewritten by routing.
    risk: float | None
    risk_organ: str | None
    #: The ceiling this target was actually compared against.
    ceiling: float | None
    #: Which of the two tolerances that ceiling came from.
    exposure: str | None
    partner: str | None
    reason: str

    @property
    def admitted(self) -> bool:
        return self.architecture in BUILT


@dataclass(frozen=True)
class Tolerances:
    """The two declared ceilings.

    Both are policy inputs, not measurements. The platform cannot derive a
    clinical risk tolerance from expression data, and a default for the
    terminable one would be this code quietly setting clinical policy — hence
    `terminable` is optional and its absence disables the adaptor row rather
    than being filled in.
    """

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
    """The first architecture that admits this target, simplest first.

    `pair_risk` is the combined risk of the best measured partner, already
    computed by Stage 4; `exclusion_marker` is the normal-restricted antigen an
    inhibitory design would need, which nothing currently supplies.
    """
    def made(architecture, ceiling, exposure, reason, prt=None):
        return Route(gene=gene, architecture=architecture, risk=risk,
                     risk_organ=risk_organ, ceiling=ceiling, exposure=exposure,
                     partner=prt, reason=reason)

    if risk is None:
        # Unmeasured is not safe. It routes nowhere and says why, rather than
        # being compared against a ceiling it has no number for.
        return made(NO_ARCHITECTURE, None, None,
                    "normal tissue risk not measured")

    # 1. One receptor, one product. Preferred whenever it is admissible.
    if risk <= tolerances.persistent:
        return made(CONVENTIONAL, tolerances.persistent, "persistent",
                    f"risk {risk:.4f} within the persistent ceiling")

    # 2. Two receptors, one product. Admissible when the *pair* clears the same
    #    persistent ceiling: an AND gate does not make the exposure terminable,
    #    it makes activation conditional, so the ceiling does not move.
    if pair_risk is not None and pair_risk <= tolerances.persistent:
        return made(AND_GATE, tolerances.persistent, "persistent",
                    f"alone {risk:.4f} exceeds the persistent ceiling; paired "
                    f"{pair_risk:.4f} clears it", partner)

    # 3. Two products, one receptor. The antigen is no safer — the adaptor still
    #    binds it — but the exposure can be stopped, so a different declared
    #    tolerance applies.
    if tolerances.terminable is None:
        adaptor_state = (NOT_CONFIGURED,
                         "no terminable ceiling declared for this project")
    elif risk <= tolerances.terminable:
        return made(ADAPTOR, tolerances.terminable, "terminable",
                    f"risk {risk:.4f} exceeds the persistent ceiling but is "
                    f"within the terminable one; exposure is dose-limited")
    else:
        adaptor_state = None

    # 4. Rows this stage does not build, named rather than silently dropped.
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
    """How many targets the adaptor row admits across a range of ceilings.

    Required by criterion A9. The terminable ceiling is a number this pipeline
    cannot measure, so its effect is reported rather than argued: a reader sees
    what any choice would have bought instead of having to trust one.
    """
    out: dict[float, int] = {}
    for value in values:
        probe = Tolerances(persistent=tolerances.persistent, terminable=value)
        out[value] = sum(
            1 for item in routes_input
            if route(*item[:3], probe, *item[3:]).architecture == ADAPTOR
        )
    return out


def configuration_payload(tolerances: Tolerances) -> dict:
    """What routing contributes to the Stage 4 configuration hash.

    The declared ceilings belong in the hash: two runs under different
    tolerances are different experiments and must not compare as one.
    """
    return {
        "routing_version": 1,
        "persistent_ceiling": tolerances.persistent,
        "terminable_ceiling": tolerances.terminable,
        "built": list(BUILT),
    }
