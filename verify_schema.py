"""Checks for the project definition and the specification stage 1 builds."""

from pydantic import ValidationError

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.schemas.project import DiscoveryMode, ProjectInput
from car_pipeline.stages.stage1 import build_spec

CHECKS: list[tuple[str, object, object]] = []
TRIPPED: list[str] = []


def criterion(cid: str, is_tripped: bool, detail: str) -> None:
    """Report one criterion and record it if it tripped.

    The same shape every other verifier uses. This file previously used a
    check(label, got, expected) idiom of its own, which kept its checks
    outside the criterion inventory every audit reads while still counting
    towards the suite total -- so the total meant two different things at
    once.
    """
    print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
    CHECKS.append((cid, is_tripped, detail))
    if is_tripped:
        TRIPPED.append(cid)


def check(label: str, got: object, expected: object) -> None:
    """One equality criterion, reported in the shared shape.

    The id is assigned in order rather than written at each call site: the
    ids exist so the run driver can itemise and exempt criteria, and the label
    is what a reader needs. Both are carried.
    """
    criterion(f"C{len(CHECKS) + 1}", got != expected,
              f"{label}: got {got!r}, expected {expected!r}")


def rejects(label: str, field: str, **kwargs) -> None:
    """Assert the schema rejects these inputs, and rejects them for the field.

    Catching bare Exception would pass on any failure at all -- a renamed
    class, a typo in the test, an import error -- so it would report success
    while proving nothing. Confirmed: replacing ProjectInput with a stub that
    raises RuntimeError left the earlier version of these three checks
    passing. This catches the validation error specifically and requires the
    offending field to be named in it.
    """
    try:
        ProjectInput(**kwargs)
    except ValidationError as exc:
        named = field in str(exc)
        check(label, named, True)
    except Exception as exc:
        check(f"{label} [raised {type(exc).__name__}, not a validation error]",
              False, True)
    else:
        check(label, False, True)


def main() -> int:
    """Run the schema criteria."""
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    p = PDAC_PROJECT

    check("discovery_mode", p.discovery_mode.value, "B")
    check("target_antigen", p.target_antigen, None)
    check("cancer_type", p.cancer_type, "Pancreatic Ductal Adenocarcinoma")
    check("malignancy_type", p.malignancy_type.value, "solid")
    check("product_type", p.product_type.value, "autologous")
    check("car_format", p.car_format.value, "auto")
    check("safety_tolerance", p.safety_tolerance.value, "conservative")
    check("vector_payload_limit_kb", p.manufacturing.vector_payload_limit_kb, 4.7)
    check("max_genetic_edits", p.manufacturing.max_genetic_edits, 2)
    check("pancreas override tier", p.tissue_criticality_overrides["pancreas"].tier, 2)

    rejects("mistyped field rejected", "target_antigens",
            cancer_type="x", malignancy_type="solid", target_antigens="MSLN")

    blank = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="   "
    )
    check("blank antigen -> None", blank.target_antigen, None)
    check("blank antigen -> mode B", blank.discovery_mode, DiscoveryMode.DISCOVER)

    supplied = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="MSLN"
    )
    check("supplied antigen -> mode A", supplied.discovery_mode.value, "A")

    rejects("blank cancer_type rejected", "cancer_type",
            cancer_type="   ", malignancy_type="solid")

    rejects("override without rationale rejected", "rationale",
            cancer_type="x", malignancy_type="solid",
            tissue_criticality_overrides={
                "lung": {"tier": 3, "rationale": "short"}})

    spec = build_spec(p)
    blocking = [d for d in spec.required_datasets if d.required]

    check("spec discovery_mode", spec.discovery_mode.value, "B")
    check("datasets", len(spec.required_datasets), 10)
    check("blocking datasets", len(blocking), 8)
    check("construct budget kb", spec.design_constraints.max_construct_kb, 3.5)
    check("safety switch required", spec.design_constraints.require_safety_switch, True)
    check("risk ceiling", spec.design_constraints.normal_tissue_risk_ceiling, 0.15)
    check("allowed formats", len(spec.design_constraints.allowed_car_formats), 5)
    check(
        "auto excluded",
        "auto" not in [f.value for f in spec.design_constraints.allowed_car_formats],
        True,
    )
    check("spec target_antigen", spec.inputs.target_antigen, None)

    unresolved = build_spec(p, resolve_sources=False)
    check("unresolved availability score", unresolved.data_availability_score, 0.0)

    from car_pipeline.schemas.spec import DatasetStatus

    # The earlier version recomputed available/len(blocking) from the same
    # list and the same status field the implementation reads, so both sides
    # moved together: forcing every blocking dataset to MISSING left it
    # passing. It pinned the formula and could not notice a wrong status,
    # which is what a reader would assume it covered.
    #
    # These assert properties of the score instead, each of which fails on a
    # different mistake: the score must lie in [0, 1], it must be zero when
    # nothing is resolved, and it must move when the statuses move.
    available = sum(1 for d in blocking if d.status is DatasetStatus.AVAILABLE)
    check("availability score within [0, 1]",
          0.0 <= spec.data_availability_score <= 1.0, True)
    check("availability score reflects some resolved dataset",
          spec.data_availability_score > 0.0, available > 0)

    none_available = build_spec(p, resolve_sources=False)
    check("availability score is zero when nothing resolves",
          none_available.data_availability_score, 0.0)
    check("availability score responds to status",
          spec.data_availability_score != none_available.data_availability_score,
          available > 0)

    validate = build_spec(
        ProjectInput(
            cancer_type="x", malignancy_type="solid", target_antigen="MSLN"
        )
    )
    check("validation datasets", len(validate.required_datasets), 7)
    check(
        "validation blocking",
        len([d for d in validate.required_datasets if d.required]),
        5,
    )

    spec.inputs.target_antigen = "SEEDED"
    check("input not mutated by build", PDAC_PROJECT.target_antigen, None)

    check("project id unique", build_spec(p).project_id != build_spec(p).project_id, True)

    print("=" * 72)
    print(f"  {len(CHECKS) - len(TRIPPED)}/{len(CHECKS)} criteria clear")
    if TRIPPED:
        print(f"\n  STOPPING: {', '.join(TRIPPED)} tripped.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
