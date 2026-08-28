"""Checks for the project definition and the specification stage 1 builds."""

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.schemas.project import DiscoveryMode, ProjectInput
from car_pipeline.stages.stage1 import build_spec

CHECKS: list[tuple[str, object, object]] = []


def check(label: str, got: object, expected: object) -> None:
    """Report one check and record it if it failed."""
    CHECKS.append((label, got, expected))


def main() -> int:
    """Run the schema criteria."""
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

    try:
        ProjectInput(
            cancer_type="x", malignancy_type="solid", target_antigens="MSLN"
        )
    except Exception:
        check("mistyped field rejected", True, True)
    else:
        check("mistyped field rejected", False, True)

    blank = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="   "
    )
    check("blank antigen -> None", blank.target_antigen, None)
    check("blank antigen -> mode B", blank.discovery_mode, DiscoveryMode.DISCOVER)

    supplied = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="MSLN"
    )
    check("supplied antigen -> mode A", supplied.discovery_mode.value, "A")

    try:
        ProjectInput(cancer_type="   ", malignancy_type="solid")
    except Exception:
        check("blank cancer_type rejected", True, True)
    else:
        check("blank cancer_type rejected", False, True)

    try:
        ProjectInput(
            cancer_type="x",
            malignancy_type="solid",
            tissue_criticality_overrides={"lung": {"tier": 3, "rationale": "short"}},
        )
    except Exception:
        check("override without rationale rejected", True, True)
    else:
        check("override without rationale rejected", False, True)

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

    available = sum(
        1
        for d in blocking
        if d.status is DatasetStatus.AVAILABLE
    )
    check(
        "resolved availability score",
        round(spec.data_availability_score, 3),
        round(available / len(blocking), 3),
    )

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

    failed = 0
    for label, got, expected in CHECKS:
        ok = got == expected
        failed += 0 if ok else 1
        print(f"  {'ok  ' if ok else 'FAIL'}  {label}: got {got!r}  expected {expected!r}")

    print()
    print(f"discovery mode: {p.discovery_mode.value}   expected: B")
    print(f"checks passed: {len(CHECKS) - failed}/{len(CHECKS)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
