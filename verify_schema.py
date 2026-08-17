"""Step 2 check: the indication config parses and lands in discovery mode."""

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.schemas.project import DiscoveryMode, ProjectInput

CHECKS: list[tuple[str, object, object]] = []


def check(label: str, got: object, expected: object) -> None:
    CHECKS.append((label, got, expected))


def main() -> int:
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

    # Unknown fields must be rejected, not absorbed.
    try:
        ProjectInput(
            cancer_type="x", malignancy_type="solid", target_antigens="MSLN"
        )
    except Exception:
        check("mistyped field rejected", True, True)
    else:
        check("mistyped field rejected", False, True)

    # A blank antigen from a form is an absent one.
    blank = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="   "
    )
    check("blank antigen -> None", blank.target_antigen, None)
    check("blank antigen -> mode B", blank.discovery_mode, DiscoveryMode.DISCOVER)

    # A supplied antigen flips the mode.
    supplied = ProjectInput(
        cancer_type="x", malignancy_type="solid", target_antigen="MSLN"
    )
    check("supplied antigen -> mode A", supplied.discovery_mode.value, "A")

    # A blank cancer type is refused.
    try:
        ProjectInput(cancer_type="   ", malignancy_type="solid")
    except Exception:
        check("blank cancer_type rejected", True, True)
    else:
        check("blank cancer_type rejected", False, True)

    # An override without a usable rationale is refused.
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
