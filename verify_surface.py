"""Surface filter check against the two validation sets.

Twenty-eight proteins with clinical or trial precedent must survive the filter.
Ten deliberately chosen negatives must not. The negatives are split between two
different reasons for rejection, and the split is asserted: if every negative
failed on the anchor rule alone, the topology gate could be broken and nothing
here would notice.
"""

from car_pipeline.data.uniprot import UniProtSource, summarise

# These figures are now measured against the pinned proteome release rather than
# reconstructed, so they are reproducible and a difference means the filter or
# the data moved. The tolerance is kept because the pin is enforced on the
# response rather than on the request: a release bump is caught by the fetch,
# but nothing here should fail on the last digit of a count. The validation sets
# below, which are exact, remain the real test.
COUNT_TOLERANCE_PCT = 3.0

KNOWN_TARGETS = [
    "MSLN", "CD19", "CD22", "EGFR", "ERBB2", "CEACAM5", "EPCAM", "PSCA",
    "FOLH1", "MUC1", "ROR1", "GPC3", "CD33", "TNFRSF17", "IL13RA2", "L1CAM",
    "CLDN18", "TSHR", "DLL3", "CD70", "TACSTD2", "NCAM1", "CD274", "MET",
    "ALPP", "CD52", "CD59", "THY1",
]

NEGATIVE_CONTROLS = [
    # nuclear and cytoskeletal
    "TP53", "ACTB", "GAPDH", "TUBB", "LMNA",
    # membrane, but of internal compartments
    "CANX", "CALR", "GOLGA2", "SEC61A1", "RPN1",
]

# Rejected because the topology gate found no outward face, despite being
# anchored in a membrane.
TOPOLOGY_REJECTS = ["CANX", "SEC61A1", "RPN1"]


def main() -> int:
    records = UniProtSource().load()
    by_gene = {}
    for rec in records:
        if rec.gene and rec.gene not in by_gene:
            by_gene[rec.gene] = rec

    stats = summarise(records)

    # What the two gates decide. These are gated on.
    decisions = {
        "entries": 20431,
        "surface": 3480,
        "single_pass": 1455,
        "multi_pass": 1889,
        "gpi_anchored": 136,
    }
    # How the withheld set is subdivided for reporting. The boundary between
    # these two rests on a compartment vocabulary that had to be inferred rather
    # than read, and no protein's fate depends on which side it lands. Reported,
    # not gated.
    subdivision = {
        "internal_anchored": 1349,
        "compartment_unresolved": 533,
    }

    drift_ok = True
    print("filter decisions")
    for key, exp in decisions.items():
        got = stats[key]
        pct = abs(got - exp) / exp * 100
        drift_ok = drift_ok and pct <= COUNT_TOLERANCE_PCT
        flag = "ok  " if pct <= COUNT_TOLERANCE_PCT else "WIDE"
        print(f"  {flag}  {key}: {got:,}  reference {exp:,}  ({pct:.2f}% off)")

    print("\nwithheld subdivision (reported, not gated)")
    for key, exp in subdivision.items():
        got = stats[key]
        pct = abs(got - exp) / exp * 100
        print(f"        {key}: {got:,}  reference {exp:,}  ({pct:.2f}% off)")

    print("\nknown targets that must survive")
    passed = 0
    for gene in KNOWN_TARGETS:
        rec = by_gene.get(gene)
        ok = rec is not None and rec.is_surface
        passed += 1 if ok else 0
        if rec is None:
            print(f"  FAIL  {gene:10s} absent from the proteome")
        else:
            ecto = (
                "unannotated"
                if rec.extracellular_residues is None
                else f"{rec.extracellular_residues} aa"
            )
            print(
                f"  {'ok  ' if ok else 'FAIL'}  {gene:10s} {rec.membrane_class or '-':14s}"
                f" tm={rec.transmem_count:<2d} anchor={'yes' if rec.gpi_anchored else 'no ':3s}"
                f" ectodomain={ecto}"
            )
    print(f"  {passed}/{len(KNOWN_TARGETS)} survived   expected 28")

    print("\nnegative controls that must not survive")
    rejected = 0
    for gene in NEGATIVE_CONTROLS:
        rec = by_gene.get(gene)
        ok = rec is None or not rec.is_surface
        rejected += 1 if ok else 0
        if rec is None:
            print(f"  ok    {gene:10s} absent from the proteome")
        else:
            if not rec.attached:
                reason = "no anchor"
            elif not rec.outward:
                reason = "no outward face (topology)"
            else:
                reason = "SURVIVED"
            print(
                f"  {'ok  ' if ok else 'FAIL'}  {gene:10s} attached={str(rec.attached):5s}"
                f" outward={str(rec.outward):5s}  rejected on: {reason}"
            )
    print(f"  {rejected}/{len(NEGATIVE_CONTROLS)} rejected   expected 10")

    print("\nrejection reasons")
    on_topology = [
        g
        for g in TOPOLOGY_REJECTS
        if (r := by_gene.get(g)) is not None and r.attached and not r.outward
    ]
    print(
        f"  rejected on topology rather than anchor: {sorted(on_topology)}"
        f"   expected at least one of {TOPOLOGY_REJECTS}"
    )

    calr = by_gene.get("CALR")
    if calr is not None:
        print(
            f"  CALR: outward={calr.outward} attached={calr.attached}"
            "   expected outward=True attached=False"
        )

    sets_ok = (
        passed == len(KNOWN_TARGETS)
        and rejected == len(NEGATIVE_CONTROLS)
        and bool(on_topology)
        and calr is not None
        and calr.outward
        and not calr.attached
    )
    print(f"\nvalidation sets: {'pass' if sets_ok else 'FAIL'}")
    print(
        f"filter decisions within {COUNT_TOLERANCE_PCT:.0f}% of the reference: "
        f"{'yes' if drift_ok else 'NO'}"
    )
    return 0 if (drift_ok and sets_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
