"""Surface filter check against the two validation sets.

Twenty-eight proteins with clinical or trial precedent must survive the filter.
Ten deliberately chosen negatives must not. The negatives are split between two
different reasons for rejection, and the split is asserted: if every negative
failed on the anchor rule alone, the topology gate could be broken and nothing
here would notice.
"""

from car_pipeline.data.uniprot import UniProtSource, summarise

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

    print("counts")
    expected = {
        "entries": 20431,
        "surface": 3496,
        "single_pass": 1464,
        "multi_pass": 1894,
        "gpi_anchored": 138,
        "internal_anchored": 1322,
        "compartment_unresolved": 550,
    }
    for key, exp in expected.items():
        got = stats[key]
        flag = "ok  " if got == exp else "DIFF"
        print(f"  {flag}  {key}: {got:,}  expected {exp:,}")

    print("\nknown targets that must survive")
    passed = 0
    for gene in KNOWN_TARGETS:
        rec = by_gene.get(gene)
        ok = rec is not None and rec.is_surface
        passed += 1 if ok else 0
        if not ok:
            state = "absent" if rec is None else (
                f"attached={rec.attached} outward={rec.outward}"
            )
            print(f"  FAIL  {gene}: {state}")
    print(f"  {passed}/{len(KNOWN_TARGETS)} survived   expected 28")

    print("\nnegative controls that must not survive")
    rejected = 0
    for gene in NEGATIVE_CONTROLS:
        rec = by_gene.get(gene)
        ok = rec is None or not rec.is_surface
        rejected += 1 if ok else 0
        if not ok:
            print(f"  FAIL  {gene} survived the filter")
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

    count_ok = all(stats[k] == v for k, v in expected.items())
    sets_ok = (
        passed == len(KNOWN_TARGETS)
        and rejected == len(NEGATIVE_CONTROLS)
        and bool(on_topology)
        and calr is not None
        and calr.outward
        and not calr.attached
    )
    print(f"\nvalidation sets: {'pass' if sets_ok else 'FAIL'}")
    print(f"counts match spec: {'yes' if count_ok else 'no'}")
    return 0 if (count_ok and sets_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
