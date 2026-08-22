"""Surface filter check against the two validation sets.

Twenty-eight proteins with clinical or trial precedent must survive the filter.
Twelve deliberately chosen negatives must not. The negatives are split across
three different reasons for rejection, and the split is asserted: if every
negative failed on the anchor rule alone, the other two gates could be broken and
nothing here would notice.
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
    # anchored in an organelle membrane, and admitted for a while because the
    # phrase "cell membrane" or "cell surface" appeared in a free-text note
    # rather than in a location statement. See NOTE_ONLY_REJECTS.
    "GOLM1", "MTLN",
]

# Rejected because the topology gate found no outward face, despite being
# anchored in a membrane.
TOPOLOGY_REJECTS = ["CANX", "SEC61A1", "RPN1"]

# Rejected because the localisation gate now reads location statements only.
# Both of these carry an organelle location and a note that happens to contain a
# plasma-membrane phrase: one describes transient trafficking ("Cycles via the
# cell surface and endosomes upon lumenal pH disruption"), the other describes
# lipid binding ("Preferentially binds to cardiolipin relative to other common
# cell membrane lipids"). Neither is a statement about where the protein rests.
# Asserted by name, in the same shape as the CALR check below, because a
# regression here is silent: the count would move by fourteen and nothing else
# would look wrong.
NOTE_ONLY_REJECTS = ["GOLM1", "MTLN"]


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
        "surface": 3466,
        "single_pass": 1446,
        "multi_pass": 1884,
        "gpi_anchored": 136,
    }
    # How the withheld set is subdivided for reporting. The boundary between
    # these two rests on a compartment vocabulary that had to be inferred rather
    # than read, and no protein's fate depends on which side it lands. Reported,
    # not gated.
    subdivision = {
        "internal_anchored": 1362,
        "compartment_unresolved": 534,
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
            elif rec.outward_note_only:
                reason = "outward evidence only in a free-text note"
            elif not rec.outward:
                reason = "no outward face (topology)"
            else:
                reason = "SURVIVED"
            print(
                f"  {'ok  ' if ok else 'FAIL'}  {gene:10s} attached={str(rec.attached):5s}"
                f" outward={str(rec.outward):5s}  rejected on: {reason}"
            )
    print(f"  {rejected}/{len(NEGATIVE_CONTROLS)} rejected   expected 12")

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

    note_only = [
        g
        for g in NOTE_ONLY_REJECTS
        if (r := by_gene.get(g)) is not None and r.outward_note_only
    ]
    print(
        f"  rejected on the location-statement rule: {sorted(note_only)}"
        f"   expected {sorted(NOTE_ONLY_REJECTS)}"
    )

    # Every entry held out in the third state, by name. Fourteen is small enough
    # to read, and a silent change in this set is exactly what a count alone
    # would hide.
    held = sorted(r.gene for r in records if r.outward_note_only)
    print(f"\n  held out, plasma-membrane evidence only in a note ({len(held)}):")
    print(f"    {', '.join(held)}")
    print("    ASTN2 is the one to keep in view: its note reads \"Integral membrane")
    print("    protein not detected at the cell membrane\", so the previous rule")
    print("    admitted it on a sentence denying the thing it was matching.")

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
        and sorted(note_only) == sorted(NOTE_ONLY_REJECTS)
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
