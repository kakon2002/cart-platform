"""Runs the safety gate and tests it against the criteria fixed in the spec."""

from __future__ import annotations

import random
import sys

from car_pipeline.configs.pdac import PDAC_PROJECT
from car_pipeline.data.antibodies import AntibodySource
from car_pipeline.data.coverage import build_coverage
from car_pipeline.data.depmap import DepMapSource, gene_index
from car_pipeline.data.gtex import GTExSource
from car_pipeline.data.hpa import HPASource, index as atlas_index
from car_pipeline.data.singlecell import SingleCellSource, match_surface as cell_match
from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
from car_pipeline.data.trials import TrialSource
from car_pipeline.data.uniprot import load_surface
from car_pipeline.stages import (
    construct_safety, stage3, stage4, stage5, stage6, stage9)
from car_pipeline.stages.stage1 import build_spec


PINNED_TRIALS = ["MSLN", "CLDN18"]

PINNED_ORIGIN = {"Amatuximab": "chimeric", "Zolbetuximab": "chimeric"}


ALT_CODON = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA",
    "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTA", "K": "AAA",
    "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG",
    "Y": "TAT", "V": "GTA",
}


def re_encode(protein: str) -> str:
    """The same protein under a second, genuinely synonymous codon assignment."""
    return "".join(ALT_CODON[r] for r in protein) + stage6.STOP


class FakeSegment:
    def __init__(self, name, accession, start, end, residues):
        """A domain map entry, for the controls that need one."""
        self.name, self.accession = name, accession
        self.start_residue, self.end_residue = start, end
        self.residues = residues
        self.bp_start, self.bp_end = 0, residues * 3
        self.provenance = "proteome"


def main() -> int:
    """Run the safety-gate criteria."""
    print("loading upstream", flush=True)
    decisions, manifest = stage4.read_decisions(allow_unusable=True)

    records = stage5.load_or_retrieve(
        decisions, AntibodySource(), manifest["stage4_hash"])
    binders = {r.gene: r for r in records}

    surface, _ = load_surface()
    atlas = HPASource().load()
    by_acc, by_sym = atlas_index(atlas)
    gtex_profiles, gtex_tissues, _ = GTExSource().match_surface(surface, by_acc)
    cohort = TCGASource().load()
    cohort_join = tcga_match(cohort, surface, by_acc)
    cells = SingleCellSource().load()
    cell_index = cell_match(cells, surface, by_acc)
    dependency, _ = DepMapSource().load()
    dep_index = gene_index(dependency)
    coverage_rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, cohort_join)
    spec = build_spec(PDAC_PROJECT)
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    overrides = {o: ov.tier
                 for o, ov in spec.inputs.tissue_criticality_overrides.items()}
    calibration = stage3.calibrate_atlas_levels(
        surface, by_acc, by_sym, gtex_profiles, gtex_tissues,
        stage3.RiskModel(overrides=overrides))
    ranked, _model, _ = stage3.rank(
        coverage_rows, {r.accession: r for r in surface}, by_acc, by_sym,
        cells, cell_index, gtex_profiles, gtex_tissues, cohort, cohort_join,
        dependency, dep_index, overrides, ceiling, calibration)
    risks = {r.gene: (r.risk, r.risk_organ) for r in ranked if r.gene}

    stage5_hash = stage5.configuration_hash(
        manifest["stage4_hash"], [r.gene for r in records])
    stage6_hash = stage6.configuration_hash(
        stage5_hash, [d["gene"] for d in decisions])
    genes = [d["gene"] for d in decisions]
    trials = TrialSource(antigens=genes).load()
    constructs = stage6.build(decisions, binders)
    gated = stage9.gate(decisions, binders, risks, trials, ceiling,
                        constructs=constructs)
    by_gene = {g.gene: g for g in gated}

    print()
    print("=" * 72)
    print("REJECTION CRITERIA")
    print("=" * 72)
    tripped: list[str] = []

    checked: list[str] = []
    def criterion(cid: str, is_tripped: bool, detail: str) -> None:
        """Report one criterion and record it if it tripped."""
        print(f"  {'TRIPPED ' if is_tripped else 'clear   '} {cid}: {detail}")
        checked.append(cid)
        if is_tripped:
            tripped.append(cid)

    empty = [g for g in PINNED_TRIALS
             if g not in trials or trials[g].total == 0]
    criterion("S1", bool(empty),
              "registry returns trials for both pins ("
              + ", ".join(f"{g}={trials[g].total}" for g in PINNED_TRIALS
                          if g in trials) + ")"
              if not empty else f"no trials returned for {empty}")

    origin_bad = [f"{n}->{stage9.binder_origin(n)}"
                  for n, want in PINNED_ORIGIN.items()
                  if stage9.binder_origin(n) != want]
    criterion("S2", bool(origin_bad),
              "pinned names classify as expected ("
              + ", ".join(f"{n}={stage9.binder_origin(n)}" for n in PINNED_ORIGIN)
              + ")" if not origin_bad else "; ".join(origin_bad))

    missing_risk = [g.gene for g in gated
                    if binders.get(g.gene)
                    and (binders[g.gene].sequence or binders[g.gene].structure)
                    and g.risk is None]
    criterion("S3", bool(missing_risk),
              "every target with a binder carries a Stage 3 risk"
              if not missing_risk else f"no risk for {missing_risk[:5]}")

    terminable = spec.design_constraints.terminable_risk_ceiling
    by_row = {d["gene"]: d for d in decisions}
    admitted = [g for g in gated
                if g.risk is not None and g.verdict != stage9.BLOCKED]
    contradicts = [f"{g.gene} risk {g.risk:.4f} over its applied {g.ceiling}"
                   for g in admitted if g.risk > g.ceiling]

    unearned = [
        f"{g.gene} at {g.risk:.4f} against ceiling {g.ceiling}, route exposure "
        f"{by_row.get(g.gene, {}).get('route_exposure')!r}"
        for g in admitted
        if g.risk > ceiling
        and not (by_row.get(g.gene, {}).get("route_exposure") == "terminable"
                 and g.ceiling == terminable)
    ]
    on_terminable = [g.gene for g in admitted if g.ceiling == terminable]
    criterion("S4", bool(contradicts) or bool(unearned),
              f"no target escapes the ceiling applied to it: "
              f"{len(admitted) - len(on_terminable)} admitted against the "
              f"persistent {ceiling}, {len(on_terminable)} against the "
              f"terminable {terminable}, each on a route declaring the exposure "
              f"stoppable"
              if not (contradicts or unearned)
              else "; ".join((contradicts + unearned)[:5]))

    leaked = [g.gene for g in gated
              if g.epitope_immunogenicity != stage9.NOT_CONNECTED]
    criterion("S5", bool(leaked),
              "epitope immunogenicity is NOT_CONNECTED on every row"
              if not leaked else f"{len(leaked)} rows carry a value")

    expected_rows = manifest["pool_size"]
    out_genes = {g.gene for g in gated}
    criterion("S6",
              len(gated) != expected_rows or len(out_genes) != expected_rows,
              f"{len(gated)} rows and {len(out_genes)} distinct genes against the "
              f"{expected_rows} the Stage 4 manifest records")

    unflagged = [g.gene for g in gated
                 if g.trials_stopped and g.verdict == stage9.PASSES]
    criterion("S7", bool(unflagged),
              "every target with stopped trials is flagged or blocked"
              if not unflagged else f"{unflagged[:5]} passed with stopped trials")

    planted = "ATGC" * 6
    with_repeat = "GGCA" + planted + ("TTAC" * 5) + planted + "CCTA"
    noise = random.Random(20260904)
    clean = "".join(noise.choice("ACGT") for _ in range(400))
    found = construct_safety.direct_repeats(with_repeat)
    at = [f.at for f in found]
    criterion(
        "S8", len(found) != 1 or "5" not in at[0]
        or bool(construct_safety.direct_repeats(clean)),
        f"a planted {len(planted)} bp repeat is reported at {at[0] if at else 'nothing'}, "
        f"and a control with none reports "
        f"{len(construct_safety.direct_repeats(clean))}")

    donor = construct_safety.splice_sites("AAAAAA" + "GTAAGT" + "CCCCCC")
    acceptor = construct_safety.splice_sites("AAA" + "CTCTCTCTCT" + "CAG" + "GGG")
    quiet = construct_safety.splice_sites("AAAAAAAAAAAAAAAA")
    kinds = {f.kind for f in donor} | {f.kind for f in acceptor}
    criterion(
        "S9",
        "splice_donor" not in {f.kind for f in donor}
        or "splice_acceptor" not in {f.kind for f in acceptor}
        or bool(quiet),
        f"planted donor and acceptor are both found ({sorted(kinds)}), and a "
        f"control with neither reports {len(quiet)}")

    body = "ATG" + "GCC" * (construct_safety.ORF_MIN + 8) + "TAA"
    in_frame_two = "GG" + body
    in_frame_zero = body
    two = construct_safety.alternate_orfs(in_frame_two)
    zero = construct_safety.alternate_orfs(in_frame_zero)
    criterion(
        "S10", not two or bool(zero),
        f"a {construct_safety.ORF_MIN + 9}-codon reading frame planted out of "
        f"frame is reported ({len(two)}), and the same frame planted as the "
        f"coding frame is not ({len(zero)})")

    built = [c for c in constructs if c.amino_acid_sequence]
    drifted, unmoved = [], []
    for c in built:
        first = construct_safety.findings(c.amino_acid_sequence, c.dna, c.segments)
        second = construct_safety.findings(
            c.amino_acid_sequence, re_encode(c.amino_acid_sequence), c.segments)

        def key(items, basis):
            """The findings of one basis, as comparable tuples."""
            return sorted((f.kind, f.at, f.detail)
                          for f in items if f.basis == basis)

        if key(first, construct_safety.CODON_INVARIANT) != key(
                second, construct_safety.CODON_INVARIANT):
            drifted.append(c.gene)
        if key(first, construct_safety.MAP_SPECIFIC) == key(
                second, construct_safety.MAP_SPECIFIC):
            unmoved.append(c.gene)
    unlabelled = [f.kind for c in built
                  for f in construct_safety.findings(
                      c.amino_acid_sequence, c.dna, c.segments)
                  if f.basis not in (construct_safety.CODON_INVARIANT,
                                     construct_safety.MAP_SPECIFIC)]
    criterion(
        "S11", bool(drifted) or bool(unmoved) or bool(unlabelled) or not built,
        f"re-encoded under a synonymous codon table, all {len(built)} "
        f"constructs keep every CODON_INVARIANT finding and every one of them "
        f"moves at least one MAP_SPECIFIC finding; no finding is unlabelled"
        if built and not (drifted or unmoved or unlabelled) else
        f"invariant drifted on {drifted[:3]}; map-specific unmoved on "
        f"{unmoved[:3]}; unlabelled {unlabelled[:3]}"
        if built else "no construct assembled, so the labels were not exercised")

    twice = [FakeSegment("CD8A leader", "P01732", 1, 21, 21),
             FakeSegment("binder", "", None, None, 100),
             FakeSegment("CD8A leader", "P01732", 1, 21, 21)]
    once = [FakeSegment("CD8A leader", "P01732", 1, 21, 21),
            FakeSegment("CD8A hinge", "P01732", 138, 182, 45)]
    duplicated = construct_safety.repeated_parts(twice)
    shipped = [c.gene for c in built
               if construct_safety.repeated_parts(c.segments)]
    criterion(
        "S12",
        len(duplicated) != 1 or bool(construct_safety.repeated_parts(once)),
        f"a domain map repeating one part reports it once and a map repeating "
        f"none reports nothing; across {len(built)} shipping design(s) the "
        f"detector reports a repeated part on {shipped or 'none'}, which is "
        f"reported rather than required either way - a dual design repeats "
        f"the leader, hinge and transmembrane by construction, and this arm "
        f"gates nothing")

    before = [(c.gene, c.amino_acid_sequence, c.dna, len(c.segments), c.verdict)
              for c in constructs]
    for c in built:
        construct_safety.analyse(c.amino_acid_sequence, c.dna, c.segments)
    after = [(c.gene, c.amino_acid_sequence, c.dna, len(c.segments), c.verdict)
             for c in constructs]
    criterion("S13", before != after,
              f"the arm changed no sequence, no domain map and no verdict "
              f"across all {len(constructs)} constructs")

    assembled = {c.gene for c in built}
    carried = {g.gene for g in gated if g.construct_safety is not None}
    criterion("S14", carried != assembled,
              f"{len(carried)} safety record(s) carry a construct-safety "
              f"report, exactly the {len(assembled)} construct(s) that "
              f"assembled")

    print("=" * 72)
    print(f"  {len(checked) - len(tripped)}/{len(checked)} criteria clear")
    if tripped:
        print(f"\n  STOPPING: {', '.join(tripped)} tripped.")
        return 2

    print()
    print("=" * 72)
    print("WHAT THE GATE SAYS")
    print("=" * 72)
    counts: dict[str, int] = {}
    for g in gated:
        counts[g.verdict] = counts.get(g.verdict, 0) + 1
    for verdict in (stage9.PASSES, stage9.FLAGGED, stage9.BLOCKED, stage9.NO_GATE):
        n = counts.get(verdict, 0)
        print(f"    {verdict:22s} {n:4d}  ({n / len(gated):.0%})")

    reached = [g for g in gated if g.verdict in (stage9.PASSES, stage9.FLAGGED)]
    print()
    print("  WHICH QUESTIONS THIS RUN ACTUALLY EXERCISED")
    print("  " + "-" * 68)
    print(f"    reached the immunogenicity and trials questions   {len(reached)} "
          f"of {len(gated)}")
    print(f"    stopped at carried Stage 3 risk                   "
          f"{sum(1 for g in gated if g.verdict == stage9.BLOCKED)}")
    print(f"    stopped for want of a binder                      "
          f"{sum(1 for g in gated if g.verdict == stage9.NO_GATE)}")
    if not reached:
        print()
        print("    So this gate's own two questions decided nothing here. Every")
        print("    target was settled by the risk Stage 3 already measured, which")
        print("    is the correct precedence and leaves the new logic untested on")
        print("    live data. S7 passes vacuously and should be read that way.")
        print("    S2 does not: it pins the origin rule against known names rather")
        print("    than against this run's output, which is why it still has force.")
        print()
        print("    The origin table below is computed for every target carrying a")
        print("    binder, including blocked ones, so it is measured rather than")
        print("    vacuous — it simply did not change any verdict.")

    print()
    print("    PASSES_STATED_CHECKS is not a safety claim. It means three named")
    print("    questions failed to show a problem, and one of them — epitope-level")
    print("    immunogenicity — was not asked at all.")

    print()
    print("  The five known targets")
    print("  " + "-" * 68)
    for gene in ("CEACAM5", "CEACAM6", "CLDN18", "MSLN", "MUC1"):
        g = by_gene.get(gene)
        if g is None:
            continue
        risk = f"{g.risk:.4f}" if g.risk is not None else "n/a"
        print(f"    {gene:9s} {g.verdict:22s} risk {risk:>7s} ({g.risk_organ}) "
              f"binder {g.binder_name or '-'} [{g.binder_origin}]")
        print(f"      trials {g.trials_total:4d}, stopped {g.trials_stopped}")
        for reason in g.reasons:
            print(f"      - {reason}")

    origins: dict[str, int] = {}
    for g in gated:
        if g.binder_name:
            origins[g.binder_origin] = origins.get(g.binder_origin, 0) + 1
    print()
    print("  Binder origin across targets that have one, read from the name stem")
    for origin, n in sorted(origins.items(), key=lambda kv: -kv[1]):
        print(f"    {origin:22s} {n:4d}")
    print("    A naming convention, not a sequence measurement. Stage 5 records")
    print("    humanisation_state NOT_CONNECTED from the deposited taxonomy; this")
    print("    is a second, equally indirect reading reported beside it.")

    print()
    print(f"  configuration hash "
          f"{stage9.configuration_hash(stage6_hash, [g.gene for g in gated], ceiling)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
