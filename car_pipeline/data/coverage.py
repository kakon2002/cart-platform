"""Joins every source onto the surface set and records what is known about each."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

PROTEIN_CONFIRMED = "PROTEIN_CONFIRMED"
RNA_SUPPORTED = "RNA_SUPPORTED"
DATA_INSUFFICIENT = "DATA_INSUFFICIENT"

JOIN_ACCESSION = "accession"
JOIN_SYMBOL = "symbol"
JOIN_ENSEMBL = "ensembl"
JOIN_ENSEMBL_BRIDGE = "ensembl_bridge"

SILENT_TPM = 1.0


@dataclass
class CoverageRow:
    accession: str
    gene: str
    membrane_class: str | None
    evidence_class: str = DATA_INSUFFICIENT
    join_paths: dict[str, str] = field(default_factory=dict)
    has_staining: bool = False
    at_plasma_membrane: bool = False
    has_normal_baseline: bool = False
    has_tumour_expression: bool = False
    silent_in_every_normal_tissue: bool = False

    @property
    def bridged(self) -> bool:
        """True when some source was reachable only after the symbol failed."""
        return JOIN_ENSEMBL_BRIDGE in self.join_paths.values()


_IMMUNE_LOCUS = re.compile(
    r"""^(
        HLA-.+          # histocompatibility loci
      | KIR[23]D\w*     # killer immunoglobulin-like receptors
      | TR[ABGD]        # antigen receptor loci, unsegmented
      | TR[ABGD][CVDJ]\d*   # and their constant, variable, diversity, joining segments
    )$""",
    re.VERBOSE,
)


def _immune_locus(gene: str) -> bool:
    """Whether the symbol belongs to an immunoglobulin or MHC locus."""
    return bool(_IMMUNE_LOCUS.match(gene))


def _endogenous_retroviral(gene: str, protein_name: str) -> bool:
    """Whether the entry is an endogenous retroviral product."""
    return gene.startswith("ERV") or "endogenous retrovirus" in protein_name.lower()


def build_coverage(
    surface,
    atlas_by_accession: dict,
    atlas_by_symbol: dict,
    gtex_profiles: dict,
    tcga_join: dict,
) -> list[CoverageRow]:
    """One row per surface protein. The list is the same length as the input."""
    rows: list[CoverageRow] = []

    for rec in surface:
        row = CoverageRow(
            accession=rec.accession,
            gene=rec.gene,
            membrane_class=rec.membrane_class,
        )

        entry = atlas_by_accession.get(rec.accession)
        if entry is not None:
            row.join_paths["tissue_atlas"] = JOIN_ACCESSION
        elif rec.gene:
            entry = atlas_by_symbol.get(rec.gene)
            if entry is not None:
                row.join_paths["tissue_atlas"] = JOIN_SYMBOL
        if entry is not None:
            row.has_staining = entry.has_staining
            row.at_plasma_membrane = entry.at_plasma_membrane

        profile = gtex_profiles.get(rec.accession)
        if profile is not None:
            row.has_normal_baseline = True
            row.join_paths["normal_baseline"] = profile.join_path
            row.silent_in_every_normal_tissue = profile.silent_everywhere(SILENT_TPM)

        tumour = tcga_join.get(rec.accession)
        if tumour is not None:
            row.has_tumour_expression = True
            row.join_paths["tumour_expression"] = tumour[1]

        if row.has_staining:
            row.evidence_class = PROTEIN_CONFIRMED
        elif row.has_normal_baseline or row.has_tumour_expression:
            row.evidence_class = RNA_SUPPORTED
        else:
            row.evidence_class = DATA_INSUFFICIENT

        rows.append(row)

    return rows


def summarise(rows: Iterable[CoverageRow]) -> dict:
    """Counts per evidence class and join path across the coverage rows."""
    rows = list(rows)
    counts = {
        PROTEIN_CONFIRMED: 0,
        RNA_SUPPORTED: 0,
        DATA_INSUFFICIENT: 0,
    }
    for r in rows:
        counts[r.evidence_class] += 1
    silent = sum(
        1
        for r in rows
        if r.evidence_class == RNA_SUPPORTED and r.silent_in_every_normal_tissue
    )
    return {
        "total": len(rows),
        **counts,
        "rna_supported_but_silent": silent,
        "bridged": sum(1 for r in rows if r.bridged),
    }


def categorise_insufficient(rows: Iterable[CoverageRow], protein_names: dict) -> dict:
    """Break down what the unresolved group actually contains."""
    immune = ervs = other = 0
    examples: list[str] = []
    for r in rows:
        if r.evidence_class != DATA_INSUFFICIENT:
            continue
        name = protein_names.get(r.accession, "")
        if _immune_locus(r.gene):
            immune += 1
        elif _endogenous_retroviral(r.gene, name):
            ervs += 1
        else:
            other += 1
            if len(examples) < 10:
                examples.append(r.gene or r.accession)
    return {
        "polymorphic_immune_loci": immune,
        "retroviral_envelopes": ervs,
        "other": other,
        "other_examples": examples,
    }


if __name__ == "__main__":
    from car_pipeline.data.gtex import GTExSource
    from car_pipeline.data.hpa import HPASource, index as atlas_index
    from car_pipeline.data.tcga import TCGASource, match_surface as tcga_match
    from car_pipeline.data.uniprot import load_surface

    surface, _ = load_surface()
    by_acc, by_sym = atlas_index(HPASource().load())
    gtex_profiles, _, _ = GTExSource().match_surface(surface, by_acc)
    cohort = TCGASource().load()
    tcga_join = tcga_match(cohort, surface, by_acc)

    rows = build_coverage(surface, by_acc, by_sym, gtex_profiles, tcga_join)
    stats = summarise(rows)

    expected = {
        PROTEIN_CONFIRMED: 1945,
        RNA_SUPPORTED: 1491,
        DATA_INSUFFICIENT: 60,
    }
    print("evidence classes")
    for key, exp in expected.items():
        got = stats[key]
        pct = abs(got - exp) / exp * 100
        print(f"  {key}: {got:,}  expected {exp:,}  ({pct:.2f}% off)")

    total = stats["total"]
    summed = sum(stats[k] for k in expected)
    print(f"\n  sum of classes: {summed:,}   surface set: {total:,}"
          f"   {'nothing dropped' if summed == total else 'MISMATCH'}")

    print(f"\n  measured and silent in every normal tissue: "
          f"{stats['rna_supported_but_silent']:,}  expected 472")
    print(f"  reached only after a symbol failed: {stats['bridged']:,}")

    names = {r.accession: r.protein_name for r in surface}
    breakdown = categorise_insufficient(rows, names)
    print("\nunresolved group")
    for key, exp in [
        ("polymorphic_immune_loci", 14),
        ("retroviral_envelopes", 13),
        ("other", 33),
    ]:
        print(f"  {key}: {breakdown[key]}  expected {exp}")
    print(f"  a few of the others: {', '.join(breakdown['other_examples'])}")
