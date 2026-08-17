"""Stage 3 — screen the surface universe and rank it.

Implements `specs/stage3-target-discovery.md`. Weights, thresholds and
saturation points are fixed there and are read from here, not tuned.

Three numbers leave this stage per protein and they are never combined:
attractiveness, normal tissue risk, and evidence confidence. The first is
tumour-side only. The second is a gate. The third describes how much
measurement stands behind the other two.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from car_pipeline.data.coverage import (
    DATA_INSUFFICIENT,
    PROTEIN_CONFIRMED,
    RNA_SUPPORTED,
    CoverageRow,
)
from car_pipeline.data.singlecell import (
    DROPOUT_EPSILON,
    ENDOTHELIAL,
    FIBROBLAST,
    IMMUNE,
    MALIGNANT,
)
from car_pipeline.data.tcga import PRIMARY_TUMOUR, SOLID_NORMAL

WEIGHT_SET_VERSION = "1.0.0"

C1 = "malignant_expression"
C2 = "malignant_vs_stroma"
C3 = "tumour_vs_normal"
C4 = "patient_prevalence"
C5 = "surface_accessibility"
C6 = "escape_resistance"

WEIGHTS: dict[str, float] = {
    C1: 0.25,
    C2: 0.20,
    C3: 0.25,
    C4: 0.15,
    C5: 0.10,
    C6: 0.05,
}

MINIMUM_MEASURED_WEIGHT = 0.40

#: Free parameters that set where each component stops rewarding more. Named
#: here so they can be perturbed as a set by the twelfth rejection criterion.
SATURATION: dict[str, float] = {
    "c1_expression": 100.0,
    "c2_ratio": 50.0,
    "c3_fold": 64.0,
    "c4_prevalence_tpm": 10.0,
    "c5_ectodomain_residues": 200.0,
    "c6_gene_effect": 1.0,
}

#: Beyond this ratio the two tumour-versus-normal denominators are treated as
#: telling different stories, and the row says so rather than picking one.
FOLD_DISAGREEMENT = 2.0

STROMAL_COMPARTMENTS = (FIBROBLAST, IMMUNE, ENDOTHELIAL)

# --------------------------------------------------------------------------
# tissue criticality
# --------------------------------------------------------------------------

TIER_WEIGHTS = {1: 1.0, 2: 0.6, 3: 0.3}

ORGAN_TIERS: dict[str, int] = {
    # tier 1
    "brain": 1, "heart": 1, "lung": 1, "liver": 1, "kidney": 1,
    "pancreas": 1, "vascular": 1,
    # tier 2
    "gi_tract": 2, "marrow_and_blood": 2, "bladder": 2, "endocrine": 2,
    "muscle": 2, "nerve": 2, "eye": 2, "mucosa": 2,
    # tier 3
    "skin": 3, "adipose": 3, "breast": 3, "reproductive": 3, "salivary": 3,
    "connective": 3,
}

#: Not normal tissue. Excluded from risk rather than mapped to an organ.
EXCLUDED_LABELS = frozenset(
    {"Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "N/A"}
)

#: Exact-match assignments. No substring or keyword test is used anywhere here:
#: substring matching is what put the adrenal at kidney criticality and folded
#: the kidney into the brain, and an exhaustive table cannot make either mistake.
BASELINE_ORGANS: dict[str, str] = {
    "Adipose_Subcutaneous": "adipose",
    "Adipose_Visceral_Omentum": "adipose",
    "Adrenal_Gland": "endocrine",
    "Artery_Aorta": "vascular",
    "Artery_Coronary": "vascular",
    "Artery_Tibial": "vascular",
    "Bladder": "bladder",
    "Brain_Amygdala": "brain",
    "Brain_Anterior_cingulate_cortex_BA24": "brain",
    "Brain_Caudate_basal_ganglia": "brain",
    "Brain_Cerebellar_Hemisphere": "brain",
    "Brain_Cerebellum": "brain",
    "Brain_Cortex": "brain",
    "Brain_Frontal_Cortex_BA9": "brain",
    "Brain_Hippocampus": "brain",
    "Brain_Hypothalamus": "brain",
    "Brain_Nucleus_accumbens_basal_ganglia": "brain",
    "Brain_Putamen_basal_ganglia": "brain",
    "Brain_Spinal_cord_cervical_c-1": "brain",
    "Brain_Substantia_nigra": "brain",
    "Breast_Mammary_Tissue": "breast",
    "Cervix_Ectocervix": "reproductive",
    "Cervix_Endocervix": "reproductive",
    "Colon_Sigmoid": "gi_tract",
    "Colon_Transverse": "gi_tract",
    "Colon_Transverse_Mixed_Cell": "gi_tract",
    "Colon_Transverse_Mucosa": "gi_tract",
    "Colon_Transverse_Muscularis": "gi_tract",
    "Esophagus_Gastroesophageal_Junction": "gi_tract",
    "Esophagus_Mucosa": "gi_tract",
    "Esophagus_Muscularis": "gi_tract",
    "Fallopian_Tube": "reproductive",
    "Heart_Atrial_Appendage": "heart",
    "Heart_Left_Ventricle": "heart",
    "Kidney_Cortex": "kidney",
    "Kidney_Medulla": "kidney",
    "Liver": "liver",
    "Liver_Hepatocyte": "liver",
    "Liver_Mixed_Cell": "liver",
    "Liver_Portal_Tract": "liver",
    "Lung": "lung",
    "Minor_Salivary_Gland": "salivary",
    "Muscle_Skeletal": "muscle",
    "Nerve_Tibial": "nerve",
    "Ovary": "reproductive",
    "Pancreas": "pancreas",
    "Pancreas_Acini": "pancreas",
    "Pancreas_Islets": "pancreas",
    "Pancreas_Mixed_Cell": "pancreas",
    "Pituitary": "endocrine",
    "Prostate": "reproductive",
    "Skin_Not_Sun_Exposed_Suprapubic": "skin",
    "Skin_Sun_Exposed_Lower_leg": "skin",
    "Small_Intestine_Terminal_Ileum": "gi_tract",
    "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate": "gi_tract",
    "Small_Intestine_Terminal_Ileum_Mixed_Cell": "gi_tract",
    "Spleen": "marrow_and_blood",
    "Stomach": "gi_tract",
    "Stomach_Mixed_Cell": "gi_tract",
    "Stomach_Mucosa": "gi_tract",
    "Stomach_Muscularis": "gi_tract",
    "Testis": "reproductive",
    "Thyroid": "endocrine",
    "Uterus": "reproductive",
    "Vagina": "reproductive",
    "Whole_Blood": "marrow_and_blood",
}

ATLAS_ORGANS: dict[str, str] = {
    "adipose tissue": "adipose",
    "adrenal gland": "endocrine",
    "appendix": "gi_tract",
    "bone marrow": "marrow_and_blood",
    "breast": "breast",
    "bronchus": "lung",
    "cartilage": "connective",
    "caudate": "brain",
    "cerebellum": "brain",
    "cerebral cortex": "brain",
    "cervix": "reproductive",
    "choroid plexus": "brain",
    "colon": "gi_tract",
    "dorsal raphe": "brain",
    "duodenum": "gi_tract",
    "endometrium": "reproductive",
    "endometrium 1": "reproductive",
    "endometrium 2": "reproductive",
    "epididymis": "reproductive",
    "esophagus": "gi_tract",
    "eye": "eye",
    "fallopian tube": "reproductive",
    "gallbladder": "gi_tract",
    "hair": "skin",
    "heart muscle": "heart",
    "hippocampus": "brain",
    "hypothalamus": "brain",
    "kidney": "kidney",
    "lactating breast": "breast",
    "liver": "liver",
    "lung": "lung",
    "lymph node": "marrow_and_blood",
    "nasopharynx": "mucosa",
    "oral mucosa": "mucosa",
    "ovary": "reproductive",
    "pancreas": "pancreas",
    "parathyroid gland": "endocrine",
    "pituitary gland": "endocrine",
    "placenta": "reproductive",
    "prostate": "reproductive",
    "rectum": "gi_tract",
    "retina": "eye",
    "salivary gland": "salivary",
    "seminal vesicle": "reproductive",
    "skeletal muscle": "muscle",
    "skin": "skin",
    "skin 1": "skin",
    "skin 2": "skin",
    "small intestine": "gi_tract",
    "smooth muscle": "muscle",
    "soft tissue 1": "connective",
    "soft tissue 2": "connective",
    "sole of foot": "skin",
    "spleen": "marrow_and_blood",
    "stomach 1": "gi_tract",
    "stomach 2": "gi_tract",
    "substantia nigra": "brain",
    "testis": "reproductive",
    "thymus": "marrow_and_blood",
    "thyroid gland": "endocrine",
    "tonsil": "marrow_and_blood",
    "urinary bladder": "bladder",
    "vagina": "reproductive",
}

#: Ordinal staining levels on a 0-1 axis. Known to be miscalibrated against the
#: transcript axis; see the open problem in the specification. Carried forward
#: unchanged so the defect stays visible rather than being half-corrected here.
ATLAS_LEVEL_SCORE = {0: 0.0, 1: 1 / 3, 2: 2 / 3, 3: 1.0}

BASELINE_TPM_SATURATION = 1000.0


def _clamp(x: float) -> float:
    return 0.0 if x < 0 else (1.0 if x > 1 else x)


@dataclass
class Component:
    value: float | None          # the score, or None when not measured
    raw: float | None = None     # the measurement it came from
    note: str = ""

    @property
    def measured(self) -> bool:
        return self.value is not None


@dataclass
class Ranked:
    accession: str
    gene: str
    evidence_class: str
    components: dict[str, Component]
    measured_weight: float
    composite: float | None
    below_evidence_floor: bool
    risk: float | None
    risk_organ: str | None
    cleared: bool
    confidence: float
    sources_disagree: bool
    fold_baseline: float | None = None
    fold_cohort: float | None = None
    bridged: bool = False
    tier_rank: int = 0

    def component_value(self, key: str) -> float | None:
        c = self.components.get(key)
        return c.value if c else None


# --------------------------------------------------------------------------
# component scoring
# --------------------------------------------------------------------------


def _score_c1(malignant: float | None, sat: float) -> Component:
    if malignant is None or math.isnan(malignant) or malignant <= DROPOUT_EPSILON:
        # At or below the silence threshold this assay is reporting its own
        # capture failure, not the protein. Unmeasured, never zero.
        return Component(None, malignant, "below capture threshold")
    return Component(
        _clamp(math.log10(1 + malignant) / math.log10(1 + sat)), malignant
    )


def _score_c2(malignant: float | None, stromal_peak: float | None, sat: float) -> Component:
    if malignant is None or stromal_peak is None:
        return Component(None, None, "no row")
    if math.isnan(malignant) or math.isnan(stromal_peak):
        return Component(None, None, "no row")
    if malignant <= DROPOUT_EPSILON or stromal_peak <= DROPOUT_EPSILON:
        # A ratio against a denominator this assay failed to capture is not a
        # specificity measurement, however large it comes out.
        return Component(None, None, "below capture threshold")
    ratio = malignant / stromal_peak
    return Component(_clamp(math.log10(ratio) / math.log10(sat)), ratio)


def _score_c3(
    tumour: float | None,
    baseline_pancreas: float | None,
    cohort_normal: float | None,
    sat: float,
) -> tuple[Component, float | None, float | None, bool]:
    if tumour is None or baseline_pancreas is None:
        # No row on one side. Genuinely unmeasured.
        return Component(None, None, "no row"), None, None, False

    def _fold(normal: float | None) -> float | None:
        """Tumour over normal, with a measured absence handled as such.

        A normal value of zero is a measurement, not a gap: it says the protein
        was looked for in normal tissue and not found, which is the single most
        favourable thing this comparison can report. Returning unmeasured here
        would discard the best evidence a target can have and, through the
        evidence floor, drop it out of the ranking entirely.
        """
        if normal is None:
            return None
        if normal > 0:
            return tumour / normal
        if tumour > 0:
            return math.inf  # present in tumour, below detection in normal
        return 0.0           # absent on both sides: no margin, still measured

    fold_baseline = _fold(baseline_pancreas)
    fold_cohort = _fold(cohort_normal)

    disagree = False
    if fold_cohort is not None and fold_baseline is not None:
        if math.isinf(fold_baseline) != math.isinf(fold_cohort):
            # One denominator sits below detection and the other does not, so
            # the two are describing different normal tissue.
            disagree = True
        elif math.isfinite(fold_baseline) and fold_baseline > 0:
            spread = fold_cohort / fold_baseline
            disagree = spread > FOLD_DISAGREEMENT or spread < 1 / FOLD_DISAGREEMENT

    if math.isinf(fold_baseline):
        score = 1.0
    elif fold_baseline > 0:
        score = _clamp(math.log2(fold_baseline) / math.log2(sat))
    else:
        score = 0.0
    return Component(score, fold_baseline), fold_baseline, fold_cohort, disagree


def _score_c4(values: np.ndarray | None, threshold: float) -> Component:
    if values is None or values.size == 0:
        return Component(None, None, "no row")
    frac = float(np.mean(values >= threshold))
    return Component(_clamp(frac), frac)


def _score_c5(
    residues: int | None, at_plasma_membrane: bool, lipid_anchored: bool, sat: float
) -> Component:
    """Accessibility.

    The lipid-anchored class reports zero annotated outward residues by
    construction — there is no transmembrane segment to annotate around — so for
    that class the residue term is excluded rather than read as zero. Otherwise
    an unannotated ectodomain is unmeasured, never imputed: a protein nobody
    annotated must not outrank one measured and found small.
    """
    if lipid_anchored:
        return Component(1.0 if at_plasma_membrane else 0.6, None, "anchor class")
    if residues is None:
        if at_plasma_membrane:
            return Component(None, None, "ectodomain unannotated")
        return Component(None, None, "ectodomain unannotated, no localisation")
    size_term = _clamp(residues / sat)
    confirm = 1.0 if at_plasma_membrane else 0.7
    return Component(_clamp(size_term * confirm), float(residues))


def _score_c6(effect: float | None, screened: int, sat: float) -> Component:
    if effect is None or math.isnan(effect) or screened == 0:
        return Component(None, None, "unscreened")
    return Component(_clamp(-effect / sat), effect)


# --------------------------------------------------------------------------
# risk
# --------------------------------------------------------------------------


def _baseline_score(tpm: float) -> float:
    return _clamp(math.log10(1 + tpm) / math.log10(1 + BASELINE_TPM_SATURATION))


@dataclass
class RiskModel:
    overrides: dict[str, int] = field(default_factory=dict)
    fall_through: set[str] = field(default_factory=set)

    def tier(self, organ: str) -> int:
        return self.overrides.get(organ, ORGAN_TIERS[organ])

    def weight(self, organ: str) -> float:
        return TIER_WEIGHTS[self.tier(organ)]

    def organ_for_baseline(self, label: str) -> str | None:
        if label in EXCLUDED_LABELS:
            return None
        organ = BASELINE_ORGANS.get(label)
        if organ is None:
            self.fall_through.add(f"baseline:{label}")
        return organ

    def organ_for_atlas(self, label: str) -> str | None:
        if label in EXCLUDED_LABELS:
            return None
        organ = ATLAS_ORGANS.get(label)
        if organ is None:
            self.fall_through.add(f"atlas:{label}")
        return organ


def compute_risk(
    model: RiskModel,
    atlas_gene,
    baseline_values: np.ndarray | None,
    baseline_tissues: list[str],
) -> tuple[float | None, str | None]:
    """Worst organ, weighted by criticality.

    Maximum rather than mean: one unmanageable organ disqualifies a target
    however clean the rest are. Returns None when nothing measured it at all,
    which fails the gate rather than reading as no risk.
    """
    per_organ: dict[str, float] = {}

    if atlas_gene is not None:
        for tissue, _cell_type, level in atlas_gene.staining:
            organ = model.organ_for_atlas(tissue)
            if organ is None:
                continue
            score = ATLAS_LEVEL_SCORE[level]
            if score > per_organ.get(organ, -1.0):
                per_organ[organ] = score

    if baseline_values is not None:
        for label, tpm in zip(baseline_tissues, baseline_values):
            organ = model.organ_for_baseline(label)
            if organ is None:
                continue
            score = _baseline_score(float(tpm))
            # Maximum of the two sources, never the minimum: a "not detected"
            # staining call must not cancel a positive transcript measurement
            # for the same organ. That would understate risk.
            if score > per_organ.get(organ, -1.0):
                per_organ[organ] = score

    if not per_organ:
        return None, None

    worst_organ = max(per_organ, key=lambda o: per_organ[o] * model.weight(o))
    return per_organ[worst_organ] * model.weight(worst_organ), worst_organ


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------


def _confidence(components: dict[str, Component], row: CoverageRow) -> float:
    measured = sum(WEIGHTS[k] for k, c in components.items() if c.measured)
    tier_bonus = {
        PROTEIN_CONFIRMED: 0.3,
        RNA_SUPPORTED: 0.15,
        DATA_INSUFFICIENT: 0.0,
    }[row.evidence_class]
    return round(min(1.0, measured * 0.7 + tier_bonus), 4)


def rank(
    coverage_rows: list[CoverageRow],
    surface_by_accession: dict,
    atlas_by_accession: dict,
    atlas_by_symbol: dict,
    cell_atlas,
    cell_index: dict[str, int],
    gtex_profiles: dict,
    gtex_tissues: list[str],
    cohort,
    cohort_join: dict,
    dependency,
    dependency_index: dict[str, int],
    overrides: dict[str, int],
    ceiling: float,
    saturation: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> tuple[list[Ranked], RiskModel]:
    sat = dict(SATURATION if saturation is None else saturation)
    wts = dict(WEIGHTS if weights is None else weights)
    model = RiskModel(overrides=overrides)

    pancreas_cols = [
        i for i, t in enumerate(gtex_tissues) if BASELINE_ORGANS.get(t) == "pancreas"
    ]

    primary_mask = cohort.sample_types == PRIMARY_TUMOUR
    normal_mask = cohort.sample_types == SOLID_NORMAL

    out: list[Ranked] = []
    for row in coverage_rows:
        rec = surface_by_accession[row.accession]
        gene = row.gene

        # cell atlas
        ci = cell_index.get(row.accession)
        malignant = stromal_peak = None
        if ci is not None:
            malignant = cell_atlas.compartment_value(MALIGNANT, ci)
            peaks = [cell_atlas.compartment_value(c, ci) for c in STROMAL_COMPARTMENTS]
            peaks = [p for p in peaks if not math.isnan(p)]
            stromal_peak = max(peaks) if peaks else None

        # cohort
        tumour_median = cohort_normal_median = None
        prevalence_values = None
        cj = cohort_join.get(row.accession)
        if cj is not None:
            col = cohort.values[:, cj[0]]
            prevalence_values = col[primary_mask]
            tumour_median = float(np.median(prevalence_values))
            if normal_mask.any():
                cohort_normal_median = float(np.median(col[normal_mask]))

        # baseline
        profile = gtex_profiles.get(row.accession)
        baseline_pancreas = None
        if profile is not None and pancreas_cols:
            baseline_pancreas = float(np.median(profile.values[pancreas_cols]))

        c3, fold_b, fold_c, disagree = _score_c3(
            tumour_median, baseline_pancreas, cohort_normal_median, sat["c3_fold"]
        )

        di = dependency_index.get(gene) if gene else None
        components = {
            C1: _score_c1(malignant, sat["c1_expression"]),
            C2: _score_c2(malignant, stromal_peak, sat["c2_ratio"]),
            C3: c3,
            C4: _score_c4(prevalence_values, sat["c4_prevalence_tpm"]),
            C5: _score_c5(
                rec.extracellular_residues,
                row.at_plasma_membrane,
                rec.gpi_anchored,
                sat["c5_ectodomain_residues"],
            ),
            C6: (
                _score_c6(
                    dependency.median_effect(di),
                    dependency.screened_lines(di),
                    sat["c6_gene_effect"],
                )
                if di is not None
                else Component(None, None, "no row")
            ),
        }

        measured_weight = sum(wts[k] for k, c in components.items() if c.measured)
        if measured_weight >= MINIMUM_MEASURED_WEIGHT:
            composite = sum(
                wts[k] * c.value for k, c in components.items() if c.measured
            ) / measured_weight
            below_floor = False
        else:
            composite = None
            below_floor = True

        atlas_gene = atlas_by_accession.get(row.accession) or (
            atlas_by_symbol.get(gene) if gene else None
        )
        risk, organ = compute_risk(
            model,
            atlas_gene,
            profile.values if profile is not None else None,
            gtex_tissues,
        )
        # Undefined risk is not low risk. It fails.
        cleared = risk is not None and risk <= ceiling

        out.append(
            Ranked(
                accession=row.accession,
                gene=gene,
                evidence_class=row.evidence_class,
                components=components,
                measured_weight=round(measured_weight, 4),
                composite=None if composite is None else round(composite, 4),
                below_evidence_floor=below_floor,
                risk=None if risk is None else round(risk, 4),
                risk_organ=organ,
                cleared=cleared,
                confidence=_confidence(components, row),
                sources_disagree=disagree,
                fold_baseline=fold_b,
                fold_cohort=fold_c,
                bridged=row.bridged,
            )
        )

    # rank within evidence class, scored rows first
    for tier in (PROTEIN_CONFIRMED, RNA_SUPPORTED, DATA_INSUFFICIENT):
        members = [r for r in out if r.evidence_class == tier]
        members.sort(
            key=lambda r: (r.composite is None, -(r.composite or 0.0), r.gene or "")
        )
        for i, r in enumerate(members, start=1):
            r.tier_rank = i

    return out, model


# --------------------------------------------------------------------------
# reproducibility header
# --------------------------------------------------------------------------


def _revision() -> str:
    try:
        rev = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=_repo_root(), timeout=30,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True, text=True, cwd=_repo_root(), timeout=30,
        ).stdout.strip()
    except Exception:
        return "unknown"
    return f"{rev}{'-dirty' if dirty else ''}"


def _repo_root() -> str:
    from pathlib import Path

    return str(Path(__file__).resolve().parents[2])


def configuration_hash(overrides: dict[str, int], ceiling: float) -> str:
    """Covers the tissue tables themselves, not only the tier assignments.

    Adding or moving a single label shifts thousands of risk scores. A hash that
    did not move with it would let two materially different experiments compare
    as the same one, which is worse than publishing no hash at all.
    """
    payload = {
        "weight_set_version": WEIGHT_SET_VERSION,
        "weights": WEIGHTS,
        "minimum_measured_weight": MINIMUM_MEASURED_WEIGHT,
        "saturation": SATURATION,
        "fold_disagreement": FOLD_DISAGREEMENT,
        "dropout_epsilon": DROPOUT_EPSILON,
        "tier_weights": TIER_WEIGHTS,
        "organ_tiers": ORGAN_TIERS,
        "baseline_organs": BASELINE_ORGANS,
        "atlas_organs": ATLAS_ORGANS,
        "excluded_labels": sorted(EXCLUDED_LABELS),
        "atlas_level_score": {str(k): v for k, v in ATLAS_LEVEL_SCORE.items()},
        "baseline_tpm_saturation": BASELINE_TPM_SATURATION,
        "overrides": overrides,
        "ceiling": ceiling,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def header(spec, model: RiskModel, universe: int, pins: dict) -> str:
    ceiling = spec.design_constraints.normal_tissue_risk_ceiling
    lines = [
        "=" * 72,
        "TARGET DISCOVERY - STAGE 3",
        "=" * 72,
        f"  project              {spec.project_id}",
        f"  indication           {spec.inputs.cancer_type}",
        f"  discovery mode       {spec.discovery_mode.value}",
        f"  target supplied      {spec.inputs.target_antigen!r}",
        f"  weight set           {WEIGHT_SET_VERSION}",
        f"  weights              " + ", ".join(f"{k}={v}" for k, v in WEIGHTS.items()),
        f"  minimum measured     {MINIMUM_MEASURED_WEIGHT}",
        f"  saturation           " + ", ".join(f"{k}={v}" for k, v in SATURATION.items()),
        f"  dropout epsilon      {DROPOUT_EPSILON}",
        f"  fold disagreement    {FOLD_DISAGREEMENT}x",
        f"  risk ceiling         {ceiling}",
        f"  universe             {universe}",
        f"  revision             {_revision()}",
        f"  configuration hash   {configuration_hash(model.overrides, ceiling)}",
        "  criticality tiers",
    ]
    for tier in (1, 2, 3):
        organs = sorted(o for o, t in ORGAN_TIERS.items() if t == tier)
        lines.append(f"    tier {tier} (w={TIER_WEIGHTS[tier]}): {', '.join(organs)}")
    if model.overrides:
        lines.append("  overrides")
        for organ, tier in model.overrides.items():
            lines.append(f"    {organ} -> tier {tier} (w={TIER_WEIGHTS[tier]})")
            rationale = spec.inputs.tissue_criticality_overrides.get(organ)
            if rationale is not None:
                lines.append(f"      rationale: {rationale.rationale}")
    lines.append("  source pins")
    for name, pin in pins.items():
        lines.append(f"    {name}: {pin}")
    lines.append("=" * 72)
    return "\n".join(lines)
