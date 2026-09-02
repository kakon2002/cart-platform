"""Stage 3 — screen the surface universe and rank it."""

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
    JOIN_ENSEMBL_BRIDGE as CELL_JOIN_ENSEMBL_BRIDGE,
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

TUMOUR_DOMINANT = "TUMOUR_DOMINANT"
STROMA_DOMINANT = "STROMA_DOMINANT"
STROMA_UNRESOLVED = "STROMA_UNRESOLVED"

STROMA_RATIO_FLOOR = 1.0


SATURATION: dict[str, float] = {
    "c1_expression": 100.0,
    "c2_ratio": 50.0,
    "c3_fold": 64.0,
    "c4_prevalence_tpm": 10.0,
    "c5_ectodomain_residues": 200.0,
    "c6_gene_effect": 1.0,
    "c3_baseline_floor_tpm": 0.1,
}


FOLD_DISAGREEMENT = 2.0


FIELD_ELEVATION = "field_elevation"

STROMAL_COMPARTMENTS = (FIBROBLAST, IMMUNE, ENDOTHELIAL)


TIER_WEIGHTS = {1: 1.0, 2: 0.6, 3: 0.3}

ORGAN_TIERS: dict[str, int] = {
    "brain": 1, "heart": 1, "lung": 1, "liver": 1, "kidney": 1,
    "pancreas": 1, "vascular": 1, "eye": 1,
    "gi_tract": 2, "marrow_and_blood": 2, "bladder": 2, "endocrine": 2,
    "muscle": 2, "nerve": 2, "mucosa": 2,
    "skin": 3, "adipose": 3, "breast": 3, "reproductive": 3, "salivary": 3,
    "connective": 3,
}


PLATFORM_ADDED_ORGANS = frozenset({"vascular", "eye", "mucosa", "connective"})


EXCLUDED_LABELS = frozenset(
    {"Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "N/A"}
)


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

LEVEL_NAMES = {0: "Not detected", 1: "Low", 2: "Medium", 3: "High"}


@dataclass
class CalibrationCurve:
    """Staining levels placed on the transcript axis by measurement."""

    tpm: dict[int, float]
    quartiles: dict[int, tuple[float, float, float]]
    counts: dict[int, int]
    separations: dict[int, float]
    monotonic: bool
    observations: int

    def score(self, level: int) -> float:
        """This tissue's contribution at the given staining level."""
        if level == 0:
            return 0.0
        return _baseline_score(self.tpm[level])

    def as_payload(self) -> dict:
        """The per-level values, rounded for storage."""
        return {str(k): round(v, 6) for k, v in sorted(self.tpm.items())}


def _separation(lower: np.ndarray, upper: np.ndarray) -> float:
    """Probability that a draw from the upper level exceeds one from the lower."""
    combined = np.concatenate([lower, upper])
    order = combined.argsort(kind="mergesort")
    ranks = np.empty(len(combined), dtype=float)
    ranks[order] = np.arange(1, len(combined) + 1)
    ordered = combined[order]
    i = 0
    while i < len(ordered):
        j = i
        while j + 1 < len(ordered) and ordered[j + 1] == ordered[i]:
            j += 1
        if j > i:
            ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    u = ranks[len(lower) :].sum() - len(upper) * (len(upper) + 1) / 2
    return float(u / (len(lower) * len(upper)))


def calibrate_atlas_levels(
    surface,
    atlas_by_accession: dict,
    atlas_by_symbol: dict,
    gtex_profiles: dict,
    gtex_tissues: list[str],
    model: "RiskModel",
) -> CalibrationCurve:
    """Measure what each staining level is worth in transcript terms."""
    populations: dict[int, list[float]] = {k: [] for k in LEVEL_NAMES}

    for rec in surface:
        gene = atlas_by_accession.get(rec.accession) or (
            atlas_by_symbol.get(rec.gene) if rec.gene else None
        )
        profile = gtex_profiles.get(rec.accession)
        if gene is None or profile is None or not gene.staining:
            continue

        atlas_organ: dict[str, int] = {}
        for tissue, _cell_type, level in gene.staining:
            organ = model.organ_for_atlas(tissue)
            if organ is not None and level > atlas_organ.get(organ, -1):
                atlas_organ[organ] = level

        base_organ: dict[str, float] = {}
        for label, tpm in zip(gtex_tissues, profile.values):
            organ = model.organ_for_baseline(label)
            if organ is not None:
                value = float(tpm)
                if value > base_organ.get(organ, -1.0):
                    base_organ[organ] = value

        for organ, level in atlas_organ.items():
            if organ in base_organ:
                populations[level].append(base_organ[organ])

    empty = [LEVEL_NAMES[k] for k, v in populations.items() if not v]
    if empty:
        raise ValueError(
            "no paired observations for staining level(s): " + ", ".join(empty)
        )

    arrays = {k: np.asarray(v, dtype=float) for k, v in populations.items()}
    tpm: dict[int, float] = {}
    quartiles: dict[int, tuple[float, float, float]] = {}
    counts: dict[int, int] = {}
    for k, v in arrays.items():
        q1, med, q3 = (float(x) for x in np.percentile(v, [25, 50, 75]))
        tpm[k] = med
        quartiles[k] = (q1, med, q3)
        counts[k] = int(v.size)

    separations = {k: _separation(arrays[k], arrays[k + 1]) for k in range(3)}
    ordered_medians = [tpm[k] for k in sorted(tpm)]
    monotonic = all(
        ordered_medians[i] < ordered_medians[i + 1]
        for i in range(len(ordered_medians) - 1)
    )

    return CalibrationCurve(
        tpm=tpm,
        quartiles=quartiles,
        counts=counts,
        separations=separations,
        monotonic=monotonic,
        observations=int(sum(counts.values())),
    )

BASELINE_TPM_SATURATION = 1000.0


DEFAULT_MARGIN_LABEL = "Pancreas"


def _clamp(x: float) -> float:
    """Hold a score inside the unit interval."""
    return 0.0 if x < 0 else (1.0 if x > 1 else x)


def _finite(x: float | None) -> float | None:
    """Pass a real number through; turn anything else into an absent one."""
    if x is None or math.isnan(x) or math.isinf(x):
        return None
    return x


@dataclass
class Component:
    value: float | None
    raw: float | None = None
    note: str = ""

    @property
    def measured(self) -> bool:
        """Whether a value was measured at all."""
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

    field_elevation: float | None = None
    bridged: bool = False
    tier_rank: int = 0

    tumour_side_verdict: str = STROMA_UNRESOLVED
    protein_arm_measured: bool = True
    risk_basis: str = "staining and transcript"
    risk_is_lower_bound: bool = False
    composite_supported: float | None = None

    def component_value(self, key: str) -> float | None:
        """One component's score, or None where it was not measured."""
        c = self.components.get(key)
        return c.value if c else None


def _score_c1(malignant: float | None, sat: float) -> Component:
    """Malignant expression, unmeasured below the capture threshold."""
    if malignant is None or math.isnan(malignant) or malignant <= DROPOUT_EPSILON:
        return Component(None, malignant, "below capture threshold")
    return Component(
        _clamp(math.log10(1 + malignant) / math.log10(1 + sat)), malignant
    )


def _score_c2(malignant: float | None, stromal_peak: float | None, sat: float) -> Component:
    """Malignant signal against the strongest stromal compartment."""
    if malignant is None or stromal_peak is None:
        return Component(None, None, "no row")
    if math.isnan(malignant) or math.isnan(stromal_peak):
        return Component(None, None, "no row")
    if malignant <= DROPOUT_EPSILON or stromal_peak <= DROPOUT_EPSILON:
        return Component(None, None, "below capture threshold")
    ratio = malignant / stromal_peak
    return Component(_clamp(math.log10(ratio) / math.log10(sat)), ratio)


def _score_c3(
    tumour: float | None,
    baseline_pancreas: float | None,
    cohort_normal: float | None,
    sat: dict[str, float],
) -> tuple[Component, float | None, float | None]:
    """Tumour level against the matched normal baseline."""
    if tumour is None or baseline_pancreas is None:
        return Component(None, None, "no row"), None, None

    def _fold(normal: float | None) -> float | None:
        """Tumour over normal, with the denominator floored at detection."""
        if normal is None:
            return None
        return tumour / max(normal, floor)

    floor = sat["c3_baseline_floor_tpm"]
    fold_baseline = _fold(baseline_pancreas)
    fold_cohort = _fold(cohort_normal)
    below_detection = baseline_pancreas < floor

    score = (
        _clamp(math.log2(fold_baseline) / math.log2(sat["c3_fold"]))
        if fold_baseline > 0
        else 0.0
    )
    component = Component(
        score,
        fold_baseline,
        "normal below detection" if below_detection else "",
    )
    return component, fold_baseline, fold_cohort


def _score_c4(values: np.ndarray | None, threshold: float) -> Component:
    """The fraction of malignant cells at or above the detection threshold."""
    if values is None or values.size == 0:
        return Component(None, None, "no row")
    frac = float(np.mean(values >= threshold))
    return Component(_clamp(frac), frac)


def _score_c5(
    residues: int | None, at_plasma_membrane: bool, lipid_anchored: bool, sat: float
) -> Component:
    """Accessibility."""

    if lipid_anchored and residues is None:
        return Component(1.0 if at_plasma_membrane else 0.6, None, "anchor class")
    if residues is None:
        if at_plasma_membrane:
            return Component(None, None, "ectodomain unannotated")
        return Component(None, None, "ectodomain unannotated, no localisation")
    size_term = _clamp(residues / sat)
    confirm = 1.0 if at_plasma_membrane else 0.7
    return Component(_clamp(size_term * confirm), float(residues))


def _score_c6(effect: float | None, screened: int, sat: float) -> Component:
    """Dependency of the tumour lineage on the gene, unmeasured if unscreened."""
    if effect is None or math.isnan(effect) or screened == 0:
        return Component(None, None, "unscreened")
    return Component(_clamp(-effect / sat), effect)


def _baseline_score(tpm: float) -> float:
    """Put a normal-tissue level on the same saturating scale as the rest."""
    return _clamp(math.log10(1 + tpm) / math.log10(1 + BASELINE_TPM_SATURATION))


@dataclass
class RiskModel:
    overrides: dict[str, int] = field(default_factory=dict)
    fall_through: set[str] = field(default_factory=set)

    def tier(self, organ: str) -> int:
        """The organ's consequence tier, honouring any override."""
        return self.overrides.get(organ, ORGAN_TIERS[organ])

    def weight(self, organ: str) -> float:
        """The weight the organ's tier carries."""
        return TIER_WEIGHTS[self.tier(organ)]

    def organ_for_baseline(self, label: str) -> str | None:
        """Map a normal-tissue label to an organ, recording labels that fall through."""
        if label in EXCLUDED_LABELS:
            return None
        organ = BASELINE_ORGANS.get(label)
        if organ is None:
            self.fall_through.add(f"baseline:{label}")
        return organ

    def organ_for_atlas(self, label: str) -> str | None:
        """Map an atlas label to an organ, recording labels that fall through."""
        if label in EXCLUDED_LABELS:
            return None
        organ = ATLAS_ORGANS.get(label)
        if organ is None:
            self.fall_through.add(f"atlas:{label}")
        return organ


def per_organ_scores(
    model: RiskModel,
    atlas_gene,
    baseline_values: np.ndarray | None,
    baseline_tissues: list[str],
    calibration: CalibrationCurve,
) -> dict[str, float]:
    """Expression score per organ, before criticality is applied."""
    per_organ: dict[str, float] = {}

    if atlas_gene is not None:
        for tissue, _cell_type, level in atlas_gene.staining:
            organ = model.organ_for_atlas(tissue)
            if organ is None:
                continue

            score = calibration.score(level)
            if score > per_organ.get(organ, -1.0):
                per_organ[organ] = score

    if baseline_values is not None:
        for label, tpm in zip(baseline_tissues, baseline_values):
            organ = model.organ_for_baseline(label)
            if organ is None:
                continue
            score = _baseline_score(float(tpm))

            if score > per_organ.get(organ, -1.0):
                per_organ[organ] = score

    return per_organ


def protein_arm_measured(model, atlas_gene) -> bool:
    """Whether any mapped normal tissue was stained for this entry."""
    if atlas_gene is None:
        return False
    for tissue, _cell_type, _level in atlas_gene.staining:
        if model.organ_for_atlas(tissue) is not None:
            return True
    return False


NOT_MEASURED = "NOT_MEASURED"

ARM_STAINING = "STAINING"
ARM_BASELINE = "BASELINE"
ARM_TIED = "TIED"


@dataclass(frozen=True)
class ArmReading:
    """One arm's winning measurement within one organ."""

    label: str
    score: float
    level: int | None = None
    level_name: str | None = None
    tpm: float | None = None

    def as_payload(self) -> dict:
        """The reading at full precision, so the score can be recomputed from it."""
        out = {"label": self.label, "score": self.score}
        if self.level is not None:
            out["level"] = self.level
            out["level_name"] = self.level_name
        if self.tpm is not None:
            out["tpm"] = self.tpm
        return out


@dataclass(frozen=True)
class OrganRisk:
    """One organ's inputs to the outer maximum."""

    organ: str
    score: float
    tier: int
    weight: float
    weighted: float
    arm: str
    staining: ArmReading | None
    baseline: ArmReading | None

    def as_payload(self) -> dict:
        """The organ row, with an unmeasured arm named rather than zeroed."""
        return {
            "organ": self.organ,
            "score": self.score,
            "tier": self.tier,
            "weight": self.weight,
            "weighted": self.weighted,
            "arm": self.arm,
            "staining": (NOT_MEASURED if self.staining is None
                         else self.staining.as_payload()),
            "baseline": (NOT_MEASURED if self.baseline is None
                         else self.baseline.as_payload()),
        }


@dataclass(frozen=True)
class RiskAttribution:
    """Every measurement behind a target's risk, and what each contributed."""

    organs: list[OrganRisk]
    risk: float | None
    winners: list[str]
    margin: float | None

    def as_payload(self) -> dict:
        """The attribution as the evidence trail serves it."""
        return {
            "risk": self.risk,
            "winning_organs": self.winners,
            "margin": self.margin,
            "organs_scored": len(self.organs),
            "organs": [o.as_payload() for o in self.organs],
        }


def attribute_risk(
    model: RiskModel,
    atlas_gene,
    baseline_values: np.ndarray | None,
    baseline_tissues: list[str],
    calibration: CalibrationCurve,
) -> RiskAttribution:
    """The inputs to all three reductions, per organ, reconstructing the risk."""
    staining: dict[str, ArmReading] = {}
    baseline: dict[str, ArmReading] = {}

    if atlas_gene is not None:
        for tissue, _cell_type, level in atlas_gene.staining:
            organ = model.organ_for_atlas(tissue)
            if organ is None:
                continue
            score = calibration.score(level)
            best = staining.get(organ)
            if best is None or score > best.score:
                staining[organ] = ArmReading(
                    label=tissue, score=score, level=level,
                    level_name=LEVEL_NAMES.get(level, str(level)))

    if baseline_values is not None:
        for label, tpm in zip(baseline_tissues, baseline_values):
            organ = model.organ_for_baseline(label)
            if organ is None:
                continue
            value = float(tpm)
            score = _baseline_score(value)
            best = baseline.get(organ)
            if best is None or score > best.score:
                baseline[organ] = ArmReading(label=label, score=score, tpm=value)

    organs: list[OrganRisk] = []
    for organ in sorted(set(staining) | set(baseline)):
        stained, based = staining.get(organ), baseline.get(organ)
        if stained is not None and based is not None:
            score = max(stained.score, based.score)
            arm = (ARM_TIED if stained.score == based.score
                   else ARM_STAINING if stained.score > based.score
                   else ARM_BASELINE)
        elif stained is not None:
            score, arm = stained.score, ARM_STAINING
        else:
            score, arm = based.score, ARM_BASELINE
        weight = model.weight(organ)
        organs.append(OrganRisk(
            organ=organ, score=score, tier=model.tier(organ), weight=weight,
            weighted=score * weight, arm=arm,
            staining=stained, baseline=based))

    organs.sort(key=lambda o: (-o.weighted, o.organ))
    if not organs:
        return RiskAttribution([], None, [], None)

    top = organs[0].weighted
    winners = [o.organ for o in organs if o.weighted == top]
    margin = top - organs[1].weighted if len(organs) > 1 else None
    return RiskAttribution(organs, top, winners, margin)


def worst_organ(
    model: RiskModel, per_organ: dict[str, float]
) -> tuple[float | None, str | None]:
    """Criticality-weighted maximum over organs."""
    if not per_organ:
        return None, None
    organ = max(per_organ, key=lambda o: per_organ[o] * model.weight(o))
    return per_organ[organ] * model.weight(organ), organ


@dataclass(frozen=True)
class RiskInputs:
    """What attribution needs to reconstruct any target's risk on demand."""

    model: RiskModel
    calibration: CalibrationCurve
    by_accession: dict
    by_symbol: dict
    profiles: dict
    tissues: list[str]

    def attribute(self, accession: str, gene: str) -> RiskAttribution:
        """Reconstruct one target's risk from the measurements underneath it."""
        profile = self.profiles.get(accession)
        return attribute_risk(
            self.model,
            self.by_accession.get(accession) or self.by_symbol.get(gene),
            profile.values if profile is not None else None,
            self.tissues,
            self.calibration,
        )


def compute_risk(
    model: RiskModel,
    atlas_gene,
    baseline_values: np.ndarray | None,
    baseline_tissues: list[str],
    calibration: CalibrationCurve,
) -> tuple[float | None, str | None]:
    """The worst organ's normal-tissue risk, and which organ it was."""
    return worst_organ(
        model,
        per_organ_scores(
            model, atlas_gene, baseline_values, baseline_tissues, calibration
        ),
    )


def _confidence(
    components: dict[str, Component], row: CoverageRow, wts: dict[str, float]
) -> float:
    """Evidence confidence, from measured weight and the annotation tier."""
    measured = sum(wts[k] for k, c in components.items() if c.measured)
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
    calibration: CalibrationCurve,
    saturation: dict[str, float] | None = None,
    margin_label: str | None = None,
    weights: dict[str, float] | None = None,
) -> tuple[list[Ranked], RiskModel, dict]:
    """Score every surface protein and return the ranking."""
    sat = dict(SATURATION if saturation is None else saturation)
    wts = dict(WEIGHTS if weights is None else weights)

    unknown = sorted(set(overrides) - set(ORGAN_TIERS))
    if unknown:
        raise KeyError(
            "criticality override names organs that do not exist: "
            + ", ".join(unknown)
            + f"; known organs are {', '.join(sorted(ORGAN_TIERS))}"
        )
    model = RiskModel(overrides=overrides)

    margin_label = margin_label or DEFAULT_MARGIN_LABEL
    if margin_label not in gtex_tissues:
        raise KeyError(
            f"the baseline has no {margin_label!r} column; "
            "the margin component has no denominator"
        )
    bulk_pancreas_col = gtex_tissues.index(margin_label)

    primary_mask = cohort.sample_types == PRIMARY_TUMOUR
    normal_mask = cohort.sample_types == SOLID_NORMAL

    out: list[Ranked] = []
    for row in coverage_rows:
        rec = surface_by_accession[row.accession]
        gene = row.gene

        cell_hit = cell_index.get(row.accession)
        malignant = stromal_peak = None
        cell_bridged = False
        if cell_hit is not None:
            ci, cell_route = cell_hit
            cell_bridged = cell_route == CELL_JOIN_ENSEMBL_BRIDGE
            malignant = cell_atlas.compartment_value(MALIGNANT, ci)
            peaks = [cell_atlas.compartment_value(c, ci) for c in STROMAL_COMPARTMENTS]
            peaks = [p for p in peaks if not math.isnan(p)]
            stromal_peak = max(peaks) if peaks else None

        tumour_median = cohort_normal_median = None
        prevalence_values = None
        cj = cohort_join.get(row.accession)
        if cj is not None:
            col = cohort.values[:, cj[0]]
            candidate = col[primary_mask]

            if candidate.size:
                prevalence_values = candidate
                tumour_median = _finite(float(np.median(candidate)))
            if normal_mask.any():
                cohort_normal_median = _finite(float(np.median(col[normal_mask])))

        profile = gtex_profiles.get(row.accession)
        baseline_pancreas = None
        if profile is not None:
            baseline_pancreas = _finite(float(profile.values[bulk_pancreas_col]))

        c3, fold_b, fold_c = _score_c3(
            tumour_median, baseline_pancreas, cohort_normal_median, sat
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
            calibration,
        )

        stroma = components[C2]
        if not stroma.measured:
            tumour_side = STROMA_UNRESOLVED
        elif (stroma.raw or 0.0) <= STROMA_RATIO_FLOOR:
            tumour_side = STROMA_DOMINANT
        else:
            tumour_side = TUMOUR_DOMINANT

        has_protein = protein_arm_measured(model, atlas_gene)
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
                tumour_side_verdict=tumour_side,
                protein_arm_measured=has_protein,
                risk_basis=("staining and transcript" if has_protein
                            else "transcript only"),
                risk_is_lower_bound=not has_protein,
                composite_supported=(
                    None if composite is None
                    else round(composite * measured_weight, 4)),
                confidence=_confidence(components, row, wts),
                sources_disagree=False,
                fold_baseline=fold_b,
                fold_cohort=fold_c,
                bridged=row.bridged or cell_bridged,
            )
        )

    log_ratios = [
        math.log10(r.fold_cohort / r.fold_baseline)
        for r in out
        if r.fold_baseline and r.fold_cohort and r.fold_baseline > 0 and r.fold_cohort > 0
    ]
    offset = float(np.median(log_ratios)) if log_ratios else 0.0
    tolerance = math.log10(FOLD_DISAGREEMENT)

    for r in out:
        if not (r.fold_baseline and r.fold_cohort):
            continue
        if r.fold_baseline <= 0 or r.fold_cohort <= 0:
            continue

        r.field_elevation = r.fold_baseline / r.fold_cohort
        residual = math.log10(r.fold_cohort / r.fold_baseline) - offset
        r.sources_disagree = abs(residual) > tolerance

    stats = {
        "field_offset_log10": offset,
        "field_offset_fold": 10.0**-offset,
        "proteins_with_both_folds": len(log_ratios),
        "tolerance_fold": FOLD_DISAGREEMENT,
    }

    for tier in (PROTEIN_CONFIRMED, RNA_SUPPORTED, DATA_INSUFFICIENT):
        members = [r for r in out if r.evidence_class == tier]
        members.sort(
            key=lambda r: (r.composite_supported is None,
                           -(r.composite_supported or 0.0), r.gene or "")
        )
        for i, r in enumerate(members, start=1):
            r.tier_rank = i

    return out, model, stats


def _revision() -> str:
    """The current revision, marked dirty when the tree has uncommitted changes."""
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
    """The repository root."""
    from pathlib import Path

    return str(Path(__file__).resolve().parents[2])


def configuration_hash(
    overrides: dict[str, int],
    ceiling: float,
    saturation: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
    calibration: CalibrationCurve | None = None,
    margin_label: str | None = None,
) -> str:
    """Covers the tissue tables themselves, not only the tier assignments."""
    payload = {
        "weight_set_version": WEIGHT_SET_VERSION,
        "weights": WEIGHTS if weights is None else weights,
        "minimum_measured_weight": MINIMUM_MEASURED_WEIGHT,
        "saturation": SATURATION if saturation is None else saturation,
        "fold_disagreement": FOLD_DISAGREEMENT,
        "dropout_epsilon": DROPOUT_EPSILON,
        "margin_denominator": margin_label or DEFAULT_MARGIN_LABEL,
        "tier_weights": TIER_WEIGHTS,
        "organ_tiers": ORGAN_TIERS,
        "baseline_organs": BASELINE_ORGANS,
        "atlas_organs": ATLAS_ORGANS,
        "excluded_labels": sorted(EXCLUDED_LABELS),
        "atlas_level_calibration": (
            calibration.as_payload() if calibration is not None else None
        ),
        "baseline_tpm_saturation": BASELINE_TPM_SATURATION,
        "overrides": overrides,
        "ceiling": ceiling,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def header(
    spec,
    model: RiskModel,
    universe: int,
    pins: dict,
    saturation: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
    stats: dict | None = None,
    calibration: CalibrationCurve | None = None,
    margin_label: str | None = None,
) -> str:
    """Describes the run that produced the output, not the module defaults."""
    sat = SATURATION if saturation is None else saturation
    wts = WEIGHTS if weights is None else weights
    margin_label = margin_label or DEFAULT_MARGIN_LABEL
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
        f"  weights              " + ", ".join(f"{k}={v}" for k, v in wts.items()),
        f"  minimum measured     {MINIMUM_MEASURED_WEIGHT}",
        f"  free parameters      " + ", ".join(f"{k}={v}" for k, v in sat.items()),
        f"  dropout epsilon      {DROPOUT_EPSILON}",
        f"  fold disagreement    {FOLD_DISAGREEMENT}x",
        f"  risk ceiling         {ceiling}",
        f"  universe             {universe}",
        f"  revision             {_revision()}",
        f"  configuration hash   "
        f"{configuration_hash(model.overrides, ceiling, sat, wts, calibration, margin_label)}",
        f"  margin denominator   {margin_label} (bulk only)",
        f"  field offset         "
        + (
            f"the cohort fold runs {stats['field_offset_fold']:.1f}x below the "
            f"baseline fold at the median, over "
            f"{stats['proteins_with_both_folds']:,} proteins with both"
            if stats
            else "not measured"
        ),
        f"                       the two normal tissues differ by about this "
        f"much by construction;",
        f"                       {FIELD_ELEVATION} is each protein's own "
        f"version of that distance,",
        f"                       and the disagreement flag marks departures "
        f"from the offset, not from parity",
        "  criticality tiers    [+] marks a platform addition, not from the "
        "reference table",
    ]
    for tier in (1, 2, 3):
        organs = sorted(o for o, t in ORGAN_TIERS.items() if t == tier)
        marked = [
            f"{o}[+]" if o in PLATFORM_ADDED_ORGANS else o for o in organs
        ]
        lines.append(f"    tier {tier} (w={TIER_WEIGHTS[tier]}): {', '.join(marked)}")
    if model.overrides:
        lines.append("  overrides")
        for organ, tier in model.overrides.items():
            lines.append(f"    {organ} -> tier {tier} (w={TIER_WEIGHTS[tier]})")
            rationale = spec.inputs.tissue_criticality_overrides.get(organ)
            if rationale is not None:
                lines.append(f"      rationale: {rationale.rationale}")
    if calibration is not None:
        lines.append(
            "  staining calibration measured against the transcript axis, "
            f"{calibration.observations:,} paired organs"
        )
        lines.append(
            f"    {'level':14s} {'n':>7s} {'Q1':>9s} {'median':>9s} "
            f"{'Q3':>9s} {'score':>8s}"
        )
        for k in sorted(calibration.tpm):
            q1, med, q3 = calibration.quartiles[k]
            lines.append(
                f"    {LEVEL_NAMES[k]:14s} {calibration.counts[k]:>7,} "
                f"{q1:>9.3f} {med:>9.3f} {q3:>9.3f} "
                f"{calibration.score(k):>8.4f}"
            )
        sep = "  ".join(
            f"{LEVEL_NAMES[k][:4]}->{LEVEL_NAMES[k + 1][:4]} "
            f"{calibration.separations[k]:.3f}"
            for k in sorted(calibration.separations)
        )
        lines.append(f"    monotonic: {'yes' if calibration.monotonic else 'NO'}")
        lines.append(f"    separation (0.50 = no information): {sep}")
        lines.append(
            "    the scale is real but weak; a one-level difference is not decisive"
        )
    lines.append("  source pins")
    for name, pin in pins.items():
        lines.append(f"    {name}: {pin}")
    lines.append("=" * 72)
    return "\n".join(lines)
