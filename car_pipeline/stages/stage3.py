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

#: Free parameters that set where each component stops rewarding more. Named
#: here so they can be perturbed as a set by the twelfth rejection criterion.
SATURATION: dict[str, float] = {
    "c1_expression": 100.0,
    "c2_ratio": 50.0,
    "c3_fold": 64.0,
    "c4_prevalence_tpm": 10.0,
    "c5_ectodomain_residues": 200.0,
    "c6_gene_effect": 1.0,
    # Below this, a normal tissue reading cannot be told apart from zero. It
    # floors the margin denominator instead of the ratio being taken as
    # unbounded. Listed here with the other free parameters so the twelfth
    # criterion perturbs it too: it moves the top of the ranking as hard as any
    # saturation point does.
    "c3_baseline_floor_tpm": 0.1,
}

#: How far a protein must depart from the *systematic* offset between the two
#: denominators before it is flagged. Not deviation from parity: the two normal
#: tissues differ by construction, so a parity test flags everything and
#: discriminates nothing.
FOLD_DISAGREEMENT = 2.0

#: Ratio of the two normal tissues for one protein: how much more of it sits in
#: tumour-adjacent pancreas than in healthy population pancreas.
#:
#: The population baseline is pristine acinar tissue. The cohort's adjacent
#: normal comes from a resected pancreas, which is typically inflamed and often
#: already carries precursor lesions, so it expresses tumour antigens before any
#: malignancy is involved. The distance between the two therefore measures how
#: much of a target's apparent margin comes from the diseased field rather than
#: from the tumour.
#:
#: That distinction is clinical, not academic: what remains in the patient after
#: resection is field tissue, so an antigen elevated here is an antigen the
#: therapy will meet in tissue that is staying put. Reported on every row, never
#: scored — it describes where a margin comes from, not how large it is.
#:
#: Named for the magnitude rather than the association, since it is a measured
#: ratio and the direction matters.
FIELD_ELEVATION = "field_elevation"

STROMAL_COMPARTMENTS = (FIBROBLAST, IMMUNE, ENDOTHELIAL)

# --------------------------------------------------------------------------
# tissue criticality
# --------------------------------------------------------------------------

TIER_WEIGHTS = {1: 1.0, 2: 0.6, 3: 0.3}

ORGAN_TIERS: dict[str, int] = {
    # tier 1
    "brain": 1, "heart": 1, "lung": 1, "liver": 1, "kidney": 1,
    "pancreas": 1, "vascular": 1, "eye": 1,
    # tier 2
    "gi_tract": 2, "marrow_and_blood": 2, "bladder": 2, "endocrine": 2,
    "muscle": 2, "nerve": 2, "mucosa": 2,
    # tier 3
    "skin": 3, "adipose": 3, "breast": 3, "reproductive": 3, "salivary": 3,
    "connective": 3,
}

#: Organs not present in the reference criticality table, assigned here. Marked
#: in the output so a reader can tell which tiers were inherited and which were
#: supplied by the platform, rather than having to trust the whole table equally.
PLATFORM_ADDED_ORGANS = frozenset({"vascular", "eye", "mucosa", "connective"})

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

LEVEL_NAMES = {0: "Not detected", 1: "Low", 2: "Medium", 3: "High"}


@dataclass
class CalibrationCurve:
    """Staining levels placed on the transcript axis by measurement.

    The two sources describe the same organs in different units, so until they
    sit on one axis the worst-organ maximum is comparing incomparable numbers.
    Rather than asserting an evenly spaced ordinal, every organ carrying both a
    staining call and a transcript value contributes one observation, and each
    level is represented by the median of its own population. That value is then
    scored through the same continuous function the transcript side uses, so the
    two axes are commensurable by construction rather than by assertion.
    """

    tpm: dict[int, float]
    quartiles: dict[int, tuple[float, float, float]]
    counts: dict[int, int]
    separations: dict[int, float]
    monotonic: bool
    observations: int

    def score(self, level: int) -> float:
        # A non-detection scores zero, not the transcript level its population
        # happens to sit at. The calibrated value for that level is a central
        # estimate of a wide distribution, and using it as a risk term means an
        # observation of *absence* raises measured risk — which it cannot.
        # Scoring it zero still leaves the organ measured, so the maximum below
        # lets a positive transcript reading stand rather than being cancelled.
        if level == 0:
            return 0.0
        return _baseline_score(self.tpm[level])

    def as_payload(self) -> dict:
        return {str(k): round(v, 6) for k, v in sorted(self.tpm.items())}


def _separation(lower: np.ndarray, upper: np.ndarray) -> float:
    """Probability that a draw from the upper level exceeds one from the lower.

    0.50 means the two levels say nothing about each other; 1.00 means they
    separate perfectly. Ties count as half.
    """
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
    """Measure what each staining level is worth in transcript terms.

    Pairs are formed exactly the way risk is computed: the level for an organ is
    the maximum across its cell types, and the transcript value is the maximum
    across the tissues mapping to that organ. Calibrating the same quantity that
    gets scored is what lets the curve absorb whatever inflation that
    aggregation introduces, instead of leaving it to be argued about.
    """
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

#: The margin component's denominator, and only its denominator. The baseline
#: names four pancreas entries; three are cell-sorted fractions and one is bulk.
#: The tumour side of the comparison is bulk, so a median across all four mixes
#: measurement types on one side of a ratio. The risk gate deliberately keeps
#: reading all four and taking the worst of them: for safety the question is
#: whether any pancreatic compartment carries the antigen, not what the organ
#: averages.
BULK_PANCREAS_LABEL = "Pancreas"


def _clamp(x: float) -> float:
    return 0.0 if x < 0 else (1.0 if x > 1 else x)


def _finite(x: float | None) -> float | None:
    """Pass a real number through; turn anything else into an absent one."""
    if x is None or math.isnan(x) or math.isinf(x):
        return None
    return x


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
    #: How much more of this antigen sits in tumour-adjacent pancreas than in
    #: healthy pancreas. An annotation, never a score term.
    field_elevation: float | None = None
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
    sat: dict[str, float],
) -> tuple[Component, float | None, float | None]:
    if tumour is None or baseline_pancreas is None:
        # No row on one side. Genuinely unmeasured.
        return Component(None, None, "no row"), None, None

    def _fold(normal: float | None) -> float | None:
        """Tumour over normal, with the denominator floored at detection.

        A normal reading of zero is a measurement, not a gap — the protein was
        looked for and not found, which is the most favourable thing this
        comparison can say. But it cannot be treated as an unbounded ratio
        either: dividing by it awards a perfect margin to proteins absent from
        the tumour as well, where there is no margin because there is nothing
        there. Flooring the denominator at the detection limit keeps the
        favourable reading favourable while leaving the numerator to decide
        whether anything is actually present.
        """
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
    # The anchor special case exists because that class has no transmembrane
    # segment to annotate topological domains around, so it reports zero
    # residues by construction. That reasoning only holds while there is no
    # annotation. A handful of proteins carry both an anchor and a segment, and
    # for those the ectodomain really was measured — discarding it in favour of
    # a class default would throw away the better evidence of the two.
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
    calibration: CalibrationCurve,
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
            # Converted to its measured transcript equivalent and scored by the
            # same function the transcript side uses, so the maximum below is
            # comparing like with like.
            score = calibration.score(level)
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


def _confidence(
    components: dict[str, Component], row: CoverageRow, wts: dict[str, float]
) -> float:
    # Uses the weights this run actually scored with, not the module defaults,
    # so a perturbed run cannot report a confidence describing a different
    # weight set from its own composites.
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
    weights: dict[str, float] | None = None,
) -> tuple[list[Ranked], RiskModel, dict]:
    sat = dict(SATURATION if saturation is None else saturation)
    wts = dict(WEIGHTS if weights is None else weights)

    # An override naming an organ that does not exist would do nothing at all,
    # while the output header went on reporting it as an applied relaxation
    # complete with its rationale. A safety default that only appears to have
    # been changed is worse than one that was never touched.
    unknown = sorted(set(overrides) - set(ORGAN_TIERS))
    if unknown:
        raise KeyError(
            "criticality override names organs that do not exist: "
            + ", ".join(unknown)
            + f"; known organs are {', '.join(sorted(ORGAN_TIERS))}"
        )
    model = RiskModel(overrides=overrides)

    if BULK_PANCREAS_LABEL not in gtex_tissues:
        # Failing here is the point. Silently having no denominator would make
        # the margin component unmeasured for every protein at once, and the
        # evidence floor would then quietly drop a third of the universe.
        raise KeyError(
            f"the baseline has no {BULK_PANCREAS_LABEL!r} column; "
            "the margin component has no denominator"
        )
    bulk_pancreas_col = gtex_tissues.index(BULK_PANCREAS_LABEL)

    primary_mask = cohort.sample_types == PRIMARY_TUMOUR
    normal_mask = cohort.sample_types == SOLID_NORMAL

    out: list[Ranked] = []
    for row in coverage_rows:
        rec = surface_by_accession[row.accession]
        gene = row.gene

        # cell atlas
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

        # cohort
        tumour_median = cohort_normal_median = None
        prevalence_values = None
        cj = cohort_join.get(row.accession)
        if cj is not None:
            col = cohort.values[:, cj[0]]
            candidate = col[primary_mask]
            # An empty group has no median, and a not-a-number median is not a
            # measurement of zero. Either would otherwise flow into the margin
            # as a real reading and be scored as a target with no enrichment.
            if candidate.size:
                prevalence_values = candidate
                tumour_median = _finite(float(np.median(candidate)))
            if normal_mask.any():
                cohort_normal_median = _finite(float(np.median(col[normal_mask])))

        # baseline
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
                confidence=_confidence(components, row, wts),
                sources_disagree=False,  # set in the second pass, below
                fold_baseline=fold_b,
                fold_cohort=fold_c,
                # Any source reached only after its symbol failed counts, the
                # cell atlas included. Its join feeds the two heaviest
                # components, so leaving it out of the audit would exempt the
                # most consequential bridge from the criterion that polices them.
                bridged=row.bridged or cell_bridged,
            )
        )

    # -- second pass: the systematic offset, then departures from it ------
    #
    # The two denominators describe different tissue, so they disagree by
    # construction and by a fairly consistent amount. Testing each protein
    # against parity therefore flags every one of them and identifies none. The
    # informative question is which proteins depart from the pattern, so the
    # pattern is measured first and each protein tested against that.
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
        # How far the two normal tissues sit apart for this protein.
        r.field_elevation = r.fold_baseline / r.fold_cohort
        residual = math.log10(r.fold_cohort / r.fold_baseline) - offset
        r.sources_disagree = abs(residual) > tolerance

    stats = {
        "field_offset_log10": offset,
        "field_offset_fold": 10.0**-offset,
        "proteins_with_both_folds": len(log_ratios),
        "tolerance_fold": FOLD_DISAGREEMENT,
    }

    # rank within evidence class, scored rows first
    for tier in (PROTEIN_CONFIRMED, RNA_SUPPORTED, DATA_INSUFFICIENT):
        members = [r for r in out if r.evidence_class == tier]
        members.sort(
            key=lambda r: (r.composite is None, -(r.composite or 0.0), r.gene or "")
        )
        for i, r in enumerate(members, start=1):
            r.tier_rank = i

    return out, model, stats


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


def configuration_hash(
    overrides: dict[str, int],
    ceiling: float,
    saturation: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
    calibration: CalibrationCurve | None = None,
) -> str:
    """Covers the tissue tables themselves, not only the tier assignments.

    Adding or moving a single label shifts thousands of risk scores. A hash that
    did not move with it would let two materially different experiments compare
    as the same one, which is worse than publishing no hash at all.
    """
    payload = {
        "weight_set_version": WEIGHT_SET_VERSION,
        "weights": WEIGHTS if weights is None else weights,
        "minimum_measured_weight": MINIMUM_MEASURED_WEIGHT,
        "saturation": SATURATION if saturation is None else saturation,
        "fold_disagreement": FOLD_DISAGREEMENT,
        "dropout_epsilon": DROPOUT_EPSILON,
        # The column the margin denominator is drawn from. Swapping it moves
        # every fold change and therefore the whole ranking, so a hash that
        # ignored it would let two different experiments compare as one.
        "margin_denominator": BULK_PANCREAS_LABEL,
        "tier_weights": TIER_WEIGHTS,
        "organ_tiers": ORGAN_TIERS,
        "baseline_organs": BASELINE_ORGANS,
        "atlas_organs": ATLAS_ORGANS,
        "excluded_labels": sorted(EXCLUDED_LABELS),
        # The measured curve, not an assumed scale. It is derived from the
        # pinned releases, so it belongs to the experiment: two runs that
        # calibrated differently must not compare as the same one.
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
) -> str:
    """Describes the run that produced the output, not the module defaults."""
    sat = SATURATION if saturation is None else saturation
    wts = WEIGHTS if weights is None else weights
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
        f"{configuration_hash(model.overrides, ceiling, sat, wts, calibration)}",
        f"  margin denominator   {BULK_PANCREAS_LABEL} (bulk only)",
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
