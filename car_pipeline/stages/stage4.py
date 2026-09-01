"""Stage 4 — target pairing and the single-versus-dual decision."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

from car_pipeline.stages import stage3
from car_pipeline.stages.stage3 import RiskModel, Ranked


from car_pipeline.stages import routing

POOL_SIZE = 200


DETECTION_COUNTS = 1
SENSITIVITY_COUNTS = (1, 2, 3)


COVERAGE_FLOOR = 0.02
PATIENT_FRACTION_FLOOR = 0.60


MIN_MALIGNANT_CELLS = 100


MIN_DETECTED_CELLS = 10

SINGLE = "SINGLE"
DUAL = "DUAL"
NO_DESIGN = "NO_DESIGN"


ADAPTOR = "ADAPTOR"
UNRESOLVED = "UNRESOLVED"


@dataclass
class PairRisk:
    combined: float | None
    optimistic: float | None
    independence: float | None

    pessimistic: float | None
    organ: str | None
    unresolved_organs: list[str] = field(default_factory=list)

    @property
    def risk_unresolved(self) -> bool:
        """Whether any organ was measured for only one member."""
        return bool(self.unresolved_organs)


def pair_risk(
    model: RiskModel, a: dict[str, float], b: dict[str, float]
) -> PairRisk:
    """Minimum per organ, then maximum across organs."""
    conservative: dict[str, float] = {}
    optimistic: dict[str, float] = {}
    independence: dict[str, float] = {}
    pessimistic: dict[str, float] = {}
    unresolved: list[str] = []

    for organ in set(a) | set(b):
        sa, sb = a.get(organ), b.get(organ)
        if sa is not None and sb is not None:
            conservative[organ] = optimistic[organ] = min(sa, sb)
            independence[organ] = sa * sb
            pessimistic[organ] = min(sa, sb)
            continue

        unresolved.append(organ)
        known = sa if sa is not None else sb
        conservative[organ] = known
        independence[organ] = known
        pessimistic[organ] = 1.0

    combined, organ = stage3.worst_organ(model, conservative)
    opt, _ = stage3.worst_organ(model, optimistic)
    ind, _ = stage3.worst_organ(model, independence)
    pess, _ = stage3.worst_organ(model, pessimistic)

    return PairRisk(
        _round(combined),
        _round(opt),
        _round(ind),
        _round(pess),
        organ,
        sorted(unresolved),
    )


def _round(value: float | None) -> float | None:
    """Round to the precision risk is stored at, passing None through."""
    return None if value is None else round(value, 4)


@dataclass
class Coverage:
    measured: bool
    f_a: float | None = None
    f_b: float | None = None
    f_ab: float | None = None
    p_a_given_b: float | None = None
    p_b_given_a: float | None = None
    sacrificed_a: float | None = None
    sacrificed_b: float | None = None
    jaccard: float | None = None
    patients_at_floor: int = 0
    patients_evaluable: int = 0
    reason: str = ""

    span_geomean_kb: float | None = None
    span_percentile: float | None = None

    @property
    def patient_fraction(self) -> float:
        """The share of evaluable patients reaching the coverage floor."""
        if not self.patients_evaluable:
            return 0.0
        return self.patients_at_floor / self.patients_evaluable

    @property
    def escape(self) -> float | None:
        """Malignant cells the AND gate does not engage."""
        return None if self.f_ab is None else 1.0 - self.f_ab

    @property
    def or_gate(self) -> float | None:
        """Union: what a design firing on either antigen would reach."""
        if self.f_ab is None:
            return None
        return self.f_a + self.f_b - self.f_ab

    @property
    def best_single(self) -> float | None:
        """The better of the two single-antigen coverages."""
        return None if self.f_ab is None else max(self.f_a, self.f_b)

    @property
    def coverage_cost(self) -> float | None:
        """How many times more of the tumour the better single target reaches."""
        if not self.f_ab:
            return None
        return self.best_single / self.f_ab


def intersection_matrix(positive: np.ndarray) -> np.ndarray:
    """Pairwise double-positive counts for every gene pair at once."""
    m = positive.astype(np.float32)
    return (m.T @ m).astype(np.int64)


def coverage_from_counts(
    n_cells: int,
    n_a: int,
    n_b: int,
    n_ab: int,
    patients_at_floor: int,
    patients_evaluable: int,
) -> Coverage:
    """Build the coverage figures from the per-cell and per-patient counts."""
    if n_a < MIN_DETECTED_CELLS or n_b < MIN_DETECTED_CELLS:
        return Coverage(measured=False, reason="below detection")

    f_a, f_b = n_a / n_cells, n_b / n_cells
    f_ab = n_ab / n_cells
    union = n_a + n_b - n_ab
    return Coverage(
        measured=True,
        f_a=f_a,
        f_b=f_b,
        f_ab=f_ab,
        p_a_given_b=n_ab / n_b,
        p_b_given_a=n_ab / n_a,
        sacrificed_a=1.0 - n_ab / n_a,
        sacrificed_b=1.0 - n_ab / n_b,
        jaccard=(n_ab / union) if union else 0.0,
        patients_at_floor=patients_at_floor,
        patients_evaluable=patients_evaluable,
    )


@dataclass
class Pair:
    gene_a: str
    gene_b: str
    accession_a: str
    accession_b: str
    risk: PairRisk
    coverage: Coverage
    risk_a: float
    risk_b: float
    cleared_a: bool
    cleared_b: bool
    organ_a: str | None
    organ_b: str | None
    composite_a: float
    composite_b: float
    confidence_a: float
    confidence_b: float
    ceiling: float
    organs_total: int = 0
    organs_resolved: int = 0

    @property
    def confidence(self) -> float:
        """Third number, never combined with the other two."""
        base = min(self.confidence_a, self.confidence_b)
        resolved = (
            self.organs_resolved / self.organs_total if self.organs_total else 0.0
        )
        measured = 1.0 if self.coverage.measured else 0.75
        return round(base * resolved * measured, 4)

    @property
    def cleared(self) -> bool:
        """Conservative value, so an unresolved organ cannot clear a pair."""
        return self.risk.combined is not None and self.risk.combined <= self.ceiling

    @property
    def delta_a(self) -> float | None:
        """How much risk the pairing removes from the first member."""
        if self.risk.combined is None:
            return None
        return self.risk_a - self.risk.combined

    @property
    def delta_b(self) -> float | None:
        """How much risk the pairing removes from the second member."""
        if self.risk.combined is None:
            return None
        return self.risk_b - self.risk.combined

    @property
    def rescued(self) -> list[str]:
        """A condition, not a statistic."""
        if not self.cleared:
            return []
        out = []
        if self.risk_a > self.ceiling:
            out.append(self.gene_a)
        if self.risk_b > self.ceiling:
            out.append(self.gene_b)
        return out

    @property
    def coverage_ok(self) -> bool:
        """Whether measured coverage clears both the pool and patient floors."""
        c = self.coverage
        return (
            c.measured
            and c.f_ab >= COVERAGE_FLOOR
            and c.patient_fraction >= PATIENT_FRACTION_FLOOR
        )

    @property
    def admissible(self) -> bool:
        """Risk-gated and measurable. Coverage does not gate."""
        return self.cleared and self.coverage.measured


def build_pool(rows: list[Ranked], size: int = POOL_SIZE) -> list[Ranked]:
    """Top of the tumour-side ranking by supported score, risk ignored entirely."""
    scored = [r for r in rows
              if r.composite_supported is not None and r.gene
              and r.tumour_side_verdict != stage3.STROMA_DOMINANT]
    scored.sort(key=lambda r: (-r.composite_supported, r.gene, r.accession))
    seen: set[str] = set()
    pool: list[Ranked] = []
    for r in scored:
        if r.gene in seen:
            continue
        seen.add(r.gene)
        pool.append(r)
        if len(pool) == size:
            break
    return pool


def evaluate(
    pool: list[Ranked],
    per_organ: dict[str, dict[str, float]],
    model: RiskModel,
    ceiling: float,
    cells,
    threshold: int = DETECTION_COUNTS,
) -> list[Pair]:
    """Every pair in the pool. No pair is excluded from computation."""
    missing_risk = [r.gene for r in pool if r.risk is None]
    if missing_risk:
        raise ValueError(
            "pool members carry no normal tissue risk: " + ", ".join(missing_risk)
        )

    column = {g: i for i, g in enumerate(cells.genes)}
    present: list[int] = []
    slot: dict[int, int] = {}
    for pos, r in enumerate(pool):
        i = column.get(r.gene)
        if i is None:
            continue
        slot[pos] = len(present)
        present.append(i)

    positive = cells.counts[:, np.asarray(present)] >= threshold
    n_cells = positive.shape[0]
    inter = intersection_matrix(positive)

    at_floor = np.zeros_like(inter)
    labels, counts = np.unique(cells.patient, return_counts=True)
    evaluable = [l for l, c in zip(labels, counts) if c >= MIN_MALIGNANT_CELLS]
    for label in evaluable:
        rows_p = cells.patient == label
        sub = positive[rows_p]
        at_floor += (intersection_matrix(sub) >= COVERAGE_FLOOR * sub.shape[0])
    n_evaluable = len(evaluable)

    out: list[Pair] = []
    for i in range(len(pool)):
        for j in range(i + 1, len(pool)):
            ra, rb = pool[i], pool[j]
            si, sj = slot.get(i), slot.get(j)
            if si is None or sj is None:
                cov = Coverage(measured=False, reason="no row in the cell atlas")
            else:
                cov = coverage_from_counts(
                    n_cells,
                    int(inter[si, si]),
                    int(inter[sj, sj]),
                    int(inter[si, sj]),
                    int(at_floor[si, sj]),
                    n_evaluable,
                )
            pr = pair_risk(
                model, per_organ[ra.accession], per_organ[rb.accession]
            )
            union = len(set(per_organ[ra.accession]) | set(per_organ[rb.accession]))
            out.append(
                Pair(
                    gene_a=ra.gene,
                    gene_b=rb.gene,
                    accession_a=ra.accession,
                    accession_b=rb.accession,
                    risk=pr,
                    coverage=cov,
                    confidence_a=ra.confidence,
                    confidence_b=rb.confidence,
                    organs_total=union,
                    organs_resolved=union - len(pr.unresolved_organs),
                    risk_a=ra.risk,
                    risk_b=rb.risk,
                    cleared_a=ra.cleared,
                    cleared_b=rb.cleared,
                    organ_a=ra.risk_organ,
                    organ_b=rb.risk_organ,
                    composite_a=ra.composite,
                    composite_b=rb.composite,
                    ceiling=ceiling,
                )
            )
    return out


@dataclass
class Decision:
    gene: str
    outcome: str
    partner: str | None = None
    pair: Pair | None = None
    failed_on: dict[str, int] = field(default_factory=dict)

    accession: str = ""
    partner_accession: str | None = None

    pool_index: int = -1
    partner_options: int | None = None
    partner_forced: bool | None = None

    architecture: str = ""
    route_reason: str = ""
    route_ceiling: float | None = None
    route_exposure: str | None = None

    evidence_class: str = ""
    tumour_side_verdict: str = ""
    normal_tissue_risk: float | None = None
    normal_tissue_risk_organ: str | None = None
    evidence_confidence: float | None = None
    protein_arm_measured: bool | None = None
    risk_basis: str = ""
    risk_is_lower_bound: bool | None = None
    composite: float | None = None
    composite_supported: float | None = None
    measured_weight: float | None = None
    components_measured: int | None = None
    components_total: int | None = None
    target_confidence: float | None = None
    risk_organs_measured: int | None = None

    partner_evidence_class: str | None = None
    partner_measured_weight: float | None = None
    pair_confidence: float | None = None
    pair_organs_total: int | None = None
    pair_organs_resolved: int | None = None
    pair_coverage_measured: bool | None = None


SPAN_BUCKETS = 5


PARTNER_MIN_TUMOUR_TPM = 5.0


SELECTION_RULE = (
    "risk-cleared-and-measured;partner_min_tumour_tpm;"
    "order:combined_risk,partner_name;v3"
)


def annotate_span_context(pairs: list[Pair], spans: dict[str, int]) -> int:
    """Attach each measured pair's span and its within-span percentile."""
    import numpy as np

    scored = []
    for pair in pairs:
        cov = pair.coverage
        if not cov.measured or cov.f_ab is None:
            continue
        a, b = spans.get(pair.gene_a), spans.get(pair.gene_b)
        if not a or not b:
            continue
        cov.span_geomean_kb = round(float(np.sqrt(float(a) * float(b))) / 1000.0, 3)
        scored.append(pair)

    if not scored:
        return 0

    geo = np.array([p.coverage.span_geomean_kb for p in scored], dtype=float)
    edges = np.percentile(geo, [100 * i / SPAN_BUCKETS for i in range(1, SPAN_BUCKETS)])
    bucket = np.digitize(geo, edges)
    fab = np.array([p.coverage.f_ab for p in scored], dtype=float)

    for b in range(SPAN_BUCKETS):
        mask = bucket == b
        if not mask.any():
            continue
        values = fab[mask]
        order = values.argsort(kind="mergesort")
        ranks = np.empty(len(values), dtype=float)

        i = 0
        srt = values[order]
        while i < len(srt):
            j = i
            while j + 1 < len(srt) and srt[j + 1] == srt[i]:
                j += 1
            ranks[order[i:j + 1]] = (i + j) / 2
            i = j + 1
        pct = ranks / max(len(values) - 1, 1)
        for pair, value in zip([p for p, m in zip(scored, mask) if m], pct):
            pair.coverage.span_percentile = round(float(value), 4)
    return len(scored)


def decide(
    pool: list[Ranked],
    pairs: list[Pair],
    tumour_tpm: dict[str, float] | None = None,
    tolerances: "routing.Tolerances | None" = None,
    per_organ: dict[str, dict[str, float]] | None = None,
) -> list[Decision]:
    """Single, dual, unresolved or no design, per pool member."""
    by_gene: dict[str, list[Pair]] = {r.gene: [] for r in pool}
    for p in pairs:
        by_gene[p.gene_a].append(p)
        by_gene[p.gene_b].append(p)

    cleared = {r.gene: r.cleared for r in pool}
    accession_of = {r.gene: r.accession for r in pool}
    ranked_of = {r.gene: r for r in pool}
    risk_of = {r.gene: (r.risk, r.risk_organ) for r in pool}

    tol = tolerances

    def eligible_partner(gene: str) -> bool:
        """Whether the gene is expressed enough in tumour to serve as a partner."""
        if tumour_tpm is None:
            return True
        value = tumour_tpm.get(gene)
        return value is not None and value >= PARTNER_MIN_TUMOUR_TPM

    out: list[Decision] = []
    for index, r in enumerate(pool):
        mine = by_gene[r.gene]
        admissible = [
            p for p in mine
            if p.admissible and eligible_partner(_other(p, r.gene))
        ]

        admissible.sort(key=lambda p: (p.risk.combined, _other(p, r.gene)))

        evidence = dict(
            evidence_class=r.evidence_class,
            tumour_side_verdict=r.tumour_side_verdict,
            normal_tissue_risk=r.risk,
            normal_tissue_risk_organ=r.risk_organ,
            evidence_confidence=r.confidence,
            protein_arm_measured=r.protein_arm_measured,
            risk_basis=r.risk_basis,
            risk_is_lower_bound=r.risk_is_lower_bound,
            composite=r.composite,
            composite_supported=r.composite_supported,
            measured_weight=r.measured_weight,
            components_measured=sum(
                1 for c in r.components.values() if c.measured),
            components_total=len(r.components),
            target_confidence=r.confidence,
            risk_organs_measured=(
                len(per_organ.get(r.accession, {}))
                if per_organ is not None else None
            ),
        )

        def pair_evidence(chosen):
            """The same fields for the pair, empty where no pair was chosen."""
            if chosen is None:
                return dict(
                    partner_evidence_class=None,
                    partner_measured_weight=None,
                    pair_confidence=None,
                    pair_organs_total=None,
                    pair_organs_resolved=None,
                    pair_coverage_measured=None,
                )
            other = ranked_of.get(_other(chosen, r.gene))
            return dict(
                partner_evidence_class=other.evidence_class if other else None,
                partner_measured_weight=other.measured_weight if other else None,
                pair_confidence=chosen.confidence,
                pair_organs_total=chosen.organs_total,
                pair_organs_resolved=chosen.organs_resolved,
                pair_coverage_measured=chosen.coverage.measured,
            )

        options = dict(partner_options=len(admissible),
                       partner_forced=len(admissible) == 1)

        risk, organ = risk_of.get(r.gene, (None, None))
        best_pair_risk = admissible[0].risk.combined if admissible else None
        if tol is None:
            decided = None
            route_fields = dict(
                architecture=routing.NOT_CONFIGURED,
                route_reason="no tolerances supplied; routing disabled",
                route_ceiling=None, route_exposure=None,
            )
        else:
            decided = routing.route(
                gene=r.gene, risk=risk, risk_organ=organ, tolerances=tol,
                pair_risk=best_pair_risk,
                partner=_other(admissible[0], r.gene) if admissible else None,
            )
            route_fields = dict(
                architecture=decided.architecture,
                route_reason=decided.reason,
                route_ceiling=decided.ceiling,
                route_exposure=decided.exposure,
            )

        if cleared[r.gene]:
            best = admissible[0] if admissible else None
            out.append(
                Decision(
                    gene=r.gene,
                    outcome=SINGLE,
                    **route_fields,
                    **options,
                    **evidence,
                    **pair_evidence(best),
                    partner=_other(best, r.gene) if best else None,
                    pair=best,
                    accession=r.accession,
                    partner_accession=(
                        accession_of.get(_other(best, r.gene)) if best else None
                    ),
                    pool_index=index,
                )
            )
            continue

        if admissible:
            best = admissible[0]
            out.append(
                Decision(
                    gene=r.gene,
                    outcome=DUAL,
                    **route_fields,
                    **options,
                    **evidence,
                    **pair_evidence(best),
                    partner=_other(best, r.gene),
                    pair=best,
                    accession=r.accession,
                    partner_accession=accession_of.get(_other(best, r.gene)),
                    pool_index=index,
                )
            )
            continue

        if decided is not None and decided.architecture == routing.ADAPTOR:
            out.append(
                Decision(
                    gene=r.gene,
                    outcome=ADAPTOR,
                    **evidence,
                    **pair_evidence(None),
                    partner=None,
                    accession=r.accession,
                    pool_index=index,
                    **route_fields,
                )
            )
            continue

        failed = {
            "risk": sum(1 for p in mine if not p.cleared),
            "coverage_below_floor": sum(
                1 for p in mine
                if p.cleared and p.coverage.measured and not p.coverage_ok
            ),
            "unmeasured": sum(1 for p in mine if not p.coverage.measured),
            "partner_ineligible": sum(
                1 for p in mine
                if p.admissible and not eligible_partner(_other(p, r.gene))
            ),
        }

        salvageable = [
            p
            for p in mine
            if p.risk.risk_unresolved
            and p.risk.optimistic is not None
            and p.risk.optimistic <= p.ceiling
            and p.coverage.measured
            and eligible_partner(_other(p, r.gene))
        ]
        out.append(
            Decision(
                gene=r.gene,
                outcome=UNRESOLVED if salvageable else NO_DESIGN,
                failed_on=failed,
                **route_fields,
                **evidence,
                **pair_evidence(salvageable[0] if salvageable else None),
                accession=r.accession,
                pool_index=index,
            )
        )
    return out


def _other(pair: Pair, gene: str) -> str:
    """The partner on the far side of a pair."""
    return pair.gene_b if pair.gene_a == gene else pair.gene_a


def configuration_hash(
    stage3_hash: str,
    pool_genes: list[str],
    tolerances: "routing.Tolerances | None" = None,
) -> str:
    """Covers the stage 3 hash as well as this stage's own parameters."""
    payload = {
        "stage3": stage3_hash,
        "routing": (routing.configuration_payload(tolerances)
                    if tolerances is not None else None),
        "pool_size": POOL_SIZE,
        "pool": pool_genes,
        "detection_counts": DETECTION_COUNTS,
        "coverage_floor": COVERAGE_FLOOR,
        "patient_fraction_floor": PATIENT_FRACTION_FLOOR,
        "min_malignant_cells": MIN_MALIGNANT_CELLS,
        "min_detected_cells": MIN_DETECTED_CELLS,
        "selection_rule": SELECTION_RULE,
        "partner_min_tumour_tpm": PARTNER_MIN_TUMOUR_TPM,
        "span_buckets": SPAN_BUCKETS,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


DECISIONS_KEY = "decisions"


DECISIONS_MANIFEST_VERSION = 1


def _decision_payload(d: Decision) -> dict:
    """One decision, flattened for storage."""
    return {
        "gene": d.gene,
        "accession": d.accession,
        "outcome": d.outcome,
        "partner": d.partner,
        "partner_accession": d.partner_accession,
        "pool_index": d.pool_index,
        "partner_options": d.partner_options,
        "partner_forced": d.partner_forced,
        "failed_on": d.failed_on or None,
        "architecture": d.architecture or None,
        "route_reason": d.route_reason or None,
        "route_ceiling": d.route_ceiling,
        "route_exposure": d.route_exposure,
        "evidence_class": d.evidence_class or None,
        "tumour_side_verdict": d.tumour_side_verdict or None,
        "normal_tissue_risk": d.normal_tissue_risk,
        "normal_tissue_risk_organ": d.normal_tissue_risk_organ,
        "evidence_confidence": d.evidence_confidence,
        "protein_arm_measured": d.protein_arm_measured,
        "risk_basis": d.risk_basis or None,
        "risk_is_lower_bound": d.risk_is_lower_bound,
        "composite": d.composite,
        "composite_supported": d.composite_supported,
        "measured_weight": d.measured_weight,
        "components_measured": d.components_measured,
        "components_total": d.components_total,
        "target_confidence": d.target_confidence,
        "risk_organs_measured": d.risk_organs_measured,
        "partner_evidence_class": d.partner_evidence_class,
        "partner_measured_weight": d.partner_measured_weight,
        "pair_confidence": d.pair_confidence,
        "pair_confidence_for": (
            f"{d.gene}+{d.partner}" if d.partner and d.pair_confidence is not None
            else None),
        "pair_organs_total": d.pair_organs_total,
        "pair_organs_resolved": d.pair_organs_resolved,
        "pair_coverage_measured": d.pair_coverage_measured,
    }


def decision_rows(decisions: list[Decision]) -> list[dict]:
    """Decisions as plain rows, the same shape `read_decisions` returns."""
    return [_decision_payload(d) for d in decisions]


def write_decisions(
    decisions: list[Decision],
    pool_genes: list[str],
    stage3_hash: str,
    criteria: dict[str, bool],
    root=None,
):
    """Persist the decisions, with the hashes and criteria that produced them."""
    from car_pipeline.data.source import CACHE_ROOT, _write_json_atomic

    base = (root or CACHE_ROOT) / "stage4"
    payload_path = base / (DECISIONS_KEY + ".json")
    manifest_path = base / (DECISIONS_KEY + ".manifest.json")

    if not criteria:
        raise ValueError(
            "criteria outcomes are required: an artifact written without them "
            "would assert usable_as_result with nothing behind it"
        )

    if manifest_path.exists():
        manifest_path.unlink()

    config_hash = configuration_hash(stage3_hash, pool_genes)
    rows = [_decision_payload(d) for d in decisions]
    _write_json_atomic(payload_path, {"decisions": rows})

    blob = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    tripped = sorted(k for k, ok in (criteria or {}).items() if not ok)
    _write_json_atomic(
        manifest_path,
        {
            "key": DECISIONS_KEY,
            "manifest_version": DECISIONS_MANIFEST_VERSION,
            "rows": len(rows),
            "digest": hashlib.sha256(blob).hexdigest(),
            "stage3_hash": stage3_hash,
            "stage4_hash": config_hash,
            "pool_size": len(pool_genes),
            "outcomes": {
                name: sum(1 for d in decisions if d.outcome == name)
                for name in (SINGLE, DUAL, NO_DESIGN, UNRESOLVED)
            },
            "criteria_tripped": tripped,
            "usable_as_result": not tripped,
        },
    )
    return payload_path


def read_decisions(
    root=None,
    allow_unusable: bool = False,
    expect_stage3_hash: str | None = None,
    expect_stage4_hash: str | None = None,
) -> tuple[list[dict], dict]:
    """Read the decisions, refusing a payload no manifest blesses."""
    from car_pipeline.data.source import CACHE_ROOT, CacheError

    base = (root or CACHE_ROOT) / "stage4"
    payload_path = base / (DECISIONS_KEY + ".json")
    manifest_path = base / (DECISIONS_KEY + ".manifest.json")
    if not manifest_path.exists():
        raise CacheError(
            "no manifest beside " + str(payload_path)
            + "; the writing run did not finish"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = json.loads(payload_path.read_text(encoding="utf-8"))["decisions"]
    blob = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    if hashlib.sha256(blob).hexdigest() != manifest["digest"]:
        raise CacheError(
            str(payload_path) + " does not match its manifest digest; refusing to read"
        )
    if len(rows) != manifest["rows"]:
        raise CacheError(
            str(payload_path) + " has " + str(len(rows)) + " rows, manifest says "
            + str(manifest["rows"])
        )
    version = manifest.get("manifest_version")
    if version != DECISIONS_MANIFEST_VERSION:
        raise CacheError(
            str(payload_path) + " was written under manifest version "
            + str(version) + "; this reader expects "
            + str(DECISIONS_MANIFEST_VERSION)
        )

    for label, expected, actual in (
        ("stage3_hash", expect_stage3_hash, manifest.get("stage3_hash")),
        ("stage4_hash", expect_stage4_hash, manifest.get("stage4_hash")),
    ):
        if expected is not None and expected != actual:
            raise CacheError(
                str(payload_path) + " was written under " + label + " "
                + str(actual) + ", but the caller is running " + str(expected)
                + "; these decisions are not from this configuration"
            )
    if not manifest.get("usable_as_result") and not allow_unusable:
        raise CacheError(
            str(payload_path) + " was written by a run that tripped "
            + ", ".join(manifest.get("criteria_tripped", []))
            + "; pass allow_unusable=True to read it as diagnostics rather than "
            "as a result"
        )
    return rows, manifest
