"""Stage 4 — target pairing and the single-versus-dual decision.

Implements `specs/stage4-target-pairing.md`. Pool, thresholds and criteria are
fixed there and are read from here, not tuned.

Two new measurements, and nothing else is re-derived:

* combined risk, the minimum per organ before the maximum across organs, which
  is the only thing an AND gate buys over either antigen alone
* co-expression, the fraction of malignant cells carrying both antigens, which
  is what the gate actually kills

The per-organ scores come from the ranking stage rather than being recomputed,
so pairing a target with itself reproduces its single-antigen risk exactly. A
second implementation of the tissue mapping would be a second place for the
mapping bugs to live.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

from car_pipeline.stages import stage3
from car_pipeline.stages.stage3 import RiskModel, Ranked

#: The pool is the top of the tumour-side ranking with risk ignored entirely.
#: Risk is the thing pairing exists to fix, so filtering on it first would leave
#: the stage unable to reach the targets a dual design is for.
POOL_SIZE = 200

#: A cell carries an antigen when at least one molecule was captured. Reported
#: at 2 and 3 as well, because the ordering is not stable across them.
DETECTION_COUNTS = 1
SENSITIVITY_COUNTS = (1, 2, 3)

#: Set against the measured range rather than chosen and then discovered to
#: admit nothing: the best pair of known targets in this atlas reaches 0.047 of
#: malignant cells, so a floor at 0.05 or 0.10 would eliminate every pair this
#: stage exists to evaluate. Deliberately low, and paired with P16, which trips
#: if even this admits nothing.
COVERAGE_FLOOR = 0.02
PATIENT_FRACTION_FLOOR = 0.60

#: A proportion cannot be estimated from three cells. 29 of 43 patients clear
#: this; the excluded ones are reported rather than dropped quietly.
MIN_MALIGNANT_CELLS = 100

#: Below this a positive fraction is not a measurement. Carrying forward the
#: rule that a single-cell zero never rejects a target: such a pair is marked
#: unmeasured, not scored zero.
MIN_DETECTED_CELLS = 10

SINGLE = "SINGLE"
DUAL = "DUAL"
NO_DESIGN = "NO_DESIGN"
UNRESOLVED = "UNRESOLVED"


# --------------------------------------------------------------------------
# combined risk
# --------------------------------------------------------------------------


@dataclass
class PairRisk:
    combined: float | None            # the gate: min per organ, unresolved filled
    optimistic: float | None          # unresolved organs contribute nothing
    independence: float | None        # score_A x score_B, the lower bound
    #: Unresolved organs charged at full criticality rather than at the measured
    #: member's score. Not the gate, and not a candidate for it: it is the
    #: question "would this pair still clear if the unmeasured antigen turned
    #: out to saturate the organ", which is what "clearance depends on an
    #: unresolved organ" means when written as a number.
    pessimistic: float | None
    organ: str | None
    unresolved_organs: list[str] = field(default_factory=list)

    @property
    def risk_unresolved(self) -> bool:
        return bool(self.unresolved_organs)


def pair_risk(
    model: RiskModel, a: dict[str, float], b: dict[str, float]
) -> PairRisk:
    """Minimum per organ, then maximum across organs.

    Minimum because an AND gate only fires where both antigens are present, so
    the organ's risk is bounded by whichever is scarcer there. Taking the
    maximum would reduce the pair to its more dangerous member and ignore the
    architecture entirely.

    The minimum assumes perfect overlap within the organ, which is maximal
    co-expression rather than minimal. That is pessimistic for safety and
    optimistic for coverage, and the two coincide, which is what makes it the
    right thing to gate on. It is a bound: neither source carries a joint
    distribution over cells, so no measurement of within-organ co-expression
    exists to be had.

    Organs measured for neither member are outside both mappings and contribute
    to nothing. The ranking stage already tolerates organs nobody measured, and
    filling them here would break the identity that pairing a target with itself
    reproduces its own risk.
    """
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
        # One member measured here and the other not. The pair's whole claim is
        # an absence, so an unobserved absence cannot be credited to it: the
        # missing member is assumed present wherever nobody looked.
        unresolved.append(organ)
        known = sa if sa is not None else sb
        conservative[organ] = known
        independence[organ] = known
        pessimistic[organ] = 1.0

    combined, organ = stage3.worst_organ(model, conservative)
    opt, _ = stage3.worst_organ(model, optimistic)
    ind, _ = stage3.worst_organ(model, independence)
    pess, _ = stage3.worst_organ(model, pessimistic)
    # Rounded to the precision the ranking stage stores its own risk at, so the
    # two are directly comparable. Comparing a full precision pair risk against
    # a rounded single risk reports a disagreement at the fifth decimal that is
    # arithmetic rather than substance.
    return PairRisk(
        _round(combined),
        _round(opt),
        _round(ind),
        _round(pess),
        organ,
        sorted(unresolved),
    )


def _round(value: float | None) -> float | None:
    return None if value is None else round(value, 4)


# --------------------------------------------------------------------------
# co-expression
# --------------------------------------------------------------------------


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

    #: Genomic span context. `f_ab` tracks how long the two genes are more
    #: strongly than how much of them is expressed (§6.5b), so the raw fraction
    #: cannot be read on its own. `span_percentile` is where this pair's `f_ab`
    #: falls among measured pairs of similar span: 0.50 means typical for genes
    #: this long, which is a different statement from the fraction itself.
    #: Both are reported and neither gates anything.
    span_geomean_kb: float | None = None
    span_percentile: float | None = None

    @property
    def patient_fraction(self) -> float:
        if not self.patients_evaluable:
            return 0.0
        return self.patients_at_floor / self.patients_evaluable

    # -- what each architecture would reach --------------------------------
    #
    # Named quantities rather than things a reader derives from two columns.
    # The AND gate's coverage is the intersection; a single target is its own
    # marginal; an OR gate is the union. The differences between those three are
    # the price of the safety the AND gate buys, and they belong in the output.

    @property
    def escape(self) -> float | None:
        """Malignant cells the AND gate does not engage.

        A floor, not an estimate. This assay drops transcripts, so the measured
        intersection understates the true one and this number overstates the
        true escape population.
        """
        return None if self.f_ab is None else 1.0 - self.f_ab

    @property
    def or_gate(self) -> float | None:
        """Union: what a design firing on either antigen would reach."""
        if self.f_ab is None:
            return None
        return self.f_a + self.f_b - self.f_ab

    @property
    def best_single(self) -> float | None:
        return None if self.f_ab is None else max(self.f_a, self.f_b)

    @property
    def coverage_cost(self) -> float | None:
        """How many times more of the tumour the better single target reaches."""
        if not self.f_ab:
            return None
        return self.best_single / self.f_ab


def intersection_matrix(positive: np.ndarray) -> np.ndarray:
    """Pairwise double-positive counts for every gene pair at once.

    `positive` is cells x genes. The product of its transpose with itself gives
    every intersection in one pass; the diagonal is the marginal count. Exact in
    float32 because the cell count is far below the point where it stops
    counting integers exactly.
    """
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
    if n_a < MIN_DETECTED_CELLS or n_b < MIN_DETECTED_CELLS:
        # This assay drops transcripts that bulk measurement finds abundantly
        # present, so silence here is the assay reporting its own capture
        # failure. Unmeasured, never zero, and never a rejection.
        return Coverage(measured=False, reason="below detection")

    f_a, f_b = n_a / n_cells, n_b / n_cells
    f_ab = n_ab / n_cells
    union = n_a + n_b - n_ab
    return Coverage(
        measured=True,
        f_a=f_a,
        f_b=f_b,
        f_ab=f_ab,
        # Named for what survives: P(B|A) is the share of A's own positive cells
        # that also carry B, so it is A's coverage that it describes.
        p_a_given_b=n_ab / n_b,
        p_b_given_a=n_ab / n_a,
        sacrificed_a=1.0 - n_ab / n_a,
        sacrificed_b=1.0 - n_ab / n_b,
        jaccard=(n_ab / union) if union else 0.0,
        patients_at_floor=patients_at_floor,
        patients_evaluable=patients_evaluable,
    )


# --------------------------------------------------------------------------
# pairs
# --------------------------------------------------------------------------


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
        """Third number, never combined with the other two.

        Bounded by the weaker member by construction: a pair cannot be better
        evidenced than the least evidenced antigen in it. The two discounts are
        the things pairing adds to the question — how much of the organ union
        was resolved for both members, and whether co-expression was measurable
        at all.
        """
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
        if self.risk.combined is None:
            return None
        return self.risk_a - self.risk.combined

    @property
    def delta_b(self) -> float | None:
        if self.risk.combined is None:
            return None
        return self.risk_b - self.risk.combined

    @property
    def rescued(self) -> list[str]:
        """A condition, not a statistic.

        A member is rescued only when its own risk is above the ceiling and the
        pair's is not. A large movement that does not cross buys nothing and is
        not a rescue; the delta still reports how far it got.
        """
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
        c = self.coverage
        return (
            c.measured
            and c.f_ab >= COVERAGE_FLOOR
            and c.patient_fraction >= PATIENT_FRACTION_FLOOR
        )

    @property
    def admissible(self) -> bool:
        """Risk-gated and measurable. Coverage does not gate.

        `f_ab` is confounded with genomic span: over the pool its rank
        correlation with span is +0.68 against +0.20 with bulk expression, the
        confound reaches the joint quantity (+0.63 for `f_ab` itself against
        +0.08 for expression), and it survives stratification by expression. A
        threshold on it therefore admits and rejects partners substantially on
        how long their genes are.

        `coverage.measured` is still required, and is a different question: it
        asks whether the co-expression was observable at all, not whether it
        cleared a number. An unmeasured pair cannot be recommended.

        What this gives up is stated rather than hidden: nothing now stops a pair
        with negligible overlap being recommended, so `f_ab` and its span
        percentile are reported per pair and a reader has to look at them.
        """
        return self.cleared and self.coverage.measured


def build_pool(rows: list[Ranked], size: int = POOL_SIZE) -> list[Ranked]:
    """Top of the tumour-side ranking, risk ignored entirely.

    One entry per symbol: several symbols carry more than one accession, and two
    pool members naming the same gene would pair with each other and report a
    perfect intersection that means nothing.
    """
    scored = [r for r in rows if r.composite is not None and r.gene]
    scored.sort(key=lambda r: (-r.composite, r.gene, r.accession))
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
        # The pool ignores the *value* of risk, not whether it exists. A member
        # with no risk at all would make every pair containing it unresolvable,
        # and silently so.
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

    # Per patient, accumulated as a count of patients at or above the floor for
    # each pair. A pair double-positive in half the patients and absent in the
    # rest pools identically to one that is uniform, and those are different
    # products.
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


# --------------------------------------------------------------------------
# the decision
# --------------------------------------------------------------------------


@dataclass
class Decision:
    gene: str
    outcome: str
    partner: str | None = None
    pair: Pair | None = None
    failed_on: dict[str, int] = field(default_factory=dict)

    #: Carried so the decision can be read without the pool beside it. A symbol
    #: is not an identity here — several symbols in this proteome carry more
    #: than one accession, and `build_pool` keeps exactly one of them. A
    #: consumer re-deriving the accession from the symbol would sometimes pick
    #: the other one and would never be told.
    accession: str = ""
    partner_accession: str | None = None
    #: Position in the pool as Stage 4 ordered it, so the order survives the
    #: round trip and can be checked rather than trusted.
    pool_index: int = -1


#: Span buckets for the percentile. Quintiles: enough to separate a 7 kb gene
#: from a 1.1 Mb one without slicing the pool so finely that a bucket holds too
#: few pairs to rank within.
SPAN_BUCKETS = 5

#: A partner must carry this much of the antigen in the tumour itself. Applied to
#: the partner only, and the asymmetry is the point: a target earns its place
#: through the tumour-side composite, which already scores expression and
#: prevalence, while a partner is chosen purely on risk and would otherwise be
#: rewarded for being absent. An AND gate fires only where both antigens are
#: present, so a partner absent from the tumour contributes nothing to killing it
#: while contributing everything to the pair looking safe.
#:
#: Measured on bulk tumour transcript level, which is the axis that is NOT
#: confounded with genomic span — the per-cell measure is (§6.5b), and using it
#: here would reintroduce the artefact this exists to work around.
#:
#: 5.0 rather than 3.0, and the reason is not margin. The concentration this
#: addresses came from one protein sitting far below every other candidate:
#: 0.0277 against 0.2272 for the next lowest. At a 3 TPM threshold the lowest
#: eligible partner still leads the next by 0.0489. At 5 TPM the leaders cluster
#: within 0.0036 of each other, so no single gene can win for every target. The
#: threshold sits at roughly the pool's 8th percentile and retains 182 of 200; it
#: is not fitted to exclude two named genes.
PARTNER_MIN_TUMOUR_TPM = 5.0

#: Names what admits and orders a partner, so a change to either shows up in the
#: configuration hash rather than being invisible to it.
SELECTION_RULE = (
    "risk-cleared-and-measured;partner_min_tumour_tpm;"
    "order:combined_risk,partner_name;v3"
)


def annotate_span_context(pairs: list[Pair], spans: dict[str, int]) -> int:
    """Attach each measured pair's span and its within-span percentile.

    Reporting only. The percentile answers "is this overlap high for genes this
    long", which is the question the raw fraction cannot answer while detection
    tracks span. A pair whose members have no span on record is left with both
    fields None rather than assigned a middle value.

    Returns the number of pairs annotated.
    """
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
        # Average rank across ties, so a block of identical f_AB values does not
        # get an ordering the data does not support.
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
) -> list[Decision]:
    """Single, dual, unresolved or no design, per pool member.

    `tumour_tpm` supplies bulk tumour transcript level per gene and gates partner
    eligibility at `PARTNER_MIN_TUMOUR_TPM`. A gene with no measurement is not
    eligible as a partner: absence of evidence that the partner is on the tumour
    is exactly the case the threshold exists for, and treating it as a pass would
    put the missing-is-a-third-state rule the wrong way round. Passing None
    disables the filter entirely, which is for measuring what it does, not for
    running without it.
    """
    by_gene: dict[str, list[Pair]] = {r.gene: [] for r in pool}
    for p in pairs:
        by_gene[p.gene_a].append(p)
        by_gene[p.gene_b].append(p)

    cleared = {r.gene: r.cleared for r in pool}
    accession_of = {r.gene: r.accession for r in pool}

    def eligible_partner(gene: str) -> bool:
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
        # Ordered by how far under the ceiling the pair sits, then by partner
        # name so the choice is deterministic. This used to order by co-expression
        # — how much of the tumour the gate still kills — which is the better
        # question and the one `f_ab` cannot currently answer: it is confounded
        # with genomic span (see `admissible`). Ordering on the risk margin is a
        # weaker criterion honestly measured, rather than a stronger one measured
        # on an artefact.
        admissible.sort(key=lambda p: (p.risk.combined, _other(p, r.gene)))

        if cleared[r.gene]:
            best = admissible[0] if admissible else None
            out.append(
                Decision(
                    gene=r.gene,
                    outcome=SINGLE,
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
                    partner=_other(best, r.gene),
                    pair=best,
                    accession=r.accession,
                    partner_accession=accession_of.get(_other(best, r.gene)),
                    pool_index=index,
                )
            )
            continue

        # Which wall each target hit, and the coverage row is now what it says.
        # It counts pairs that cleared risk and were measured but fell under the
        # reported floor — which no longer excludes anything, so a non-zero count
        # here is information about the pair, not a reason it was rejected. Left
        # separate from `unmeasured` rather than folded into it: a pair nobody
        # could measure and a pair measured and found thin are different facts.
        failed = {
            "risk": sum(1 for p in mine if not p.cleared),
            "coverage_below_floor": sum(
                1 for p in mine
                if p.cleared and p.coverage.measured and not p.coverage_ok
            ),
            "unmeasured": sum(1 for p in mine if not p.coverage.measured),
            # Pairs that cleared everything and were rejected only because the
            # partner is not expressed enough in the tumour. Without this row a
            # target excluded entirely on partner eligibility is persisted with
            # every counter at zero — a rejection with no stated reason.
            "partner_ineligible": sum(
                1 for p in mine
                if p.admissible and not eligible_partner(_other(p, r.gene))
            ),
        }
        # `coverage.measured`, not `coverage_ok`: the coverage floor no longer
        # selects (§6.5b), and leaving it here would decide NO_DESIGN against
        # UNRESOLVED on a threshold the stage has stopped applying anywhere else.
        # Eligibility applies here as well. A target whose only salvageable
        # pairs run through partners the tumour-expression gate rejects is not
        # salvageable: resolving the organ would not make that partner usable,
        # so reporting UNRESOLVED would promise a design that cannot follow.
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
                accession=r.accession,
                pool_index=index,
            )
        )
    return out


def _other(pair: Pair, gene: str) -> str:
    return pair.gene_b if pair.gene_a == gene else pair.gene_a


# --------------------------------------------------------------------------
# reproducibility
# --------------------------------------------------------------------------


def configuration_hash(stage3_hash: str, pool_genes: list[str]) -> str:
    """Covers the stage 3 hash as well as this stage's own parameters.

    A pairing result is not interpretable without knowing which ranking produced
    its inputs, so the two travel together.
    """
    payload = {
        "stage3": stage3_hash,
        "pool_size": POOL_SIZE,
        "pool": pool_genes,
        "detection_counts": DETECTION_COUNTS,
        # Retained because both still bound the reported coverage numbers, even
        # though neither selects any more.
        "coverage_floor": COVERAGE_FLOOR,
        "patient_fraction_floor": PATIENT_FRACTION_FLOOR,
        "min_malignant_cells": MIN_MALIGNANT_CELLS,
        "min_detected_cells": MIN_DETECTED_CELLS,
        # What admits and orders a partner. Without this a run from before
        # coverage was removed from selection hashes identically to one after,
        # and `read_decisions(expect_stage4_hash=...)` would accept the old
        # artifact as current — the one thing carrying the hash is meant to stop.
        "selection_rule": SELECTION_RULE,
        "partner_min_tumour_tpm": PARTNER_MIN_TUMOUR_TPM,
        "span_buckets": SPAN_BUCKETS,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


# --------------------------------------------------------------------------
# persistence
# --------------------------------------------------------------------------
#
# Until now this stage returned its decisions and printed a summary, so nothing
# downstream could read what it decided without re-running it — and re-running
# it re-derives numbers that are themselves under question. The artifact is
# written under the same discipline as a source cache: payload first, manifest
# second, both atomic, the manifest acting as the commit marker. A payload with
# no manifest beside it is a run that died mid-write and must not be read.

DECISIONS_KEY = "decisions"

#: Bumped when the payload's shape changes. Without it an artifact written by an
#: older layout reads as current and fails somewhere further away.
DECISIONS_MANIFEST_VERSION = 1


def _decision_payload(d: Decision) -> dict:
    return {
        "gene": d.gene,
        "accession": d.accession,
        "outcome": d.outcome,
        "partner": d.partner,
        "partner_accession": d.partner_accession,
        "pool_index": d.pool_index,
        # Present only where the outcome is terminal. An absent mapping and a
        # mapping of zeros mean different things and are kept apart.
        "failed_on": d.failed_on or None,
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
    """Persist the decisions, with the hashes and criteria that produced them.

    The criteria mapping is written into the manifest rather than left to the
    reader's memory. These decisions are currently produced by a run that stops
    on five tripped criteria, and an artifact that did not say so would be read
    as a result. A consumer is expected to refuse the payload when anything is
    tripped, which is why the outcome is stored beside the data and not in a log.
    """
    from car_pipeline.data.source import CACHE_ROOT, _write_json_atomic

    base = (root or CACHE_ROOT) / "stage4"
    payload_path = base / (DECISIONS_KEY + ".json")
    manifest_path = base / (DECISIONS_KEY + ".manifest.json")

    # Validated before anything on disk is touched. Raising after the unlink
    # would destroy a previously valid artifact to punish a bad call, which
    # turns a caller's mistake into data loss.
    if not criteria:
        raise ValueError(
            "criteria outcomes are required: an artifact written without them "
            "would assert usable_as_result with nothing behind it"
        )

    # A stale manifest must never bless a new payload. Removed first, so a crash
    # between the two writes leaves an unblessed payload rather than a blessed
    # mismatch, and the reader below refuses the first and cannot detect the
    # second.
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
    """Read the decisions, refusing a payload no manifest blesses.

    Returns the rows *and* the manifest, because the rows alone do not say
    whether they may be read as a result and a caller handed only the rows
    cannot find out. Today they may not: the writing run stops on five tripped
    criteria. A caller that genuinely wants the decisions anyway — to inspect
    why the run stopped, which is a legitimate thing to want — has to say so.

    The digest is re-derived rather than trusted. A truncated or hand-edited
    payload that still parses as JSON is exactly the failure this guards, and it
    is silent without the check.
    """
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
    # The hashes are recorded so a consumer can tell whether these decisions came
    # from the configuration it is holding. Recording them and never checking
    # them would let an artifact from a different ranking be read as current,
    # which is precisely what carrying the hashes was supposed to prevent.
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
