# Stage 11 — multi-objective ranking

Written before `stages/stage11.py`. Combines numbers the earlier stages already
produced. **No new measurement.**

Status: **implemented.** Criteria below were fixed before the run.

---

## 1. No weighted sum

The objectives are not commensurable and are not made so. Tumour attractiveness,
safety margin, binder availability and sequence liability are reported as a
**Pareto front**: a design is on the front when nothing else is at least as good
on every objective and better on one.

A weighted sum would require weights this project has no basis to invent, and
would let a design with an unmanageable safety margin be rescued by a strong
tumour score — the exact failure Stage 3 §1 refuses by keeping its three numbers
apart.

Objectives, all carried:

| objective | source | direction |
| --- | --- | --- |
| tumour attractiveness | Stage 3 composite | higher better |
| safety margin | ceiling − Stage 3 risk | higher better |
| binder availability | Stage 5 candidate count | higher better |
| sequence cleanliness | negated Stage 10 flag count | higher better |

## 2. The attrition chain is the primary output

Every one of the 200 pool members is attributed to **the first gate it fails**, in
pipeline order: safety, then recommendation, then binder, then construct, then
budget. The counts sum to 200 and a criterion checks that they do.

This is the primary output because for this indication it is the only one with
content. The ranking is a table of whatever survives; the chain is the account of
why so little does.

## 3. When nothing survives

If no design reaches the end, this stage emits
**`NO_DESIGN_REACHES_THE_END`** with the full attrition chain and the reason at
each step.

**It does not emit an empty ranking.** An empty table reads as "nothing ranked
highly"; the true statement is "nothing arrived to be ranked", and the difference
is the entire result of this pipeline for this indication. The status is the
finding, not the absence of rows.

## 4. Criteria

| id | criterion |
| --- | --- |
| **N1** | **a synthetic dominated point appears on the Pareto front** — pins the front's logic to a hand-checkable answer rather than to the data |
| **N2** | **a synthetic non-dominated point is missing from the front** — the same pin in the other direction |
| N3 | the attrition counts do not sum to the pool size |
| N4 | any weighted or summed score across objectives is emitted |
| N5 | an empty result is emitted as an empty ranking rather than as a named status |
| N6 | the output row count differs from the pool size |
