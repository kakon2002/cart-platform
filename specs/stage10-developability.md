# Stage 10 — developability

Written before `stages/stage10.py`. Mechanical: arithmetic over sequences that
already exist. **No new sources.**

Status: **implemented.** Criteria below were fixed before the run.

---

## 1. What it scores, and what it cannot

Stage 6 produces **zero buildable constructs**, so scoring constructs would score
nothing. This stage scores the **binder variable regions** retrieved by Stage 5,
which do exist — 30 of the 200 pool targets carry one.

Every measure is computed from the amino acid sequence alone. None is a
prediction of manufacturing failure, and none is called one. They are **sequence
liabilities**: properties known to correlate with difficulty, reported as counts
and values so a reader can weigh them rather than receive a verdict.

| measure | what it is | why it is a liability |
| --- | --- | --- |
| isoelectric point | pH of zero net charge, by bisection over a fixed pKa table | a pI near the formulation pH tends to lower solubility |
| net charge at pH 7.4 | same table | very low absolute charge tracks aggregation |
| unpaired cysteine | count and parity | an odd count guarantees one cysteine unpaired |
| N-glycosylation sequons | `N-X-S/T`, X not proline | added heterogeneity, and a sequon in a binding loop can block binding |
| aggregation-prone regions | windows of 7 with mean hydropathy ≥ 1.0 | contiguous hydrophobic surface drives self-association |
| hydropathy (GRAVY) | mean Kyte–Doolittle | whole-sequence hydrophobicity |

**Cysteine pairing cannot be determined from sequence.** An even count is not
evidence of full pairing; it is only the absence of a guarantee of the opposite.
The output says `parity`, never `unpaired: 0`.

## 2. Thresholds, fixed here

`pI within 1.0 of 7.4` · `|net charge| < 1.0` · `odd cysteine count` ·
`≥ 1 sequon` · `≥ 1 aggregation-prone region`. Each raises one flag. **Flags are
counted and listed, never summed into a score** — a single number would let a
strong liability be averaged away by four weak absences.

## 3. When the input is empty

If no binder carries a sequence, this stage emits **`NOTHING_TO_SCORE`** with the
count it expected and the reason, and reports zero rows. **It does not emit an
empty table.** An empty table is read as "nothing had liabilities", which is the
opposite of "nothing was examined".

## 4. Criteria

Positive pins on synthetic controls with hand-checkable answers, so the arithmetic
is pinned independently of the data.

| id | criterion |
| --- | --- |
| **D1** | **a poly-lysine control does not score basic, or poly-glutamate does not score acidic** — pins the charge model to a known answer |
| **D2** | **`NST` does not yield one sequon, or `NPT` yields any** — pins the proline exclusion, which is the rule most easily dropped |
| **D3** | a 3-cysteine control does not report odd parity |
| **D4** | any scored binder has a pI outside 1–14, or a flag count outside 0–5 |
| **D5** | the row count differs from the number of binders carrying a sequence |
| **D6** | any liability is summed into a single score |
