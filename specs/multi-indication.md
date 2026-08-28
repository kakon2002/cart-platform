# Multi-indication

The reference document asks for a platform that accepts any cancer type. The API
accepted one, and every cache underneath it was pinned to pancreas.

## 1. The priority defect: caches collide in place

Fingerprints were indication-aware. **Paths were not.**

```
data/tcga/cohort.npz              fingerprint {"project": "TCGA-PAAD", ...}
data/tcga/cohort.manifest.json    manifest_path = root / f"{entry.key}.manifest.json"
data/singlecell/group_means.npz   fingerprint {"series": "GSE202051", ...}
```

Two indications share one disk slot. Switching invalidates the fingerprint, so
`ensure()` rebuilds and `commit()` overwrites unconditionally through
`os.replace`. Nothing errors.

> **This is the worst failure shape available to this project.** After running
> breast, a PDAC run would find its own fingerprint invalid, rebuild against
> whatever it could reach, and — if a rebuild were skipped or a stale payload
> survived — screen pancreas against the breast atlas and produce a ranked list
> that looks entirely plausible. No exception, no warning, wrong answer.

**The fix is that every indication-specific cache entry carries its indication in
both its filename and its key.** The key matters as much as the filename,
because the manifest path is derived from the key alone: namespacing one and not
the other leaves the manifests colliding, which is the same bug with a longer
fuse.

| entry | was | becomes |
| --- | --- | --- |
| cohort | `cohort.npz` | `cohort__TCGA-BRCA.npz` |
| file index | `file_index.json` | `file_index__TCGA-BRCA.json` |
| atlas archive | `GSE202051_….h5ad.gz` | unchanged — already carries the accession |
| atlas matrix | `totaldata-final-toshare.h5ad` | `matrix__GSE176078.h5ad` |
| group means | `group_means.npz` | `group_means__GSE176078.npz` |
| malignant cells | `malignant_cells_{digest}.npz` | `malignant_cells__GSE176078_{digest}.npz` |
| DepMap lineage | `lineage_matrix.npz` | `lineage_matrix__Breast.npz` |

Sources that describe the human body rather than a tumour — UniProt, GTEx, HPA,
GENCODE, SAbDab, the construct parts — are **not** namespaced. They are the same
measurement whatever the indication, and duplicating them would be a claim that
normal-tissue biology changes with the cancer being screened.

**Criterion M1 exists for exactly this**: run both indications, then assert that
neither one's artifacts moved. It records every indication-specific payload's
digest before and after and fails on any difference.

## 2. Degradation must refuse, not caveat

Measured, by emptying each join and re-running the ranking:

| missing | weight lost | targets still scored | usable? |
| --- | --- | --- | --- |
| DepMap | 0.05 (C6) | 3,400 of 3,466 | **yes** — ranking essentially unchanged |
| single-cell atlas | 0.45 (C1+C2) | 3,399 of 3,466 | **no** |

The atlas row is the important one, and the reason is structural rather than
arithmetic. **C2 (`malignant_vs_stroma`) is the only component that rejects
stromal and immune genes.** Without it nothing in the score distinguishes a
tumour-cell antigen from a gene expressed by the lymphocytes sitting next to the
tumour, so the pool fills with immunoglobulin, T-cell receptor and MHC-II genes:
top-20 overlap with the reference ranking falls to **5 of 20**, and 10 of the
top 20 become Ig/TCR/MHC-II.

Losing 0.45 of weight is survivable arithmetic. Losing the only discriminator
against stroma is not, and the renormalised mean hides that by rescaling what
remains as though it were the whole score.

> **So an indication without an atlas does not get a ranking with a warning
> attached. It gets `NOT_USABLE` and the reason.** A number that looks like an
> answer is worse than a refusal, because only one of them gets checked.

`Indication.atlas is None` therefore produces `ScreenUsability.NOT_USABLE`, the
API reports it, and no ranked list is served for it.

## 3. Both Stage 3 modes

| | |
| --- | --- |
| **Mode B** — no target supplied | screen the surface proteome, propose single or paired targets. Exercised on every run to date. |
| **Mode A** — target supplied | evaluate whether that target is suitable for CAR-T. Code path exists; **never run end to end.** |

Mode A is proven on **a target the platform would not itself have chosen**, so
the test exercises validation rather than agreement with its own ranking. A
target it already ranks first would pass for the wrong reason.

## 4. Rejection criteria

Fixed before the run.

| | criterion | trips if |
| --- | --- | --- |
| **M1** | Indications coexist | any indication-specific artifact digest changes when the other indication runs |
| **M2** | Shared sources are not duplicated | any of UniProt, GTEx, HPA, GENCODE, SAbDab, domains gains a per-indication copy |
| **M3** | No indication-specific constant survives | any module-level constant still names a single cohort, atlas, lineage or organ |
| **M4** | An atlas-less indication refuses | a screen without an atlas returns a ranking rather than `NOT_USABLE` |
| **M5** | Positive pin — the atlas-less reason is structural | the refusal does not name C2 as the missing discriminator |
| **M6** | Mode A runs end to end | a supplied target does not produce a verdict |
| **M7** | Mode A is not self-agreement | the Mode A pin is a target inside the platform's own top 20 |
| **M8** | Mode A and Mode B agree on the shared evidence | a target evaluated in Mode A reports a different risk or composite than the same target in Mode B |
| **M9** | PDAC is unchanged | the reference indication's ranking, pool or end state differs from the recorded reference |
| **M10** | A missing source is named, not inferred | an unavailable component is reported without saying which source is absent |

M9 is the one the whole exercise rests on: **two indications is only real if the
first still produces what it produced before.**
