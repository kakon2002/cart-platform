# Handover: running the platform and the dashboard

Everything needed to stand this up on a clean machine. Read §1 first — it says
what the platform will and will not tell you, and every figure below assumes it.

---

## 0. What you are receiving

A cancer-agnostic CAR-T design platform: a thirteen-stage pipeline, an HTTP API
in two shapes, a binder benchmark, and twenty-three verifier scripts that test
the platform against criteria fixed before each run.

**Two indications are configured**: pancreatic ductal adenocarcinoma, which
returns five candidate designs, and invasive breast carcinoma, which returns
none and says why. Both are correct outputs.

---

## 1. Three things to know before reading any number

**Affinity is UNKNOWN on all 422 retrieved binder candidates.** No connected
evidence release carries an affinity, dissociation-constant or free-energy
column. Any comparison that needs a predicted and an experimental value — fold
error, rank correlation — has no value on either side. Those are *structurally
uncomputable*, not failed.

**83 of 315 structure-route binders are annotated against a different protein
than the target they were retrieved for.** The binder search queries a
structural database by the target's accession, so a hit means the entry
*contains* the target, not that the antibody in it was raised against it. The
retrieval counts did not change; what they mean did. Every binder row now
carries a target-match verdict.

**Nothing here has been measured on cells.** Normal-tissue risk is derived from
expression data. Every validation figure reads `TO_BE_MEASURED`, and every
acceptance threshold reads `TO_BE_SET_BEFORE_THE_RUN`.

---

## 2. Setup

Python **3.13**. Three dependencies.

```
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt     # Windows
./.venv/bin/python -m pip install -r requirements.txt           # macOS, Linux
```

`requirements.txt` is `numpy`, `h5py` and `pydantic`. Everything else is the
standard library, including the HTTP server.

**Use the virtual environment's interpreter for everything.** A bare `python`
or `py -3.13` is the system interpreter and has none of the three.

---

## 3. Data provisioning

The cache is roughly **12 GB**, of which 11 GB is a single-cell atlas used at
build time and not at serve time. `bootstrap.py` provisions it three ways:

```
.venv\Scripts\python.exe bootstrap.py                    # report what is present
.venv\Scripts\python.exe bootstrap.py --from-archive PATH # unpack a cache you were given
.venv\Scripts\python.exe bootstrap.py --from-release      # download the prepared cache
.venv\Scripts\python.exe bootstrap.py --from-sources      # rebuild from origin, ~3 hours
```

**Start with `--from-archive`** if a cache archive came with this package: it is
the fastest path and it excludes the heavy build-time inputs. Run the bare
command first either way — it prints what is present, what is missing, and the
release each source is pinned to.

`--from-sources` re-fetches from the public sources and takes about three hours.
It is the fallback when no archive is available, not the normal path.

---

## 4. Running it

### The pipeline, end to end

```
.venv\Scripts\python.exe run_all.py
```

Runs all stages for both indications and writes to `reports/`. Roughly 22
minutes on a warm cache.

### The server

```
.venv\Scripts\python.exe serve.py --port 8080
```

Serves two shapes on the same port:

- `/v1/cart/...` — the nineteen endpoints the dashboard expects.
- `/projects/...` — the platform's own surface, which the verifiers use.

It binds `127.0.0.1` by default. **Keep it on loopback.** The server sends
open cross-origin headers so a dashboard opened from disk can reach it, which
is safe on a loopback binding and is not an origin policy for anything exposed
to a network.

### The verifiers

```
.venv\Scripts\python.exe verify_adapter.py      # the dashboard surface, 14 criteria
.venv\Scripts\python.exe verify_api.py          # the HTTP contract, 12
.venv\Scripts\python.exe verify_benchmark.py    # the binder benchmark, 19
.venv\Scripts\python.exe verify_package.py      # the candidate package, 9
```

Each prints its criteria and stops on the first that trips. A verifier that
reports `TRIPPED` is telling you something is wrong, and the message says what.

---

## 5. Pointing the dashboard at the server

1. Start the server on port 8080 as above.
2. Open `specs/dashboard.html` in a browser — double-clicking the file works.
3. Open its settings drawer. **API base URL** defaults to
   `http://localhost:8080`, which is already correct if you used the port
   above. Change it if you started the server elsewhere, then save — it is
   remembered in the browser.
4. Use the console tab to call an endpoint.

**The order the endpoints expect:**

```
POST /v1/cart/projects        {"indication": "pancreatic ductal adenocarcinoma",
                               "target_mode": "DISCOVER"}
POST /v1/cart/runs            {"project_id": "<what the first call returned>"}
GET  /v1/cart/runs/{run_id}   until status is COMPLETE, about 20 minutes cold
POST /v1/cart/candidates/rank {"project_id": "..."}
```

**If a project call returns 400, read the body.** The platform refuses fields it
does not honour rather than accepting and ignoring them, and the refusal names
the field and says what to remove. `objective` and `delivery_mode` are refused
for that reason: nothing in the platform reads either, so accepting them would
tell you an instruction was followed when it was discarded. The dashboard's own
example request sends both.

**What is honoured:** `indication` (or `cancer_type`), `target_mode`,
`project_id` as your own reference, `max_final_candidates`, `architecture_mode`
where it maps, and `binder_mode: RETRIEVAL_FIRST`.

**Endpoints that answer with a refusal rather than a result:**
`structure/evaluate` returns `UNKNOWN`, `function/predict` returns
`INSUFFICIENT_VALIDATED_DATA`, `learning/calibrate` returns
`BLOCKED_INSUFFICIENT_DATA`, and `binders/rank`, `binders/optimize` and
`experiments/results` return `NOT_IMPLEMENTED`. Each carries its reason and
what would change it. They are complete answers, not stubs.

---

## 6. The benchmark

```
.venv\Scripts\python.exe run_benchmark.py       # produce and freeze the blind output
.venv\Scripts\python.exe compare_benchmark.py   # compare against the literature panel
```

The blind rule is enforced in four layers, not by intention: the algorithm's
import graph is walked and refused if it reaches the answers, the held-out
values are scanned for in every reachable module, the output is fingerprinted
and frozen, and the comparison refuses to run unless the frozen output was
committed *before* the answers were. The commit graph is the evidence.

Current result: SS1 and 15B6 were both recovered. M5 and M11 are named by no
mesothelin record in either connected source, so there was nothing to retrieve
them by — source coverage, not a retrieval failure.

---

## 7. Where things are

| | |
| --- | --- |
| `car_pipeline/stages/` | the thirteen stages and the scoring frame |
| `car_pipeline/api/` | the server, the `/v1/cart` adapter, the pipeline runner |
| `specs/` | every specification and design decision, including the two client documents |
| `specs/decision-binder-count-objective.md` | **an open decision that needs an answer** |
| `reports/packages/` | one document per candidate design |
| `benchmark/` | the blind harness and the frozen output |
| `data/` | the provisioned cache, from `bootstrap.py` |

---

## 8. Known limits, stated plainly

- **No structural stage.** Nothing predicts geometry. The cache holds no
  coordinate file, and no folding model fits ordinary hardware.
- **No functional stage.** Its inputs are experimental measurements that are not
  connected. This is the one gap a decision alone cannot close.
- **No binder generation.** Binder discovery retrieves and stops.
- **No cross-reactivity screen.** The proteome cache carries no sequence column
  and no alignment library is installed. The cost of connecting it is recorded
  in the gaps section of every candidate package.
- **Every surviving pancreatic design is an adaptor design**, which is two
  biologics rather than one. The second is an antibody the platform never names.
- **The anti-tag binder is murine and is emitted as deposited**, crystallisation
  artifacts included. It is not the molecule to order.

Every one of these is reported in the platform's own output, per candidate, in
the section each package titles *What this package cannot tell you*.
