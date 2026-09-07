# Runbook

From an empty directory to a running server with the dashboard connected. One
ordered sequence, with what each step should print and how long it should take.

Every timing here was measured from an empty directory, not estimated. Where a
step is untested, it says so.

If something looks wrong, check §7 before assuming it is. **A correct first run
prints several things that read like failure.**

---

## 0. Before you start

You need Python **3.13** and about **2 GB of disk** for the cache the normal
path provisions. Network access is required once, for provisioning.

You do not need a GPU, a database, or any service running.

---

## 1. Unpack and check the interpreter

```
.venv\Scripts\python.exe -c "import sys; print(sys.version); print(sys.prefix)"
```

Run this *after* step 2. It is listed first because everything below depends on
it and because a wrong interpreter fails much later, somewhere unrelated.

**Expected:** a `3.13.x` version, and a `sys.prefix` ending in `.venv` inside
the directory you unpacked into.

**Do not use a bare `python`.** What it resolves to depends on the machine —
on at least one it was already a different virtual environment.

---

## 2. Create the environment  ·  ~6 seconds

```
python -m venv .venv
```

**Expected:** no output. `.venv\Scripts\python.exe` exists afterwards.

The bare `python` here is the only place one is used, and only to create the
environment. If your system `python` is not 3.13, invoke the 3.13 you have —
the venv takes its version from whatever creates it.

---

## 3. Install the three dependencies  ·  ~10 seconds

```
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

**Expected, last line:**

```
Successfully installed annotated-types-0.8.0 h5py-3.16.0 numpy-2.5.2
pydantic-2.13.4 pydantic-core-2.46.4 typing-extensions-4.16.0
typing-inspection-0.4.4
```

Three direct dependencies, four transitive. Everything else is the standard
library, including the HTTP server.

Now go back and run §1.

---

## 4. See what the cache holds  ·  ~1 second

```
.venv\Scripts\python.exe bootstrap.py
```

**Expected on a first run: thirteen `MISSING` lines.** This is correct. The
archive ships no data. The report ends:

```
  1/8 shared sources usable
  6 per-indication cache(s) missing:
    ...
Missing: uniprot, hpa, gtex, depmap, genespan, antibodies, domains, ...
Run with --from-release (minutes) or --from-sources (~3 hours).
```

Read the per-source rebuild costs it prints; they are the real ones.

---

## 5. Provision the cache  ·  ~106 seconds

```
.venv\Scripts\python.exe bootstrap.py --from-release
```

**Expected, last three lines:**

```
  8/8 shared sources usable
  the 8.3 GB matrix and its 2.6 GB archive are build-time only and are not expected here
The cache is complete. Nothing to do.
```

Every source should read `present`, and both indications should read
`tcga=ok  singlecell=ok  depmap=ok`.

**No cache archive ships with this package**, so `--from-archive` is not the
path here — it exists for a cache someone hands you directly. `--from-sources`
rebuilds everything from the public sources in about three hours and is the
fallback when nothing else is reachable.

### If it does not complete

**This step downloads about 300 MB and can fail transiently.** It failed once
in three measured attempts, running for 89 seconds and leaving the cache
empty. **Re-run the same command.** It is safe to repeat: an unfinished file
is written beside the real one and never renamed over it, so a failed attempt
cannot corrupt a good cache.

A failed attempt does leave an unfinished download behind. `bootstrap.py`
names it at the top of its report:

```
  1 unfinished download(s), left by an attempt that did not complete:
    data/uniprot/human_reviewed.tsv.partial  11.5 MB
    They are not counted in the sizes below and are never packaged.
```

That is a note, not a problem. Re-running `--from-release` replaces it.

`bootstrap.py` provisions from release tag **`data-v2`**. If your copy is
configured for `data-v1`, it is an older checkout: that payload is missing
eight malignant-cell digests and the anti-tag binder, and §7 describes what
that costs. `data-v1` is deliberately left published so a cache provisioned
from it can still be identified.

**One thing this step does not do, if you are on `data-v1`.** It provisions
the caches, and the first run afterwards still rebuilds one single-cell
artifact from raw at a cost of about 40 minutes. §7 explains why and how to
tell it from a hang. On `data-v2` this does not happen.

**Do not proceed past this step until the report reads `8/8 shared sources
usable`.** Everything downstream depends on the cache, and a half-provisioned
one produces failures that look unrelated to provisioning — a run that reaches
`FAILED` in under two minutes, and a `candidates/rank` that answers 500.

---

## 6. Start the server  ·  ~2 seconds

```
.venv\Scripts\python.exe serve.py --port 8080
```

**Expected:**

```
  listening on http://127.0.0.1:8080
```

It binds loopback. Keep it there: the server sends open cross-origin headers so
a page opened from disk can reach it, which is safe on loopback and is not an
origin policy for anything exposed to a network.

Leave this running and use another terminal for what follows.

---

## 7. What a correct first run looks like when it looks wrong

**Read this before concluding anything is broken.** Five things a correct
system prints that read like failure.

### Thirteen `MISSING` lines from `bootstrap.py`

Correct on a first run. The archive carries code, specifications and reports,
and no data. §5 fixes it in under two minutes.

### Three criteria trip, by name, and they are not failures

| id | stage | what it means |
| --- | --- | --- |
| `R14` | 3, target discovery | The staining arm gates on presence rather than amount: every grade blocks in criticality tiers 1 and 2, none reaches the ceiling in tier 3. **An open tolerance decision nobody has taken**, priced in `reports/staining-veto-decision.md`. Not a regression. |
| `P17` | 4, target pairing | 22 of 30 dual recommendations had exactly one admissible partner, so the stage named a partner it did not choose between. **Measured deliberately** rather than hidden, after `P13` showed the concentration is structural. |
| `K2` | 6, construct assembly | No dual carries a binder on both arms, so the two-arm join is not exercised anywhere in this decision set. The criterion **reports that rather than clearing on an empty set** — which is the behaviour you want. |

Five further criteria trip in stage 4 and 4a — `P4`, `P8`, `P13`, `P15`, `A6` —
and are **accepted exemptions**. `run_all.py` prints the recorded reason beside
each. Anything *not* in this list is a genuine regression and worth stopping
for.

### The breast indication returns zero designs

Correct. Invasive breast carcinoma is configured and screens cleanly, and no
design survives its gates. **Zero candidates with reasons is a result**, not a
crash. Pancreatic returns five. Both are correct outputs, and breast is the one
that exercises the refusal paths.

### `BM5` says it cannot ask its question

From an unpacked archive, `verify_benchmark.py` prints something like:

```
clear  BM5: no repository here, so the commit ordering cannot be checked from
       this copy...
```

Correct. That criterion reads the commit graph to prove the benchmark's blind
output was frozen *before* the literature answers were revealed. An archive
carries no history, so the question cannot be asked from a copy — which is
different from the rule being broken. It is enforced where the work is done.
What travels with the archive is the frozen fingerprint, not the history that
ordered it.

### A first pipeline run that takes about 40 minutes

**This applies to caches provisioned from `data-v1`. On `data-v2` it does not
happen.** If you provisioned from the older tag, the first thing that runs the
pipeline — a verifier, `run_all.py`, or a run started over HTTP — downloads
2.6 GB and expands it to 8.3 GB before doing anything else. Measured at 2,003
seconds on one pass, and it happens for both indications.

It is not a hang and it is not the verifier. Here is the mechanism, because it
is not guessable from the symptom.

The single-cell malignant-cell cache is keyed by a **digest of the exact gene
set** a run asks for. The published release carries digests built from an older
pool, so a current run's lookup misses, and the source rebuilds the artifact
from raw counts — which needs the archive the same `bootstrap.py` report calls
*build-time only and not expected here*. The file being regenerated is about
**2.5 MB**.

**How to tell it apart from a hang:**

```
dir data\singlecell\*.partial          # Windows
ls -la data/singlecell/*.partial       # macOS, Linux
```

A `.partial` that grows between two checks means it is downloading. Nothing
else in this platform writes one.

**Why it happened and why it will not recur.** The packager now computes the
digest a standard run of each registered indication asks for and refuses to
build an archive without it, so a stale release cannot be produced by accident
again. `data-v2` was built under that check and carries both. A cache that has
already paid the 40 minutes holds the artifact and later runs are minutes.

> **Retraction.** An earlier version of this section said the wait was because
> the verifier "runs the whole pipeline through the adapter". That was written
> without being established and is wrong twice: it named the wrong cause, and
> it told a reader to expect a delay that is a property of one stale artifact
> rather than of the code.

---

## 8. Verify the installation

Two options, depending on how much time you have.

### Quick — one verifier, seconds

```
.venv\Scripts\python.exe verify_package.py
```

**Expected:** `9/9 criteria clear`, then a summary of five candidate packages
and the elements they cannot carry.

About 8 seconds once the pipeline has run before. **The first time, expect
about 40 minutes** — see §7, which explains why and how to confirm it is
working rather than stuck.

### Whole suite — 16 stages, 27 minutes

```
.venv\Scripts\python.exe run_all.py --fresh
```

Two measured runs took 28.1 and 26.8 minutes, with derived artifacts cleared
each time. Most of it is three stages: architecture routing, binder discovery
and the multi-indication check, at roughly 3 to 10 minutes each.

**Both figures are from a machine whose cache was already complete.** On a
freshly provisioned one, add the 40 minutes §7 describes, once.

**Expected:** every stage reports `clear` except stages 3, 4, 4a and 6, whose
trips are the eight named in §7. The last measured run was **217 of 225
criteria clear**, with an `unexpected` list containing exactly `3/R14`,
`4/P17`, `6/K2`, and:

```
  raw caches unchanged: 80 file(s) identical in size and modification time
  after the run (trials excluded: query-keyed, rebuilt when the screened set
  changes)
```

That line is the cache guard: every pinned payload is checked against its
manifest before the first stage, and a size-and-mtime fingerprint is compared
afterwards. `trials` is excluded because it is keyed by the screened antigen
set and rewriting it is correct.

The total and the full per-criterion detail are written to
`reports/full-run.md` and `reports/full-run.json`. **Read the total there
rather than from any number written in prose**, including this document —
a figure in prose goes stale the moment a criterion is added, and this project
has already reported one that did.

---

## 9. Drive it over HTTP

With the server from §6 running. Bodies exactly as given.

```
POST /v1/cart/projects   {"indication": "pancreatic ductal adenocarcinoma",
                          "target_mode": "DISCOVER"}
  -> 201  {"project_id": "...", "canonical_project_id": "...",
           "status": "CREATED", ...}

POST /v1/cart/runs       {"project_id": "<project_id from above>"}
  -> 202  {"run_id": "...", "status": "QUEUED"}

GET  /v1/cart/runs/{run_id}
  -> 200  {"status": "RUNNING", "stage": "...", "stages_complete": 4,
           "stages_total": 9, "progress": null}
     Poll until "COMPLETE". About 20 minutes on a cold cache.

POST /v1/cart/candidates/rank  {"project_id": "..."}
  -> 200  {"run_status": "RANKED_CANDIDATES", "eligible_candidate_count": 5, ...}
```

**`progress` is `null` on purpose.** A fraction would imply a measured
proportion of the work; what exists is an ordered stage list, and
`stages_complete` / `stages_total` carry that without inventing one.

**If a call returns 400, read the body.** The platform refuses fields it does
not honour rather than accepting and ignoring them, and the refusal names the
field. `objective` and `delivery_mode` are refused for that reason — nothing
reads either, so accepting them would tell you an instruction was followed when
it was discarded. The dashboard's own example request sends both.

---

## 10. Connect the dashboard

1. Start the server on port 8080 (§6).
2. Open `specs/dashboard.html` in a browser. Double-clicking the file works.
3. Open the settings drawer. **API base URL** defaults to
   `http://localhost:8080`, already correct if you used that port. Change it if
   you started the server elsewhere and save; it is remembered in the browser.
4. Use the console tab to issue the calls in §9.

**This step is untested.** Everything else in this runbook was executed from an
empty directory and timed. The dashboard step needs a browser and was not run,
so the instructions above are derived from reading the page's source rather
than from watching it work. The parts that *were* tested headlessly are that
the server answers on `/v1/cart/…`, that it sends
`Access-Control-Allow-Origin: *` so a page opened from disk is not blocked, and
that it answers the `OPTIONS` preflight a cross-origin POST sends first.

### Headless equivalent, if you have no browser

Also untested as a substitute for the page, though each call below was executed
directly:

```
.venv\Scripts\python.exe -c "import json,urllib.request; r=urllib.request.Request('http://127.0.0.1:8080/v1/cart/projects', data=json.dumps({'indication':'pancreatic ductal adenocarcinoma','target_mode':'DISCOVER'}).encode(), headers={'Content-Type':'application/json'}); print(urllib.request.urlopen(r).read().decode())"
```

Then follow §9 with the returned `project_id`.

---

## 11. When something is actually wrong

| symptom | likely cause |
| --- | --- |
| `ModuleNotFoundError: numpy` | the system interpreter, not the venv. Re-read §1 |
| `OSError: [WinError 10048]` on start | the port is in use. Pass a different `--port` and set the dashboard's API base to match |
| A verifier trips something not in §7 | a genuine regression. The message says what; the criterion id maps to a specification in `specs/` |
| A stage reads `no criteria, exit N` | the verifier crashed before reporting. Its full transcript is in `reports/run-logs/` |
| `bootstrap.py` reports `BROKEN` | a payload disagrees with its manifest. Re-run `--from-release`; the caches are pinned releases and a drifted one makes the run unreproducible |
| Provisioning ends without `The cache is complete` | a transient download failure. Re-run the same command; it is safe to repeat and cannot corrupt a good cache |
| A run reaches `FAILED` in under two minutes, or `candidates/rank` answers 500 | almost always an incomplete cache. Run `bootstrap.py` and check for `8/8 shared sources usable` before looking anywhere else |
| The run says `RAW CACHES WERE MODIFIED` | a stage wrote to a pinned cache. That is a defect worth reporting, not something to work around |

---

## 12. Outstanding

One thing in this document describes a state that is waiting on an action
outside the repository.

**`data-v2` is built and verified but not yet published.** Until it is, `gh
release download data-v2` fails and §5 cannot complete as written; provision
from `data-v1` and expect what §7 describes. The payload and its checksum are
built by `bootstrap.py --package`, which refuses to produce a stale one.

What `data-v2` adds over `data-v1`, measured file by file: **eighteen files,
21.6 MB**, and nothing removed. Sixteen are malignant-cell artifacts covering
eight gene-set digests, including the one each indication's standard run asks
for. Two are the anti-tag binder, which `data-v1` omits entirely — it is
fetched on demand from a public structure database when absent, so its absence
is slow rather than fatal. Five derived artifacts differ in size because they
are rebuilt per run.

---

## 13. Where to go next

- `ORIENTATION.md` — what the platform is, what it refuses to do, and the
  failure shapes this project keeps producing. Read it before changing
  anything.
- `HANDOVER.md` — the same operational ground as this document, arranged by
  topic rather than sequence.
- `specs/decision-binder-count-objective.md` — an open decision that changes
  which design the platform advances.
