# Deployment

Cloud Run, one instance. The container is built and verified; the deploy needs
a Google account this machine does not have — see **Status** at the end.

> ## The constraint to read first
>
> **This deploys one indication.** The platform's *design* is cancer-agnostic;
> this deployment is not, and the difference is a cache fingerprint.
>
> `SingleCellSource.malignant_entry()` fingerprints its artifact with a digest
> of the gene set. A second cancer type produces a different pool, a different
> digest, a cache miss, and the 8.3 GB single-cell matrix back in the serving
> path — the file this whole deployment shape exists by leaving out.
>
> That single line is the boundary between **half a day and a week**. It is not
> a tuning problem or a bigger machine; it is a different architecture, because
> reading an HDF5 file's chunked random access through a bucket mount is the
> wrong answer and the right one is a batch job that materialises the summary
> per indication ahead of time.
>
> Stated here rather than left to be discovered by whoever adds the second
> indication and finds the container will not start.

## The cache is 648 MB, not 9 GB

The number that shapes every option here is smaller than it looks, and it is
worth establishing before choosing anything.

```
data/                                       11 GB
  singlecell/totaldata-final-toshare.h5ad   8.3 GB    the expanded matrix
  singlecell/GSE202051_...h5ad.gz           2.6 GB    the archive it came from
  everything else                           648 MB
```

`SingleCellSource.load()` does not open the matrix. It calls
`build_group_means()`, which goes through `cache.ensure()` and returns the
cached `group_means.npz` — 5.7 MB — whenever the fingerprint matches. The same
holds for `load_malignant()`. **The 8.3 GB file is a build-time input, not a
serve-time one.** A deployed service that never invalidates those fingerprints
never opens it.

Measured rather than reasoned: with `totaldata-final-toshare.h5ad` renamed out
of the way, a full screen completed in **7.3 seconds** — all eight stages,
3,466 ranked, a pool of 200, and the same `NO_DESIGN_REACHES_THE_END`. Nothing
in the serving path touched it.

The fingerprint that matters is in `malignant_entry()`: it includes a digest of
the gene set. **A different indication produces a different pool, a different
digest, a cache miss, and an 8.3 GB stream.** So the matrix has to live
somewhere reachable — just not somewhere fast, and not in the serving path.

That splits cleanly:

| | where | why |
| --- | --- | --- |
| 648 MB derived + raw sources | in the container image | read on every run, chunked random access |
| 8.3 GB matrix | a bucket | read only when a fingerprint changes |
| 2.6 GB archive | nowhere | it is the download the 8.3 GB was expanded from |

## Shortest path: Cloud Run, one instance

Managed HTTPS and a public URL with no load balancer, no TLS certificate and no
machine to patch. That is what makes it shortest — not the compute.

```bash
gcloud run deploy cart-platform --source . \
  --region=us-central1 \
  --min-instances=1 \
  --max-instances=1 \
  --concurrency=1 \
  --no-cpu-throttling \
  --cpu=4 --memory=8Gi \
  --timeout=900 \
  --allow-unauthenticated
```

No inline comments in that block on purpose: in bash a `\` followed by a space
escapes the space, not the newline, so a trailing `# ...` ends the command
there. A version of this with the reasons written beside the flags would paste
as a deploy missing `--no-cpu-throttling`, `--cpu`, `--memory` and
`--allow-unauthenticated` — Cloud Run defaults of 1 CPU and 512 MiB, which
OOM-kills on the 413 MB DepMap CSV, with the background thread throttled. The
reasons are below instead, and `deploy.ps1` carries them as real comments.

**`--no-cpu-throttling` is the one that breaks it silently if omitted.** It
produces no error: a service that looks deployed and whose jobs never advance.

Every flag on that command is load-bearing:

- **`--max-instances=1`.** The job table is a dict in memory. With two
  instances a client can poll the one that is not running its job and be told
  the run does not exist.
- **`--concurrency=1`.** Necessary, and easy to miss. The application's guard
  is **per project** — `start_run` rejects a second run only for the *same*
  project id — so two projects submitted to one instance run two pipelines at
  once and both write the shared binder cache, whose writer unlinks the
  manifest before rewriting the payload. Cloud Run's default concurrency is 80,
  so pinning instances to one does nothing about this on its own. Making the
  guard global in `start_run` would be the better fix and would let this flag
  go; until then the flag is what serialises runs.
- **`--min-instances=1`.** A scale-to-zero instance loses the job table between
  the run and the poll. It also avoids pulling a ~1.3 GB image on a cold start.
- **`--no-cpu-throttling` — because the screen runs on a background thread
  after the `202`.** Cloud Run's default allocates CPU only while a request is
  in flight. `start_run` returns `202` immediately and the pipeline continues on
  a worker thread, so under the default that thread is throttled to near zero
  the instant the response is written. The symptom is not an error: the job
  sits at `running`, the stage never advances, and polling returns a valid
  answer forever. **Confirm this flag's current name at deploy time** — it has
  been spelled more than one way across releases.
- **`--concurrency=1` — a second line of defence, and not the fix.** The run
  guard *was* per project: `start_run` rejected a second run only for the same
  project id, so two projects would each start a pipeline and both would write
  the shared binder cache, whose writer replaces the payload before the
  manifest. **A concurrency cap does not fix that**, which is worth stating
  because it is the obvious assumption: a run is a detached thread and the
  request returns `202` in milliseconds, so a limit on *in-flight requests*
  never sees two runs overlap. The guard in `start_run` is now global — it
  rejects a run while any other is queued or running, and says which project
  holds it. This flag now only bounds polling load.
- **`--allow-unauthenticated`.** Open, so the URL works from a browser or a
  bare `curl` with nothing to install. That is a deliberate trade for a demo
  endpoint: it also means anyone with the URL can start a screen. Close it with
  `--no-allow-unauthenticated` once the demo is done, and callers then need
  `-H "Authorization: Bearer $(gcloud auth print-identity-token)"`.

### The container

Both done and verified.

1. `serve()` binds `0.0.0.0` and reads `PORT`. `--host`/`--port` override, with
   `HOST`/`PORT` environment fallbacks resolved *after* parsing, so an empty or
   malformed `PORT` gives a named error instead of a traceback and never blocks
   an explicit `--port`.
2. `Dockerfile` and `.dockerignore`. The ignore file is what keeps the image at
   ~1.1 GB rather than ~12 GB, and it is the only thing standing between this
   deployment and the 8.3 GB matrix.

The image ships `data/stage5` warm — the retrieved binder set — so the first
screen can skip one network call per pool member. On the build machine that
takes a full screen to **9.1 seconds** instead of six minutes.

Verified without Docker by assembling exactly the file set the `Dockerfile`
copies — `.dockerignore` applied, 86 files, 680 MB, neither single-cell archive
present — and running the server from that tree: project created, run submitted,
job complete in 9.1 s, and `/targets`, `/pairs`, `/constructs` and `/result` all
answering `200` with their named statuses.

**That 9.1 s is a Windows number and may not survive the move to Linux.**
`load_or_retrieve` reuses the cache only when the manifest's Stage 4 hash
matches, and that hash chains from `stage3.configuration_hash(...)`, whose
payload embeds calibration floats produced by numpy reductions. A different
BLAS build or SIMD dispatch can change a last bit, change the JSON, change the
SHA-256, and miss. **The miss is safe but slow**: the first screen falls back to
the real five-minute retrieval, which is the honest cost and not a failure.
Confirm the first deployed run's timing rather than quoting 9 s at anyone.

The pool digest is the same class of risk with a worse landing, and that is
what `CART_NO_MATRIX_FETCH=1` in the `Dockerfile` is for: a digest miss now
raises a named error naming this document, instead of downloading 2.6 GB and
expanding it to 8.3 GB onto a filesystem that is in memory here.

### Rough timings

| | |
| --- | --- |
| Dockerfile and `.dockerignore` | 1 h |
| First build and push of a 1.3 GB image | 30–45 min, mostly upload |
| Deploy, flags, and an IAM invoker | 1 h |
| Point `verify_api.py` at the URL and run it | 30 min |
| Upload the matrix to a bucket | 1 h, mostly upload |

**Half a day**, assuming a project, billing, and the Cloud Run and Artifact
Registry APIs already enabled. A day if any of those need setting up.

Cost is roughly **$80–110/month** for an always-on 4 vCPU / 8 GiB instance,
plus about $0.20/month of bucket storage. Confirm against current pricing;
these are the right order of magnitude, not a quote.

## What this shape cannot do, and what fixing it costs

The single instance is the whole compromise. It buys a half-day deployment and
it caps the system at one screen at a time.

Each of these breaks the moment you need more:

- **Two concurrent screens.** The job table has to leave memory — Firestore or
  Memorystore — and the run has to leave the request path, as a Cloud Run Job
  triggered through Cloud Tasks or Pub/Sub. The HTTP surface then only reads
  state.
- **A second indication.** The gene-pool digest changes, so the 8.3 GB matrix
  is back in the serving path. Reading an HDF5 file's chunked random access
  through a bucket mount is slow enough to be the wrong answer; the right one is
  a batch job that materialises the `.npz` per indication ahead of time and
  publishes it to the bucket, leaving the service reading only small derived
  artifacts. This is a real piece of work, not a flag.
- **Surviving a restart mid-run.** Today a restart loses queued and running
  jobs and loses nothing else, because every stage's real output is on disk
  under its own manifest. Making a run resumable means persisting job state,
  which is the same change as the first item.

**About a week**, and it is the same week for all three, because they are one
change: move the job table and the run out of the process.

## The alternative, and why not

A Compute Engine VM with a persistent disk is conceptually simpler and the
8.3 GB matrix could sit on the disk with proper block access. But reaching it
over the internet means a static IP, a firewall rule, a systemd unit, and then
either a load balancer or a domain plus a certificate — and afterwards you own
OS patching. Cloud Run gives HTTPS and a URL as a property of deploying. The VM
is the better answer once the matrix must be in the serving path, which is the
second-indication case above, and that rewrite is the trigger to revisit it.

## Cancer type to final result, in eight calls

Verified end to end against the container's own file set. Substitute the
deployed URL for `BASE`. Works in bash, macOS, Linux and Cloud Shell; Windows
PowerShell needs different quoting, so prefer Git Bash there.

```bash
BASE=https://cart-platform-XXXXXXXX.us-central1.run.app

# 1. A project, from a cancer type. No target is named — the null is what
#    selects discovery mode.
curl -s -X POST "$BASE/projects" \
  -H 'Content-Type: application/json' \
  -d '{"cancer_type":"Pancreatic Ductal Adenocarcinoma"}'
PID=<project_id from above>

# 2. Start the screen. Returns 202 immediately; it does not wait.
curl -s -X POST "$BASE/projects/$PID/runs"
JID=<job_id from above>

# 3. Poll. Reports the stage it is on, and finishes at "complete".
curl -s "$BASE/jobs/$JID"

# 4-8. Once complete:
curl -s "$BASE/projects/$PID/targets?limit=10"   # ranked, with breakdowns
curl -s "$BASE/projects/$PID/pairs?limit=10"     # co-expression + span percentile
curl -s "$BASE/projects/$PID/constructs"         # NO_BUILDABLE_CONSTRUCT + reasons
curl -s "$BASE/projects/$PID/result"             # the end state + attrition chain
curl -s "$BASE/projects/$PID/evidence/MSLN"      # every stage's view of one gene
```

To poll without reading the output each time. **Both terminal states**, not
just the good one — matching only `complete` spins forever on a failed job
while the answer sits in the body it is discarding:

```bash
until curl -s "$BASE/jobs/$JID" | grep -qE '"status": "(complete|failed)"'; do
  sleep 5
done
curl -s "$BASE/jobs/$JID"
```

**What the last two return is the point.** `/constructs` answers `200` with
`NO_BUILDABLE_CONSTRUCT` and the reasons computed from the run; `/result`
answers `200` with `NO_DESIGN_REACHES_THE_END` and all 200 candidates
attributed to the first gate each one failed. Neither is an error and neither
is an empty list.

## Status

| | |
| --- | --- |
| `Dockerfile` and `.dockerignore` | written; file set verified by running the server from exactly what they produce |
| `serve()` binds `0.0.0.0`, honours `PORT` | done |
| Both load-bearing flags, with reasons | in `deploy.ps1` and above |
| The curl sequence | verified end to end |
| Google Cloud SDK | installed (581.0.0) |
| **`gcloud auth login`** | **blocked — needs a browser sign-in to an account this session does not have** |
| `gcloud run deploy` | blocked behind the above |

The sign-in is the wall, and it is not one that can be worked around: it is a
Google account, its project, and its billing. After

```
gcloud auth login
gcloud config set project <PROJECT_ID>
```

`.\deploy.ps1` does the rest — enables the three APIs, deploys with the flags
above, then creates a project, submits a run, polls it and prints the attrition
chain, so a deploy that comes up wrong says so rather than looking finished.

## Stages 7 and 8

They stay untouched. Nothing in this path needs them: the deployment serves what
the nine implemented stages produce, and the end state it serves —
`NO_DESIGN_REACHES_THE_END` — is reached at Stage 3 for 199 of 200 targets and
Stage 5 for the last one. Neither 7 nor 8 is on that path.
