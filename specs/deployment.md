# Deployment

Not built. This is what it would take.

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

```
gcloud run deploy cart-platform --source . \
  --region=us-central1 \
  --min-instances=1 --max-instances=1 --concurrency=1 \
  --no-cpu-throttling \
  --cpu=4 --memory=8Gi \
  --no-allow-unauthenticated
```

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
- **`--no-cpu-throttling`.** Cloud Run's default allocates CPU only while a
  request is in flight. The screen runs on a background thread after the `202`
  has been sent, so under the default it would be throttled to near zero the
  moment the response is written. **Confirm this flag's current name at deploy
  time** — it has been spelled more than one way across releases.
- **`--no-allow-unauthenticated`.** An unauthenticated endpoint that spawns a
  multi-minute screen per call is a way to spend money on strangers.

### What has to change in the code

Two lines, and one is already done.

1. `serve()` must bind `0.0.0.0` and read `PORT`. **Done** — `--host` and
   `--port` with `HOST`/`PORT` environment fallbacks.
2. A `Dockerfile`: `python:3.13-slim`, `pip install -r requirements.txt`, copy
   the package and the 648 MB of `data/` minus the two single-cell archives,
   then `CMD ["python", "-m", "car_pipeline.api.server", "--host", "0.0.0.0"]`.
   A `.dockerignore` that excludes `*.h5ad`, `*.h5ad.gz` and `.venv` is what
   keeps the image at 1.3 GB rather than 12 GB.

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

## Stages 7 and 8

They stay untouched. Nothing in this path needs them: the deployment serves what
the nine implemented stages produce, and the end state it serves —
`NO_DESIGN_REACHES_THE_END` — is reached at Stage 3 for 199 of 200 targets and
Stage 5 for the last one. Neither 7 nor 8 is on that path.
