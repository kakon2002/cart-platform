# Deploying this platform

A cancer-agnostic CAR-T target discovery pipeline, served over HTTP. This gets
it from a fresh clone to a URL. No prior knowledge of the project assumed.

Budget about **35 minutes**, most of it waiting on a container build.

---

## 0. What you are deploying, in one paragraph

You give it a cancer type. It screens the whole human surface proteome for
targets, pairs them, routes each to a receptor architecture, retrieves antibody
binders, assembles CAR constructs, and ranks what survives.

Two indications are configured — pancreatic ductal adenocarcinoma and invasive
breast carcinoma — and the deployment holds both without rebuilding. **They do
not return the same shape of answer, and both are correct.**

For PDAC, 192 of 200 candidates are blocked on normal-tissue risk and 3 more
have no retrievable binder; **five designs reach the end**, all `BUILDABLE` and
all carrying a full sequence. Every one routes to an adaptor receptor, which
binds a tag rather than the antigen, with a separately dosed adaptor molecule
carrying the specificity. That is what admits them: an adaptor's exposure can be
stopped, so it answers to a looser risk ceiling than a receptor whose exposure
cannot be withdrawn.

For breast, **no design reaches the end**. The constructs endpoint answers
`NO_BUILDABLE_CONSTRUCT` at HTTP 200 with the reasons behind it, and the end
state is `NO_DESIGN_REACHES_THE_END`.

The API reports both as results with reasons, not as errors and not as empty
lists. `RANKED` and `NO_DESIGN_REACHES_THE_END` are both the platform working.

---

## 1. Prerequisites

| | |
| --- | --- |
| **Google Cloud account** | with a project and **billing enabled** — Cloud Run will not deploy without it |
| **Google Cloud CLI** | https://cloud.google.com/sdk/docs/install |
| **GitHub CLI** | https://cli.github.com — only to fetch the data cache, which lives in a private release |
| **Python 3** | any version, to run `bootstrap.py`. Deploying needs nothing more; **3.13 with the dependencies** only if you rebuild the cache from origin or run the pipeline locally |
| **Disk** | ~1.5 GB |

You do **not** need Docker. The image is built in the cloud from source.

> **Cost.** This deploys one always-on instance at 4 vCPU / 8 GiB, roughly
> **$80–110/month** while it is up. It is not covered by the always-free tier,
> which requires scaling to zero — and this cannot scale to zero, because the
> job table lives in memory and a poll must reach the instance running the job.
> Delete the service when you are finished:
> `gcloud run services delete cart-platform --region us-central1`

---

## 2. Clone, and get the data

**`data/` is not in git.** It is 852 MB of cached scientific sources — the human
proteome, a tumour cohort, expression atlases, a single-cell atlas, antibody
structures. A clone has none of it, and nothing runs without it.

```bash
git clone https://github.com/kakon2002/cart-platform.git
cd cart-platform
python bootstrap.py --from-release
```

That downloads a 471 MB archive from the repository's private release, verifies
its checksum, and unpacks it. **Two minutes.** It refuses to unpack on a
checksum mismatch rather than proceeding — a truncated transfer would give you a
cache that reads as present and answers with the wrong data.

Check what you have at any point:

```bash
python bootstrap.py
```

Expect `8/8 shared sources usable` and, beneath it, a per-indication line
reading `tcga=ok  singlecell=ok  depmap=ok` for **every** registered indication.

**The per-indication block is the one to read.** The shared sources describe the
human body and there is one copy of each; the cohort, the atlas and the
dependency matrix describe a tumour and there is one *per indication*. A single
averaged count would report a clone as ready while a second indication had
nothing — which is exactly what it used to do.

Two entries look unusual and are correct:

- **singlecell** shows ~18 MB, not 11 GB. The 8.3 GB matrix and the 2.6 GB
  archive it came from are build-time inputs; the derived summaries are what a
  served run reads.
- **trials** shows `deferred`. Its cache is keyed by the screened antigen list,
  so it has no meaning until a pool exists, and it is built during the first run.
- **singlecell** carries only derived summaries. The atlases themselves — 8.3 GB
  for one indication, 844 MB for the other — are build-time inputs and are
  deliberately not in the archive.

A payload named by a manifest but absent is reported `BROKEN` with the filename,
rather than counted as present.

<details>
<summary><b>If you would rather build the cache from the original sources</b></summary>

Everything is fetchable programmatically. No accounts, no registration, no
manual downloads — UniProt, the Human Protein Atlas, GTEx, the GDC, DepMap via
figshare, GENCODE, SAbDab, NCBI GEO, RCSB, ClinicalTrials.gov.

```bash
.venv\Scripts\python.exe bootstrap.py --from-sources    # Windows
.venv/bin/python bootstrap.py --from-sources            # macOS, Linux
```

**This path alone needs the project interpreter.** It imports the pipeline,
which needs `h5py` and `numpy`, so a bare `python` fails on the first import
before a byte is fetched. Create the environment first with `python -m venv
.venv` and `pip install -r requirements.txt`. Every other `bootstrap.py` command
is standard library only and runs under any Python 3.

**Allow about three hours and 12 GB of free disk.** Nearly all of it is one
step: deriving the single-cell group means streams a 8.3 GB matrix end to end
and took **2 h 19 min** on the machine that first built this cache. Every
download combined was about fifteen minutes.

Interrupting is safe — each artifact gets its manifest only once it is complete,
so a re-run resumes rather than restarts.

The deployed container does **not** carry that 8.3 GB matrix. It is a build-time
input; the 5.7 MB summary derived from it is what the served pipeline reads.
Verified by renaming the matrix away and completing a full screen in 7.3 s.
</details>

---

## 3. Sign in and pick your project

```bash
gcloud auth login
```

Opens a browser. Pick your account, **Allow**. Ends with
`You are now logged in as [you@example.com]`.

```bash
gcloud projects list
```

**If it is empty**, create one — the ID must be globally unique:

```bash
gcloud projects create cart-platform-001 --name="CAR-T Platform"
```

**Enable billing** — this part is console-only:

1. https://console.cloud.google.com/billing — add a billing account if you have
   none (needs a card; new accounts usually carry $300 of free credit)
2. `https://console.cloud.google.com/billing/linkedaccount?project=YOUR_PROJECT_ID`
3. **Link a billing account** → choose yours → **Set account**

Then point the CLI at it. **This is the one place your project ID goes** — every
script reads it from here:

```bash
gcloud config set project YOUR_PROJECT_ID
```

Confirm. This deploys nothing and bills nothing:

```bash
./deploy.sh --check-only          # macOS, Linux, WSL, Git Bash, Cloud Shell
.\deploy.ps1 -CheckOnly           # Windows PowerShell
```

You want:

```
signed in as you@example.com
project your-project-id
billing enabled

Ready to deploy. Re-run without --check-only.
```

Anything else stops here and names the fix. A billing status that cannot be
*read* also stops — unknown is not the same as enabled, and the alternative is
discovering it eight minutes into a build. That case is genuinely ambiguous:
reading it needs a permission on the *billing account* that a project Owner
often does not hold. If you know billing is on:

```bash
SKIP_BILLING_CHECK=1 ./deploy.sh      # bash
.\deploy.ps1 -SkipBillingCheck        # PowerShell
```

There is no override for billing that is genuinely *off*. That one is a wall.

---

## 4. Deploy

```bash
./deploy.sh                       # macOS, Linux, WSL, Git Bash, Cloud Shell
.\deploy.ps1                      # Windows PowerShell
```

Identical behaviour; pick whichever matches your shell. Each one enables three
APIs (`run`, `cloudbuild`, `artifactregistry`), builds the image in the cloud,
deploys it, then runs a smoke test — creating a project, submitting a run,
polling it, and printing the result. **A deploy that comes up wrong says so
rather than looking finished.**

**Roughly 10–15 minutes**, nearly all of it the first container build.

### The flags, and why

Two of them are load-bearing and produce no error when wrong:

- **`--no-cpu-throttling`** — the screen runs on a background thread *after* the
  HTTP response is sent. Cloud Run's default allocates CPU only while a request
  is in flight, which would throttle that thread to near zero. The symptom is
  not a crash: the job sits at `running` and the stage never advances.
- **`--min-instances=1 --max-instances=1`** — the job table is a dictionary in
  memory. Scale to zero and it is gone between run and poll; scale out and a
  poll can reach an instance that knows nothing about the job.

`--concurrency=1` is also set, as a second line of defence only. It does not
serialise runs — a run is a detached thread and the request returns in
milliseconds, so a request-concurrency cap never sees two runs overlap. What
actually serialises them is a global guard in the application.

### The URL

```
https://cart-platform-<hash>-uc.a.run.app
```

The script prints it. Retrieve it again any time with:

```bash
gcloud run services describe cart-platform --region us-central1 \
    --format='value(status.url)'
```

It is deployed **`--allow-unauthenticated`** so it works from a browser or a
bare `curl`. That also means anyone with the URL can start a screen. Close it
when the demo is done:

```bash
gcloud run services update cart-platform --region us-central1 \
    --no-allow-unauthenticated
```

---

## 5. Cancer type to final result

Nine calls. Substitute your URL.

```bash
BASE=https://cart-platform-<hash>-uc.a.run.app
```

### 0 — Ask what it can answer for

```bash
curl -s "$BASE/indications"
```

Returns every configured indication with the cohort, atlas, dependency lineage
and normal denominator behind it. Call this first rather than guessing a cancer
type: anything not on the list is refused at creation.

### 1 — Create a project

```bash
curl -s -X POST "$BASE/projects" \
  -H 'Content-Type: application/json' \
  -d '{"cancer_type":"Pancreatic Ductal Adenocarcinoma"}'
```

```json
{
  "project_id": "705e8ee1eced",
  "cancer_type": "Pancreatic Ductal Adenocarcinoma",
  "target_antigen": null,
  "discovery_mode": "B"
}
```

**`target_antigen` is null on purpose and stays null.** No target is ever seeded;
that null is what selects discovery mode. Keep the `project_id`.

Two indications are configured. Any other cancer type is refused at creation
with `400` naming what is available, rather than answered with these results
under a different name. Short aliases work — `PDAC`, `breast` — so a caller
need not reproduce the full oncological name.

### 2 — Start the screen

```bash
PID=705e8ee1eced
curl -s -X POST "$BASE/projects/$PID/runs"
```

Returns `202` **immediately** with a `job_id`. It does not wait — a screen scores
3,466 proteins and evaluates 19,900 pairs, taking about six minutes, so a
synchronous endpoint would time out in every client. It reads the derived
single-cell summaries, not the 8.3 GB matrix they came from.

### 3 — Poll

```bash
JID=2dbd755425cb
until curl -s "$BASE/jobs/$JID" | grep -qE '"status": "(complete|failed)"'; do
  sleep 5
done
curl -s "$BASE/jobs/$JID"
```

Match **both** terminal states. Waiting only for `complete` spins forever on a
failed job while discarding the body that explains it.

Progresses through `sources → screen → pairing → binders → constructs → safety
→ developability → ranking → package`, ending at
`"status": "complete", "note": "RANKED"`.

### 4 — Ranked targets

```bash
curl -s "$BASE/projects/$PID/targets?limit=10"
```

3,466 proteins screened, 3,400 scored. Top result is **CEACAM5**, composite
0.8769, with a six-component breakdown.

`null` components are genuinely unmeasured, never imputed as zero — CEACAM5 has
two of them, and `measured_weight: 0.55` says its composite rests on 55% of the
evidence. Render them as "not measured", never as `0`.

### 5 — Pairs

```bash
curl -s "$BASE/projects/$PID/pairs?limit=10"
```

19,900 pairs evaluated. Each carries `coverage_f_ab` **and**
`coverage_span_percentile`, with `coverage_caveat` on the number itself — that
fraction correlates with genomic span more than with expression, so the raw
value alone would mislead.

### 6 — Constructs · **the one to read**

```bash
curl -s "$BASE/projects/$PID/constructs"
```

```json
{
  "status": "BUILDABLE",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "counts": {
    "NO_CONSTRUCT": 195,
    "BUILDABLE": 5
  },
  "buildable": 5,
  "complete": 5,
  "awaiting_binder": 0,
  "over_budget": 0,
  "constructs": [
    {
      "gene": "FER1L6",
      "verdict": "BUILDABLE",
      "state": "COMPLETE",
      "design_class": "ADVANCED",
      "architecture": "adaptor, anti-tag receptor, antigen on the adaptor",
      "binder_supplied": true,
      "binder": "anti-tag binder, peptide neo-epitope, GCN4(7P-14P) (PDB 1P4B entities 1+2, antigen entity 3 excluded)",
      "total_bp": 2868,
      "budget_bp": 3500,
      "headroom_bp": 632,
      "amino_acid_sequence": "<955 residues, elided here>",
      "dna": "<2868 bases, elided here>"
    }
  ]
}
```

The two sequence fields are elided above for length only; the response carries
them in full.

> ### This is HTTP 200, and so is the empty case.
>
> Five designs fit the budget and every one carries a sequence. They route to
> an adaptor receptor, which binds a tag rather than the antigen.
> **The tag-binding sequence is a murine antibody fragment taken from a public
> structure and used exactly as deposited, crystallisation artifacts included.**
> It is not the molecule to order, nothing here has assessed its immunogenicity,
> and no developability figure in the platform describes it. Each design is also
> two manufactured biologics, not one.
>
> `AWAITING_BINDER` remains a third state the API can return, deliberately not
> merged with either neighbour: a design that fits the budget and has no
> residues is neither one that failed to fit nor a finished one. This pool
> returns none of them.
>
> **Run the same calls against breast and the answer is
> `NO_BUILDABLE_CONSTRUCT` with none buildable — still HTTP 200, still with
> reasons.** Never a 404, never a 500, never a bare `[]`. An empty list would
> read as "we looked and had nothing to say." Something specific and measured
> stopped each design, and naming it is the answer.

### 7 — The end state

```bash
curl -s "$BASE/projects/$PID/result"
```

```json
{
  "status": "RANKED",
  "usability": "USABLE",
  "unavailable": [],
  "indication": "Pancreatic Ductal Adenocarcinoma",
  "pool_size": 200,
  "reached_the_end": 5,
  "complete": 5,
  "awaiting_binder": 0,
  "attrition": [
    {
      "gate": "blocked on normal tissue risk",
      "dropped": 192,
      "remaining": 8
    },
    {
      "gate": "no design recommended",
      "dropped": 0,
      "remaining": 8
    },
    {
      "gate": "no binder retrieved",
      "dropped": 3,
      "remaining": 5
    },
    {
      "gate": "no construct assembled",
      "dropped": 0,
      "remaining": 5
    },
    {
      "gate": "construct over budget",
      "dropped": 0,
      "remaining": 5
    }
  ]
}
```

Every one of the 200 candidates is attributed to the **first** gate it failed,
so the drops sum to the pool rather than overlapping.

Breast returns the same shape with different numbers and a different end state:

```json
{
  "status": "NO_DESIGN_REACHES_THE_END",
  "indication": "Invasive Breast Carcinoma",
  "pool_size": 200,
  "reached_the_end": 0,
  "attrition": [
    {
      "gate": "blocked on normal tissue risk",
      "dropped": 196,
      "remaining": 4
    },
    {
      "gate": "no design recommended",
      "dropped": 0,
      "remaining": 4
    },
    {
      "gate": "no binder retrieved",
      "dropped": 3,
      "remaining": 1
    },
    {
      "gate": "no construct assembled",
      "dropped": 1,
      "remaining": 0
    },
    {
      "gate": "construct over budget",
      "dropped": 0,
      "remaining": 0
    }
  ]
}
```

**A deployment that only ever screens PDAC never sees this path.** The
no-survivor branch of the service — the refusal reasons, the attrition
explanation, the conservative-backup counts — is reached only by an indication
that returns nothing, and it is worth exercising once after deploying.

`usability`, `unavailable` and `indication` ride on every collection response,
not just this one — so a caller reading `/targets` never has to know to ask
elsewhere whether the ranking beneath it is supported. Full reference in
[API.md](API.md).

### 8 — Evidence for one gene

```bash
curl -s "$BASE/projects/$PID/evidence/MSLN"
```

Every stage's view of a single target — screen, pairing, binders, construct,
safety, developability, ranking — with the provenance behind each.

---

### 9 — The candidate package · **what a reader actually receives**

```bash
curl -s "$BASE/projects/$PID/package"
```

One package per surviving design, assembled from the stages that produced it
and computing nothing new. For PDAC that is `PACKAGED` with five packages; for
breast no candidate reaches the end, so the endpoint says so rather than
returning an empty list.

Each package carries eight sections — ranking, design class, construct with its
full sequence and domain map, target evidence with the per-organ risk
attribution, binders by route, safety, developability, and an experimental
validation plan — plus the release pin of every data source and the
configuration-hash chain that identifies the run.

**It also carries what it cannot tell you.** The reference specification asks
for twelve deliverables; eight have something to carry, and the other four are
named per deliverable rather than omitted — 26 elements across eight, of which
16 are checked mechanically on every run so the list cannot quietly go stale.

`GET /projects/$PID/package/FER1L6` returns one package whole.

---

## 6. What this deployment cannot do

**Only the indications that are configured.** Two are: pancreatic ductal
adenocarcinoma and invasive breast carcinoma. The deployment holds both without
rebuilding, because every tumour-side cache is namespaced by indication in both
its filename and its key. A third needs its cohort, atlas, dependency lineage
and normal denominator declared in `car_pipeline/configs/`, and its caches
provisioned by `bootstrap.py --from-sources`, which iterates the registry.

That provisioning is the slow part, not the code: deriving the single-cell
summaries streams the whole atlas, which for PDAC took **2 h 19 min**. The
container never carries those atlases — 8.3 GB for one, 844 MB for the other.
They are build-time inputs; the derived summaries are what a served run reads.

**Where a source does not exist for an indication, the pipeline degrades rather
than substituting.** It names the components it could not measure, and if the
missing evidence undermines the ranking it refuses to present one at all —
`usability: NOT_USABLE` with reasons, not numbers with a caveat attached. Every
collection endpoint carries that judgement, so a caller reading `/targets`
never has to know to ask somewhere else whether the ranking beneath it is
supported.

**Two of the twelve deliverables do not exist**, and the difference between
them matters for planning. The **structural report** (Stage 7) was never built,
but it is buildable from what is already connected — it needs only sequences the
platform emits. **Functional predictions** (Stage 8) are not buildable here at
all: the reference specification names partner-generated experimental data among
the required inputs, and no such data is connected. That is the one gap a
decision alone cannot close. Neither stage is stubbed; the numbering shows the
gap rather than renumbering around it, and each package names both.

**One run at a time**, globally, and jobs do not survive a restart. Both follow
from the same in-memory job table.

---

## 7. Running it locally instead

No deployment needed.

```bash
.venv\Scripts\python.exe -m car_pipeline.api.server     # Windows
.venv/bin/python -m car_pipeline.api.server             # macOS / Linux
```

`http://127.0.0.1:8000`, loopback only. Nothing needs to be running first — no
database, no broker, no worker. Same curl sequence against `BASE=http://127.0.0.1:8000`.

**Use the interpreter in `.venv`.** A bare `python` or `py -3.13` is the system
Python, which does not have this project's three dependencies and fails on the
first import.

To verify the whole pipeline instead of just serving it:

```bash
.venv\Scripts\python.exe run_all.py --fresh
```

Twelve stages, **118/124 criteria clear**, about 35 minutes with `--fresh` and
under ten without. The six that trip are recorded limitations, each with the
decision that accepted it; a criterion tripping *without* one is reported as a
regression and exits non-zero.

For a single page showing the whole run, `make_artifact.py` renders
`reports/full-run.html` from the run's own JSON.
