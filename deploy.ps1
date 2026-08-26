<#
    Deploy the design platform to Cloud Run.

    Everything before the sign-in is done: the container is defined and its exact
    file set has been run end to end. This script is the rest.

    Run it after:
        gcloud auth login
        gcloud config set project <PROJECT_ID>

    Those two need a browser and your own Google account, which is why they are
    not in here.

        .\deploy.ps1 -CheckOnly            # report readiness, deploy nothing
        .\deploy.ps1                       # deploy and smoke-test
        .\deploy.ps1 -Region europe-west1  # somewhere else
#>
param(
    [string]$Service = "cart-platform",
    [string]$Region  = "us-central1",
    [string]$Project = "",
    # Report readiness and stop. Deploys nothing, bills nothing.
    [switch]$CheckOnly
)

# Deliberately NOT "Stop". In Windows PowerShell 5.1 a native executable's
# stderr becomes an ErrorRecord when its output is captured, and gcloud writes
# routine notices there: "Your active configuration is: [default]", component
# update nudges, filter warnings. Under "Stop" the first of those is terminating
# and the script dies pointing at a harmless line. Exit codes are checked
# explicitly instead, which is what actually indicates failure for a native
# command.
$ErrorActionPreference = "Continue"

# Prefer whatever is on PATH; fall back to the per-user install location. Both
# are checked so this still works in a shell opened before the SDK was
# installed, which is the single most common way this looks broken.
$gcloud = (Get-Command gcloud -ErrorAction SilentlyContinue).Source
if (-not $gcloud) {
    $gcloud = "$env:LOCALAPPDATA\Google\Cloud SDK\google-cloud-sdk\bin\gcloud.cmd"
}
if (-not (Test-Path $gcloud)) {
    Write-Host "gcloud not found." -ForegroundColor Red
    Write-Host "Open a NEW terminal - an existing one keeps the PATH it started"
    Write-Host "with, from before the SDK was installed."
    exit 1
}

function Step($text) { Write-Host "`n=== $text ===" -ForegroundColor Cyan }

function Invoke-GCloud {
    # Run gcloud, return stdout, fail loudly on a nonzero exit.
    param([string[]]$GcArgs, [string]$What, [switch]$AllowFailure)
    $out = & $gcloud @GcArgs 2>$null
    if ($LASTEXITCODE -ne 0 -and -not $AllowFailure) {
        Write-Host "`n$What failed (gcloud exit $LASTEXITCODE)." -ForegroundColor Red
        Write-Host "Re-run this to see gcloud's own message:"
        Write-Host "    gcloud $($GcArgs -join ' ')"
        exit 1
    }
    return $out
}

Write-Host "gcloud: $gcloud"

Step "Checking sign-in"
# Not --filter=status:ACTIVE: with zero accounts gcloud warns that the filter
# key matched no resource, which is noise on the one path where it matters.
$accounts = Invoke-GCloud @("auth", "list", "--format=value(account)") "Reading accounts" -AllowFailure
if (-not $accounts) {
    Write-Host "Not signed in." -ForegroundColor Red
    Write-Host "Run these two, then re-run this script:"
    Write-Host "    gcloud auth login"
    Write-Host "    gcloud config set project <PROJECT_ID>"
    exit 1
}
Write-Host "signed in as $accounts"

if ($Project) { Invoke-GCloud @("config", "set", "project", $Project) "Setting project" | Out-Null }
$Project = Invoke-GCloud @("config", "get-value", "project") "Reading project" -AllowFailure
if (-not $Project -or $Project -eq "(unset)") {
    Write-Host "No project set." -ForegroundColor Red
    Write-Host "    gcloud projects list                       # what you have"
    Write-Host "    gcloud config set project <PROJECT_ID>"
    exit 1
}
Write-Host "project $Project"

Step "Checking billing"
# Cloud Run and Cloud Build both require a billing account on the project.
# Checked up front, because without it the failure surfaces partway through the
# first deploy as something less obvious.
$billing = Invoke-GCloud @("beta", "billing", "projects", "describe", $Project,
                           "--format=value(billingEnabled)") "Reading billing" -AllowFailure
if ($LASTEXITCODE -ne 0) {
    Write-Host "Could not read billing status; the Cloud Billing API may be off." -ForegroundColor Yellow
    Write-Host "  Check: https://console.cloud.google.com/billing/linkedaccount?project=$Project"
} elseif ("$billing" -ne "True") {
    Write-Host "Billing is NOT enabled on $Project. Cloud Run cannot deploy without it." -ForegroundColor Red
    Write-Host "  Link an account: https://console.cloud.google.com/billing/linkedaccount?project=$Project"
    exit 1
} else {
    Write-Host "billing enabled"
}

if ($CheckOnly) {
    Write-Host "`nReady to deploy. Re-run without -CheckOnly." -ForegroundColor Green
    exit 0
}

Step "Enabling the APIs the deploy needs"
# run: the service itself. cloudbuild + artifactregistry: --source builds the
# image remotely and pushes it, so both are required even though neither is
# named in the deploy command.
Invoke-GCloud @("services", "enable", "run.googleapis.com",
                "cloudbuild.googleapis.com", "artifactregistry.googleapis.com",
                "--project", $Project) "Enabling APIs" | Out-Null
Write-Host "run, cloudbuild, artifactregistry enabled"

Step "Deploying"
# --no-cpu-throttling  the screen runs on a background thread after the 202 has
#                      been sent; the default allocates CPU only while a request
#                      is in flight and would throttle that thread to near zero.
#                      The symptom is not an error - the job sits at "running"
#                      and the stage never advances.
# --concurrency=1      a second line of defence, not the fix. Runs are detached
#                      threads and the request returns in milliseconds, so a
#                      request-concurrency cap never sees two runs overlap. What
#                      actually serialises them is the global guard in
#                      start_run, which rejects a run while any other is in
#                      flight rather than only one for the same project. This
#                      flag bounds polling load and nothing more.
# --max-instances=1    the job table is a dict in memory, so a poll must reach
#                      the instance that is running the job.
# --min-instances=1    scale-to-zero loses the job table between run and poll,
#                      and avoids pulling a ~1.1 GB image on a cold start.
& $gcloud run deploy $Service `
    --source . `
    --region $Region `
    --project $Project `
    --min-instances=1 `
    --max-instances=1 `
    --concurrency=1 `
    --no-cpu-throttling `
    --cpu=4 `
    --memory=8Gi `
    --timeout=900 `
    --allow-unauthenticated
if ($LASTEXITCODE -ne 0) {
    Write-Host "`nDeploy failed; the build log above explains why." -ForegroundColor Red
    exit 1
}

$url = Invoke-GCloud @("run", "services", "describe", $Service, "--region", $Region,
                       "--project", $Project, "--format=value(status.url)") "Reading the service URL"
if (-not $url) {
    Write-Host "Deployed but no URL returned." -ForegroundColor Red
    exit 1
}

Step "Smoke test"
Write-Host "URL: $url"
# Not $project: PowerShell names are case-insensitive, so that would overwrite
# the $Project parameter holding the GCP project id.
$created = Invoke-RestMethod -Method Post -Uri "$url/projects" `
    -ContentType "application/json" `
    -Body '{"cancer_type":"Pancreatic Ductal Adenocarcinoma"}'
Write-Host "created project $($created.project_id), target_antigen=$($created.target_antigen)"

$job = Invoke-RestMethod -Method Post -Uri "$url/projects/$($created.project_id)/runs"
Write-Host "submitted job $($job.job_id)"

$deadline = (Get-Date).AddMinutes(20)
do {
    Start-Sleep -Seconds 5
    $state = Invoke-RestMethod -Uri "$url/jobs/$($job.job_id)"
    Write-Host "  $($state.status)  $($state.stage)"
} while ($state.status -notin @("complete", "failed") -and (Get-Date) -lt $deadline)

if ($state.status -ne "complete") {
    # $state.error is null while a job is merely slow, so printing only that
    # gives "Job did not complete: " and no diagnostic. The stage and note say
    # where it actually stopped.
    if ($state.status -eq "failed") {
        Write-Host "Job failed at stage '$($state.stage)': $($state.error)" -ForegroundColor Red
    } else {
        Write-Host "Job still '$($state.status)' at stage '$($state.stage)' ($($state.note)) after 20 minutes." -ForegroundColor Red
    }
    exit 2
}

$result = Invoke-RestMethod -Uri "$url/projects/$($created.project_id)/result"
Write-Host "`nresult: $($result.status)"
foreach ($gate in $result.attrition) {
    Write-Host ("  {0,-34} -{1,4}  {2,4} remain" -f $gate.gate, $gate.dropped, $gate.remaining)
}

Write-Host "`nDeployed: $url" -ForegroundColor Green
