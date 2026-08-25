<#
    Deploy the design platform to Cloud Run.

    Everything before the sign-in is already done: the container is defined and
    its exact file set has been run end to end. This script is the rest.

    Run it after:
        gcloud auth login
        gcloud config set project <PROJECT_ID>

    Those two need a browser and your own Google account, which is why they are
    not in here.

        .\deploy.ps1                       # deploy and smoke-test
        .\deploy.ps1 -Region europe-west1  # somewhere else
#>
param(
    [string]$Service = "cart-platform",
    [string]$Region  = "us-central1",
    [string]$Project = ""
)

$ErrorActionPreference = "Stop"

$gcloud = "$env:LOCALAPPDATA\Google\Cloud SDK\google-cloud-sdk\bin\gcloud.cmd"
if (-not (Test-Path $gcloud)) { $gcloud = "gcloud" }

function Step($text) { Write-Host "`n=== $text ===" -ForegroundColor Cyan }

Step "Checking sign-in"
$accounts = & $gcloud auth list --filter=status:ACTIVE --format="value(account)"
if (-not $accounts) {
    Write-Host "Not signed in. Run these two first, then re-run this script:"
    Write-Host "    $gcloud auth login"
    Write-Host "    $gcloud config set project <PROJECT_ID>"
    exit 1
}
Write-Host "signed in as $accounts"

if ($Project) { & $gcloud config set project $Project | Out-Null }
# No 2>$null here: redirecting a native command's stderr in PowerShell 5.1
# wraps each line in a NativeCommandError, and gcloud writes "Your active
# configuration is: [default]" to stderr on essentially every call. Under
# $ErrorActionPreference = "Stop" that is terminating, and the script would
# die here pointing at a harmless status line.
$Project = & $gcloud config get-value project
if (-not $Project -or $Project -eq "(unset)") {
    Write-Host "No project set. Run: $gcloud config set project <PROJECT_ID>"
    exit 1
}
Write-Host "project $Project"

Step "Enabling the APIs the deploy needs"
# run: the service itself. cloudbuild + artifactregistry: --source builds the
# image remotely and pushes it, so both are required even though neither is
# named in the deploy command.
& $gcloud services enable run.googleapis.com cloudbuild.googleapis.com `
    artifactregistry.googleapis.com --project $Project
# $ErrorActionPreference does not apply to a native executable's exit code in
# Windows PowerShell 5.1, so this has to be explicit. Without it a failed deploy
# runs the smoke test anyway and reports a URI parse error instead of the build
# log that explains what actually went wrong.
if ($LASTEXITCODE -ne 0) { Write-Host "Enabling APIs failed." -ForegroundColor Red; exit 1 }

Step "Deploying"
# --no-cpu-throttling  the screen runs on a background thread after the 202 has
#                      been sent; the default allocates CPU only while a request
#                      is in flight and would throttle that thread to near zero.
#                      The symptom is not an error — the job sits at "running"
#                      and the stage never advances.
# --concurrency=1      a second line of defence, not the fix. Runs are detached
#                      threads and the request returns in milliseconds, so a
#                      request-concurrency cap never sees two runs overlap. What
#                      actually serialises them is the global guard in
#                      start_run, which now rejects a run while any other is in
#                      flight rather than only one for the same project. This
#                      flag bounds the polling load and nothing more.
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
if ($LASTEXITCODE -ne 0) { Write-Host "Deploy failed; the build log above explains why." -ForegroundColor Red; exit 1 }

$url = & $gcloud run services describe $Service --region $Region `
    --project $Project --format="value(status.url)"
if ($LASTEXITCODE -ne 0 -or -not $url) { Write-Host "Deployed but no URL returned." -ForegroundColor Red; exit 1 }

Step "Smoke test"
Write-Host "URL: $url"
# Not $project: PowerShell names are case-insensitive, so that would
# overwrite the $Project parameter holding the GCP project id.
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
