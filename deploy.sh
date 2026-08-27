#!/usr/bin/env bash
# Deploy the design platform to Cloud Run. The bash equivalent of deploy.ps1,
# for macOS, Linux, WSL, Git Bash and Cloud Shell.
#
# Run it after:
#     gcloud auth login
#     gcloud config set project <PROJECT_ID>
#
# Those two need a browser and your own Google account, which is why they are
# not in here.
#
#     ./deploy.sh --check-only        # report readiness, deploy nothing
#     ./deploy.sh                     # deploy and smoke-test
#     REGION=europe-west1 ./deploy.sh # somewhere else
#
# Needs only bash, curl and gcloud. Deliberately no Python: this script must
# work on a machine that is deploying but never runs the pipeline locally.

set -euo pipefail

SERVICE="${SERVICE:-cart-platform}"
REGION="${REGION:-us-central1}"
CHECK_ONLY=0

# Every argument is checked. Falling through on an unrecognised one would
# provision an always-on 4 vCPU / 8 GiB instance because of a typo - including
# `-CheckOnly`, the PowerShell spelling, which sits one line away in DEPLOY.md.
while [ $# -gt 0 ]; do
    case "$1" in
        --check-only) CHECK_ONLY=1 ;;
        -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
        *)
            echo "Unrecognised argument: $1"
            echo "Usage: ./deploy.sh [--check-only]"
            echo "  (--check-only, not -CheckOnly; that is the PowerShell spelling)"
            exit 1
            ;;
    esac
    shift
done

step() { printf '\n=== %s ===\n' "$1"; }

# Pull one string field out of a JSON body. Enough for the four flat fields this
# script reads, and it keeps the dependency list at curl.
jsonval() { grep -o "\"$1\"[[:space:]]*:[[:space:]]*\"[^\"]*\"" | head -1 | sed -E 's/.*:[[:space:]]*"([^"]*)".*/\1/'; }

command -v gcloud >/dev/null 2>&1 || {
    echo "gcloud not found. Install the Google Cloud CLI:"
    echo "    https://cloud.google.com/sdk/docs/install"
    echo "If you just installed it, open a NEW shell - an existing one keeps"
    echo "the PATH it started with."
    exit 1
}
command -v curl >/dev/null 2>&1 || { echo "curl not found."; exit 1; }

step "Checking sign-in"
# Not --filter=status:ACTIVE: with zero accounts gcloud warns that the filter
# key matched no resource, which is noise on the one path where it matters.
ACCOUNT="$(gcloud auth list --format='value(account)' 2>/dev/null | head -1 || true)"
if [ -z "$ACCOUNT" ]; then
    echo "Not signed in. Run these two, then re-run this script:"
    echo "    gcloud auth login"
    echo "    gcloud config set project <PROJECT_ID>"
    exit 1
fi
echo "signed in as $ACCOUNT"

PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
    echo "No project set."
    echo "    gcloud projects list                  # what you have"
    echo "    gcloud config set project <PROJECT_ID>"
    exit 1
fi
echo "project $PROJECT"

step "Checking billing"
# Cloud Run and Cloud Build both require a billing account on the project.
# Checked up front: without it the failure surfaces partway through the first
# deploy as something less obvious, and the fix is a console page, not a flag.
if [ "${SKIP_BILLING_CHECK:-0}" = "1" ]; then
    echo "skipped (SKIP_BILLING_CHECK=1)"
elif ! BILLING="$(gcloud beta billing projects describe "$PROJECT" \
        --format='value(billingEnabled)' 2>/dev/null)"; then
    # Unknown is not enabled, so this stops rather than waving the deploy
    # through. But it is genuinely ambiguous: reading billing needs a permission
    # on the *billing account* that a project Owner often does not hold, and it
    # also fails when the beta component is absent. Hence the override, which
    # the enabled=false branch below deliberately does not offer.
    echo "Could not read billing status for $PROJECT."
    echo "  This is ambiguous: billing may be fine and you may simply lack"
    echo "  permission to read it, or the 'beta' component may be missing."
    echo "  Check by hand:  gcloud beta billing projects describe $PROJECT"
    echo "  Or in console:  https://console.cloud.google.com/billing/linkedaccount?project=$PROJECT"
    echo "  If you know billing is on:  SKIP_BILLING_CHECK=1 ./deploy.sh"
    exit 1
elif [ "$BILLING" != "True" ]; then
    echo "Billing is NOT enabled on $PROJECT. Cloud Run cannot deploy without it."
    echo "  Link an account: https://console.cloud.google.com/billing/linkedaccount?project=$PROJECT"
    exit 1
else
    echo "billing enabled"
fi

if [ "$CHECK_ONLY" = "1" ]; then
    printf '\nReady to deploy. Re-run without --check-only.\n'
    exit 0
fi

step "Enabling the APIs the deploy needs"
# run: the service itself. cloudbuild + artifactregistry: --source builds the
# image remotely and pushes it, so both are required even though neither is
# named in the deploy command.
gcloud services enable run.googleapis.com cloudbuild.googleapis.com \
    artifactregistry.googleapis.com --project "$PROJECT"

step "Deploying"
# --no-cpu-throttling  the screen runs on a background thread after the 202 has
#                      been sent; the default allocates CPU only while a request
#                      is in flight and would throttle that thread to near zero.
#                      The symptom is not an error - the job sits at "running"
#                      and the stage never advances.
# --concurrency=1      a second line of defence, not the fix. Runs are detached
#                      threads and the request returns in milliseconds, so a
#                      request-concurrency cap never sees two runs overlap. What
#                      serialises them is the global guard in start_run.
# --max-instances=1    the job table is a dict in memory, so a poll must reach
#                      the instance running the job.
# --min-instances=1    scale-to-zero loses the job table between run and poll.
#
# No inline comments inside the command: in bash a backslash followed by a
# space escapes the space, not the newline, so a trailing "# ..." would end the
# command there and silently drop every flag after it.
gcloud run deploy "$SERVICE" \
    --source . \
    --region "$REGION" \
    --project "$PROJECT" \
    --min-instances=1 \
    --max-instances=1 \
    --concurrency=1 \
    --no-cpu-throttling \
    --cpu=4 \
    --memory=8Gi \
    --timeout=900 \
    --allow-unauthenticated

URL="$(gcloud run services describe "$SERVICE" --region "$REGION" \
    --project "$PROJECT" --format='value(status.url)')"
[ -n "$URL" ] || { echo "Deployed but no URL returned."; exit 1; }

step "Smoke test"
echo "URL: $URL"
PID="$(curl -s -X POST "$URL/projects" -H 'Content-Type: application/json' \
    -d '{"cancer_type":"Pancreatic Ductal Adenocarcinoma"}' | jsonval project_id)"
[ -n "$PID" ] || { echo "Could not create a project."; exit 2; }
echo "created project $PID"

JID="$(curl -s -X POST "$URL/projects/$PID/runs" | jsonval job_id)"
[ -n "$JID" ] || { echo "Could not start a run."; exit 2; }
echo "submitted job $JID"

# Anything outside the known progress states is terminal. Waiting only for
# complete|failed spins the full deadline on NOT_FOUND - which is what a
# restarted instance returns, since the job table is in memory - while the
# body explaining it is discarded every five seconds.
DEADLINE=$(( $(date +%s) + 1200 ))
while true; do
    STATE="$(curl -s "$URL/jobs/$JID")"
    STATUS="$(printf '%s' "$STATE" | jsonval status)"
    STAGE="$(printf '%s' "$STATE" | jsonval stage)"
    echo "  ${STATUS:-<no status>}  $STAGE"
    case "$STATUS" in
        queued|running) ;;
        *) break ;;
    esac
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then
        echo "Still '$STATUS' at stage '$STAGE' after 20 minutes."
        exit 2
    fi
    sleep 5
done

if [ "$STATUS" != "complete" ]; then
    echo "Job did not complete. Status '$STATUS':"
    printf '%s\n' "$STATE"
    exit 2
fi

printf '\n'
curl -s "$URL/projects/$PID/result"
printf '\n\nDeployed: %s\n' "$URL"
