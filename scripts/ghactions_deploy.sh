#!/usr/bin/env bash

# Do NOT set -v or -x or your GitHub API token will be leaked!
#set -ue # exit with nonzero exit code if anything fails

echo "Parse memote.ini for values."
deployment=$(awk -F '=' '{if (! ($0 ~ /^;/) && $0 ~ /deployment/) print $2}' memote.ini | tr -d ' ')

git config --global user.email "github-actions[bot]"
git config --global user.name "41898282+github-actions[bot]@users.noreply.github.com"

# Always need the deployment branch available locally for storing results.
echo "Checking if deployment branch ${deployment} is available"
git checkout "${deployment}"
echo "Back to base ${HEAD_REF}"
git checkout "${HEAD_REF}"
echo "Tracked build."
mkdir -p results
memote run
echo "Start deploy to ${deployment}..."

# Generate the history report on the deployment branch.
output="index.html"
git checkout "${deployment}"
echo "Generating updated history report '${output}'."
memote report history --filename="${output}"

# Add, commit and push the files.
git add "${output}"
git commit -m "Github actions report # ${GITHUB_SHA}"
git push --quiet "https://github.com/${GITHUB_REPOSITORY}.git" "${deployment}" > /dev/null

echo "Your new report will be visible at https://wheeldon-lab.github.io/iyali26 in a moment."
