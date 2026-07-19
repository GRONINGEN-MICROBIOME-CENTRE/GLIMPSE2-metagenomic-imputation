#!/usr/bin/env bash
# =============================================================================
# run_full_pipeline.sh — Submit the full pipeline with correct SLURM
# dependencies: preprocessing -> imputation (array) -> merge/QC.
#
# USAGE:
#   bash submit/run_full_pipeline.sh path/to/sample_list.txt
#
# This submits everything at once using --dependency so each stage only
# starts once the previous one has fully succeeded. For iterative
# development / debugging, it's usually clearer to run each stage's script
# manually and inspect results before moving on — see docs/USAGE.md.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi
SAMPLE_LIST="$1"

# ---- Stage 1: preprocessing, one job per sample ----
PREPROCESS_JOB_IDS=()
while read -r SAMPLE; do
    [ -z "${SAMPLE}" ] && continue
    JOBID=$(sbatch --parsable "${REPO_ROOT}/scripts/01_preprocess_sample.sh" "${SAMPLE}")
    PREPROCESS_JOB_IDS+=("${JOBID}")
    echo "Submitted preprocessing for ${SAMPLE}: job ${JOBID}"
done < "${SAMPLE_LIST}"

DEPENDENCY=$(IFS=:; echo "afterok:${PREPROCESS_JOB_IDS[*]}")

# ---- Stage 2: imputation, one array task per chromosome, after all preprocessing ----
IMPUTE_JOBID=$(sbatch --parsable --dependency="${DEPENDENCY}" \
    --array=1-22 "${REPO_ROOT}/scripts/02_impute_chromosome.sh")
echo "Submitted imputation array: job ${IMPUTE_JOBID}"

# ---- Stage 3: merge + QC, after all chromosomes ----
MERGE_JOBID=$(sbatch --parsable --dependency="afterok:${IMPUTE_JOBID}" \
    "${REPO_ROOT}/scripts/03_merge_and_qc.sh")
echo "Submitted merge/QC: job ${MERGE_JOBID}"

echo ""
echo "Pipeline submitted. Monitor with: squeue -u \$USER"
echo "Run scripts/04_validate_against_array.sh manually afterward if you have array data."
