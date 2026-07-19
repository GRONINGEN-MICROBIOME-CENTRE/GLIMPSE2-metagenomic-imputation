#!/usr/bin/env bash
# =============================================================================
# submit_all_samples.sh — Submit scripts/01_preprocess_sample.sh for every
# sample in a list file (one sample ID per line).
#
# USAGE:
#   bash submit/submit_all_samples.sh path/to/sample_list.txt
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi
SAMPLE_LIST="$1"

while read -r SAMPLE; do
    [ -z "${SAMPLE}" ] && continue
    sbatch "${REPO_ROOT}/scripts/01_preprocess_sample.sh" "${SAMPLE}"
done < "${SAMPLE_LIST}"
