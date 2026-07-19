#!/usr/bin/env bash
# =============================================================================
# common.sh — Shared helpers for the GLIMPSE2 metagenomic imputation pipeline
#
# Sourced by every stage script. Provides:
#   - Checkpointing (resume-safe reruns per sample or per chromosome)
#   - Consistent logging
#   - Small validation helpers
#
# This file defines functions only — it is not meant to be executed directly.
# =============================================================================

# ---- Logging ----------------------------------------------------------------
# Usage: log "some message"
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
    exit 1
}

# ---- Checkpointing ------------------------------------------------------------
# Each stage script sets CHECKPOINT_DIR and CHECKPOINT_KEY (e.g. a sample ID or
# "chr${CHR}") before sourcing this file / calling these functions.
#
# checkpoint_done  <stage_name>   marks a stage complete
# checkpoint_exists <stage_name>  returns 0 (true) if already done
#
# IMPORTANT — lesson learned in production use:
# checkpoint_exists only checks that a marker file is present. It does NOT
# verify that the underlying output is still valid for the CURRENT input set.
# If you change the sample list, reference panel, or any upstream input
# between runs, you must explicitly clear the relevant checkpoints — see
# docs/TROUBLESHOOTING.md ("Sample set changed between runs") for a real
# incident this caused (mismatched sample counts across GLIMPSE2 chunks).

checkpoint_done() {
    local stage="$1"
    mkdir -p "${CHECKPOINT_DIR}"
    date '+%Y-%m-%d %H:%M:%S' > "${CHECKPOINT_DIR}/${CHECKPOINT_KEY}.${stage}.done"
    log "  Checkpoint written: ${stage} [${CHECKPOINT_KEY}]"
}

checkpoint_exists() {
    local stage="$1"
    local marker="${CHECKPOINT_DIR}/${CHECKPOINT_KEY}.${stage}.done"
    if [ -f "${marker}" ]; then
        log "  Already done: ${stage} [${CHECKPOINT_KEY}] ($(cat "${marker}")) — skipping"
        return 0
    fi
    return 1
}

clear_checkpoint() {
    local stage="$1"
    rm -f "${CHECKPOINT_DIR}/${CHECKPOINT_KEY}.${stage}.done"
}

# ---- Validation helpers -------------------------------------------------------

require_file() {
    local f="$1"
    local msg="${2:-Required file missing}"
    [ -f "${f}" ] || die "${msg}: ${f}"
}

require_nonempty_file() {
    local f="$1"
    local msg="${2:-Required file missing or empty}"
    [ -s "${f}" ] || die "${msg}: ${f}"
}

# require_paired_fastq_counts <R1> <R2>
# Guards against the class of bug where R1/R2 read counts silently diverge
# (e.g. one file truncated by an interrupted job) and Bowtie2 dies deep into
# an alignment run with "fewer reads in file specified with -1/-2".
require_paired_fastq_counts() {
    local r1="$1" r2="$2"
    local n1 n2
    n1=$(zcat "${r1}" 2>/dev/null | awk 'END{print NR/4}')
    n2=$(zcat "${r2}" 2>/dev/null | awk 'END{print NR/4}')
    if [ "${n1}" != "${n2}" ]; then
        die "Paired FASTQ read count mismatch: $(basename "${r1}")=${n1} reads vs $(basename "${r2}")=${n2} reads. See docs/TROUBLESHOOTING.md."
    fi
}
