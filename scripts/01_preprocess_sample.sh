#!/usr/bin/env bash
#SBATCH --job-name=preprocess_sample
#SBATCH --partition=regularmedium
#SBATCH --cpus-per-task=24
#SBATCH --mem=99G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/slurm/preprocess_%A_%a.out

# =============================================================================
# 01_preprocess_sample.sh — Per-sample preprocessing
#
#   Trimmomatic (adapter trimming)
#   -> Bowtie2 (human alignment)
#   -> Picard MarkDuplicates
#   -> Bowtie2 (bacterial decontamination)
#
# Produces: <OUTDIR>/03_decontaminated/<SAMPLE>.decontaminated.bam
#
# USAGE:
#   sbatch scripts/01_preprocess_sample.sh <SAMPLE_ID>
#
# Expects paired FASTQs named <SAMPLE_ID>_R1.fastq.gz / _R2.fastq.gz under
# FASTQ_DIR (either directly, or in a per-sample subdirectory FASTQ_DIR/<SAMPLE_ID>/).
# If your raw files use a different convention (e.g. _1/_2), rename or symlink
# them first — see docs/TROUBLESHOOTING.md ("FASTQ naming mismatch").
#
# Resumable: re-running for the same sample skips stages already completed,
# tracked in CHECKPOINT_DIR.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

source "${REPO_ROOT}/config/config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

# =============================================================================
# ARGUMENT
# =============================================================================
if [ $# -lt 1 ]; then
    die "Usage: $0 <SAMPLE_ID>"
fi
SAMPLE="$1"
CHECKPOINT_KEY="${SAMPLE}"

set +u
source "${CONDA_SH}"
conda activate "${CONDA_ENV_ALIGN}"
set -u
export BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS_ALIGN}"

log "========================================"
log " 01_preprocess_sample.sh — sample: ${SAMPLE}"
log " Host: $(hostname)"
log "========================================"

# =============================================================================
# PATHS
# =============================================================================
TRIM_DIR="${OUTDIR}/01_trimmed"
ALIGN_DIR="${OUTDIR}/02_aligned"
DECON_DIR="${OUTDIR}/03_decontaminated"

mkdir -p "${CHECKPOINT_DIR}" "${TRIM_DIR}" "${ALIGN_DIR}" "${DECON_DIR}" \
         "${LOG_DIR}/slurm" "${LOG_DIR}/preprocess"

# Assigned up front (not just inside the "run fresh" branches) so that
# post-stage cleanup blocks can always reference them safely, even when a
# stage was already checkpointed and its main branch is skipped entirely.
BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
DECON_BAM="${DECON_DIR}/${SAMPLE}.decontaminated.bam"

# =============================================================================
# STAGE 1: NORMALIZE / STAGE INPUT FASTQs
# =============================================================================
if ! checkpoint_exists "stage1_normalize"; then
    log "STAGE 1: Format normalization"

    OUT_R1="${TRIM_DIR}/${SAMPLE}_R1.fastq.gz"
    OUT_R2="${TRIM_DIR}/${SAMPLE}_R2.fastq.gz"

    if [ -s "${OUT_R1}" ] && [ -s "${OUT_R2}" ]; then
        log "[${SAMPLE}] Already normalized — skipping"
    else
        SEARCH_DIR="${FASTQ_DIR}"
        [ -d "${FASTQ_DIR}/${SAMPLE}" ] && SEARCH_DIR="${FASTQ_DIR}/${SAMPLE}"
        log "[${SAMPLE}] Searching in: ${SEARCH_DIR}"

        R1_SRC=$(find -L "${SEARCH_DIR}" -name "${SAMPLE}_R1.fastq.gz" 2>/dev/null | head -1)
        R2_SRC=$(find -L "${SEARCH_DIR}" -name "${SAMPLE}_R2.fastq.gz" 2>/dev/null | head -1)

        if [ -z "${R1_SRC}" ] || [ -z "${R2_SRC}" ]; then
            die "R1/R2 not found for ${SAMPLE} in ${SEARCH_DIR} (expected *_R1.fastq.gz / *_R2.fastq.gz — see docs/TROUBLESHOOTING.md if your files use a different naming convention)"
        fi

        cp "${R1_SRC}" "${OUT_R1}"
        cp "${R2_SRC}" "${OUT_R2}"
        require_nonempty_file "${OUT_R1}" "R1 copy failed"
        require_nonempty_file "${OUT_R2}" "R2 copy failed"

        log "[${SAMPLE}] R1: $(du -sh "${OUT_R1}" | cut -f1)  R2: $(du -sh "${OUT_R2}" | cut -f1)"
    fi

    checkpoint_done "stage1_normalize"
fi

# =============================================================================
# STAGE 2: ADAPTER TRIMMING (Trimmomatic)
# =============================================================================
if ! checkpoint_exists "stage2_trimmomatic"; then
    log "STAGE 2: Adapter trimming (Trimmomatic)"

    R1_RAW="${TRIM_DIR}/${SAMPLE}_R1.fastq.gz"
    R2_RAW="${TRIM_DIR}/${SAMPLE}_R2.fastq.gz"
    TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1_paired.fastq.gz"
    TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2_paired.fastq.gz"
    TRIM_R1_U="${TRIM_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
    TRIM_R2_U="${TRIM_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"
    LOG="${LOG_DIR}/preprocess/${SAMPLE}_trimming.log"

    if [ -s "${TRIM_R1}" ] && [ -s "${TRIM_R2}" ]; then
        log "[${SAMPLE}] Already trimmed — skipping"
    else
        rm -f "${TRIM_R1}" "${TRIM_R2}" "${TRIM_R1_U}" "${TRIM_R2_U}"

        trimmomatic PE \
            -threads "${PREPROCESS_CPUS}" -phred33 \
            "${R1_RAW}" "${R2_RAW}" \
            "${TRIM_R1}" "${TRIM_R1_U}" \
            "${TRIM_R2}" "${TRIM_R2_U}" \
            ILLUMINACLIP:"${ADAPTERS_FASTA}":2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            >> "${LOG}" 2>&1

        require_nonempty_file "${TRIM_R1}" "Trimmomatic failed for ${SAMPLE}"

        # Guard against the class of bug where a killed/interrupted prior run
        # leaves a truncated (but non-empty) trimmed FASTQ that silently
        # passes downstream and breaks Bowtie2 hours later. See
        # docs/TROUBLESHOOTING.md ("Truncated trimmed FASTQ").
        require_paired_fastq_counts "${TRIM_R1}" "${TRIM_R2}"

        rm -f "${TRIM_R1_U}" "${TRIM_R2_U}" "${R1_RAW}" "${R2_RAW}"
    fi

    checkpoint_done "stage2_trimmomatic"
fi

TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1_paired.fastq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2_paired.fastq.gz"

# =============================================================================
# STAGE 3: ALIGNMENT (Bowtie2 -> human reference)
# Filters applied at alignment time: mapped (-F 4), properly paired (-f 2),
# MAPQ >= MIN_MAPQ.
# =============================================================================
if ! checkpoint_exists "stage3_alignment"; then
    log "STAGE 3: Alignment (Bowtie2 -> human)"

    LOG="${LOG_DIR}/preprocess/${SAMPLE}_alignment.log"

    if [ -s "${BAM}" ] && [ -f "${BAM}.bai" ]; then
        log "[${SAMPLE}] BAM already exists — skipping"
    else
        rm -f "${BAM}" "${BAM}.bai"

        [ -f "${HUMAN_BT2_IDX}.1.bt2" ] || [ -f "${HUMAN_BT2_IDX}.1.bt2l" ] \
            || die "Bowtie2 human index missing: ${HUMAN_BT2_IDX}"

        bowtie2 \
            -x "${HUMAN_BT2_IDX}" -1 "${TRIM_R1}" -2 "${TRIM_R2}" \
            -p "${PREPROCESS_CPUS}" --no-unal --very-sensitive \
            --rg-id "${SAMPLE}" --rg "SM:${SAMPLE}" --rg "PL:ILLUMINA" --rg "LB:lib1" \
            2>>"${LOG}" | \
        samtools view -bS -F 4 -f 2 -q "${MIN_MAPQ}" 2>>"${LOG}" | \
        samtools sort -@ 8 -o "${BAM}" >> "${LOG}" 2>&1

        require_nonempty_file "${BAM}" "Bowtie2 human alignment failed for ${SAMPLE}"
        samtools index "${BAM}" >> "${LOG}" 2>&1

        DEPTH=$(samtools coverage "${BAM}" 2>/dev/null | awk 'NR>1{s+=$7;n++} END{printf "%.3f", s/n}')
        log "[${SAMPLE}] Mean depth: ${DEPTH}x"
    fi

    checkpoint_done "stage3_alignment"
fi

if checkpoint_exists "stage3_alignment"; then
    rm -f "${TRIM_R1}" "${TRIM_R2}"
fi

# =============================================================================
# STAGE 4: MARKDUPLICATES + BACTERIAL DECONTAMINATION
# =============================================================================
if ! checkpoint_exists "stage4_decontamination"; then
    log "STAGE 4: MarkDuplicates + Decontamination"

    MARKDUP_BAM="${DECON_DIR}/${SAMPLE}.markdup.bam"
    METRICS="${DECON_DIR}/${SAMPLE}_markdup_metrics.txt"
    HUMAN_R1="${DECON_DIR}/${SAMPLE}_human_R1.fastq.gz"
    HUMAN_R2="${DECON_DIR}/${SAMPLE}_human_R2.fastq.gz"
    BACT_BAM="${DECON_DIR}/${SAMPLE}_bacterial_hits.bam"
    BACT_READS="${DECON_DIR}/${SAMPLE}_bacterial_read_names.txt"
    LOG="${LOG_DIR}/preprocess/${SAMPLE}_decontamination.log"

    if [ -s "${DECON_BAM}" ] && [ -f "${DECON_BAM}.bai" ]; then
        log "[${SAMPLE}] Already decontaminated — skipping"
    else
        rm -f "${MARKDUP_BAM}" "${MARKDUP_BAM}.bai" "${DECON_BAM}" "${DECON_BAM}.bai" \
              "${HUMAN_R1}" "${HUMAN_R2}" "${BACT_BAM}" "${BACT_BAM}.bai" "${BACT_READS}"

        picard MarkDuplicates \
            INPUT="${BAM}" OUTPUT="${MARKDUP_BAM}" METRICS_FILE="${METRICS}" \
            REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT \
            >> "${LOG}" 2>&1
        require_nonempty_file "${MARKDUP_BAM}" "MarkDuplicates failed for ${SAMPLE}"
        samtools index "${MARKDUP_BAM}" >> "${LOG}" 2>&1

        samtools fastq -@ 8 -1 "${HUMAN_R1}" -2 "${HUMAN_R2}" -0 /dev/null -s /dev/null -n \
            "${MARKDUP_BAM}" >> "${LOG}" 2>&1
        require_nonempty_file "${HUMAN_R1}" "FASTQ extraction failed for ${SAMPLE}"

        bowtie2 \
            -x "${BACTERIAL_BT2_IDX}" -1 "${HUMAN_R1}" -2 "${HUMAN_R2}" \
            -p "${PREPROCESS_CPUS}" --no-unal --very-sensitive 2>>"${LOG}" | \
        samtools view -bS -f 3 2>>"${LOG}" | \
        samtools sort -@ 8 -o "${BACT_BAM}" >> "${LOG}" 2>&1
        samtools index "${BACT_BAM}" >> "${LOG}" 2>&1

        samtools view "${BACT_BAM}" | awk '{print $1}' | sort -u > "${BACT_READS}" 2>>"${LOG}"
        N_BACT=$(wc -l < "${BACT_READS}")

        if [ "${N_BACT}" -eq 0 ]; then
            cp "${MARKDUP_BAM}" "${DECON_BAM}"
            cp "${MARKDUP_BAM}.bai" "${DECON_BAM}.bai"
        else
            picard FilterSamReads \
                INPUT="${MARKDUP_BAM}" OUTPUT="${DECON_BAM}" \
                READ_LIST_FILE="${BACT_READS}" FILTER=excludeReadList \
                VALIDATION_STRINGENCY=LENIENT >> "${LOG}" 2>&1
            require_nonempty_file "${DECON_BAM}" "FilterSamReads failed for ${SAMPLE}"
            samtools index "${DECON_BAM}" >> "${LOG}" 2>&1
        fi

        DEPTH_POST=$(samtools coverage "${DECON_BAM}" 2>/dev/null | awk 'NR>1{s+=$7;n++} END{printf "%.4f", s/n}')
        BREADTH=$(samtools coverage "${DECON_BAM}" 2>/dev/null | awk 'NR>1{s+=$6;n++} END{printf "%.2f", s/n*100}')
        MAPPED_POST=$(samtools flagstat "${DECON_BAM}" | grep "mapped (" | head -1 | awk '{print $1}')

        rm -f "${HUMAN_R1}" "${HUMAN_R2}" "${BACT_BAM}" "${BACT_BAM}.bai" \
              "${MARKDUP_BAM}" "${MARKDUP_BAM}.bai" "${BACT_READS}"

        log "[${SAMPLE}] Decontamination complete — depth=${DEPTH_POST}x breadth=${BREADTH}% mapped=${MAPPED_POST}"

        # Flag (do not silently drop) low-depth samples so a human decides
        # whether to include them — see docs/USAGE.md "Depth QC".
        BELOW_THRESHOLD=$(awk -v d="${DEPTH_POST}" -v t="${MIN_DEPTH_X}" 'BEGIN{print (d < t) ? "1" : "0"}')
        if [ "${BELOW_THRESHOLD}" = "1" ]; then
            log "[${SAMPLE}] WARNING: depth ${DEPTH_POST}x is below MIN_DEPTH_X=${MIN_DEPTH_X}x — flagged in depth_flags.tsv"
            echo -e "${SAMPLE}\t${DEPTH_POST}\tbelow_threshold" >> "${OUTDIR}/logs/depth_flags.tsv"
        fi

        echo -e "${SAMPLE}\t${DEPTH_POST}\t${BREADTH}\t${MAPPED_POST}" >> "${OUTDIR}/logs/decon_summary.tsv"
    fi

    checkpoint_done "stage4_decontamination"
fi

if checkpoint_exists "stage4_decontamination"; then
    if [ -s "${DECON_BAM}" ] && [ -f "${DECON_BAM}.bai" ]; then
        rm -f "${BAM}" "${BAM}.bai"
    fi
fi

log "========================================"
log " 01_preprocess_sample.sh complete — ${SAMPLE}"
log "========================================"
