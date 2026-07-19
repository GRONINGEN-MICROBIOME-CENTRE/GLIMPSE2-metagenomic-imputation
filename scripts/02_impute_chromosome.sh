#!/usr/bin/env bash
#SBATCH --job-name=impute_chr
#SBATCH --partition=regularlong
#SBATCH --cpus-per-task=60
#SBATCH --mem=99G
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/slurm/impute_%A_%a.out

# =============================================================================
# 02_impute_chromosome.sh — Joint GLIMPSE2 imputation, one chromosome per job
#
#   Stage 0: Panel QC (biallelic SNPs only, one-time per chromosome)
#   Stage 1: Chunk genome (GLIMPSE2_chunk)
#   Stage 2: Split reference panel to binary (GLIMPSE2_split_reference)
#   Stage 3: Joint-phase all samples per chunk (GLIMPSE2_phase)
#   Stage 4: Ligate chunks, RAF-filter, rename samples
#
# USAGE (array job, one task per chromosome):
#   sbatch --array=1-22 scripts/02_impute_chromosome.sh
# Or a single chromosome:
#   sbatch scripts/02_impute_chromosome.sh 11
#
# IMPORTANT: the sample list (bam_samples.txt) is regenerated from whatever
# BAMs are currently in <OUTDIR>/03_decontaminated/ every time this script
# runs. If you add or remove samples between partial runs of the SAME
# chromosome, chunks imputed before and after the change will have different
# sample counts, and GLIMPSE2_ligate will refuse to stitch them together
# ("Different number of samples in ..."). If you change the sample set,
# clear that chromosome's imputation + gl_merge checkpoints AND its partial
# chunk BCFs before rerunning. See docs/TROUBLESHOOTING.md
# ("Sample set changed between runs") — this happened during development
# and is the single most important gotcha in this pipeline.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

source "${REPO_ROOT}/config/config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CHR="${SLURM_ARRAY_TASK_ID}"
elif [ $# -ge 1 ]; then
    CHR="$1"
else
    die "No chromosome specified. Usage: sbatch --array=1-22 $0  OR  sbatch $0 <CHR>"
fi
CHECKPOINT_KEY="chr${CHR}"

set +u
source "${CONDA_SH}"
conda activate "${CONDA_ENV_IMPUTE}"
set -u
export BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS_IMPUTE}"

log "========================================"
log " 02_impute_chromosome.sh — chr${CHR}"
log " Host: $(hostname)"
log "========================================"

# =============================================================================
# PATHS
# =============================================================================
DECON_DIR="${OUTDIR}/03_decontaminated"
PANEL_FIXED_DIR="${OUTDIR}/05_glimpse_panels"
PANEL_CHUNKS_DIR="${OUTDIR}/05_glimpse_chunks"
PANEL_BINARY_DIR="${OUTDIR}/05_glimpse_binary"
GLIMPSE_DIR="${OUTDIR}/06_glimpse_imputed"
IMPUTED_DIR="${OUTDIR}/06_imputed"

mkdir -p "${PANEL_FIXED_DIR}" "${PANEL_CHUNKS_DIR}" "${PANEL_BINARY_DIR}" \
         "${GLIMPSE_DIR}" "${IMPUTED_DIR}" "${LOG_DIR}/imputation" \
         "${LOG_DIR}/gl_merge" "${CHECKPOINT_DIR}"

# =============================================================================
# SAMPLE LIST — regenerated fresh every run from 03_decontaminated/
# =============================================================================
SAMPLE_LIST="${DECON_DIR}/bam_samples.txt"
ls "${DECON_DIR}"/*.decontaminated.bam 2>/dev/null | \
    xargs -I{} basename {} .decontaminated.bam > "${SAMPLE_LIST}"
sort -u "${SAMPLE_LIST}" > "${SAMPLE_LIST}.dedup"
mapfile -t SAMPLES < "${SAMPLE_LIST}.dedup"
N_SAMPLES=${#SAMPLES[@]}
[ "${N_SAMPLES}" -gt 0 ] || die "No decontaminated BAMs found in ${DECON_DIR}"
log "Samples in current run: ${N_SAMPLES}"

# Record the sample count alongside the panel/chunk checkpoints so a
# mismatch with a partially-completed chromosome is caught loudly instead
# of surfacing later as a cryptic ligate error.
SAMPLE_COUNT_FILE="${CHECKPOINT_DIR}/chr${CHR}.sample_count.txt"
if [ -f "${SAMPLE_COUNT_FILE}" ]; then
    PREV_COUNT=$(cat "${SAMPLE_COUNT_FILE}")
    if [ "${PREV_COUNT}" != "${N_SAMPLES}" ] && checkpoint_exists "imputation" >/dev/null 2>&1; then
        die "Sample count changed for chr${CHR} (was ${PREV_COUNT}, now ${N_SAMPLES}) but chr${CHR}.imputation.done already exists. Clear chr${CHR}'s imputation/gl_merge checkpoints and partial chunk BCFs in ${GLIMPSE_DIR}/chr${CHR}_*_imputed.bcf* before rerunning. See docs/TROUBLESHOOTING.md."
    fi
fi
echo "${N_SAMPLES}" > "${SAMPLE_COUNT_FILE}"

PANEL_RAW="${PANEL_DIR}/chr${CHR}.vcf.gz"
require_file "${PANEL_RAW}" "Reference panel missing"
MAP="${GENETIC_MAP_DIR}/chr${CHR}.b38.gmap.gz"
require_file "${MAP}" "Genetic map missing"

# =============================================================================
# STAGE 0: PANEL QC (biallelic SNPs only)
# =============================================================================
if ! checkpoint_exists "panel_prep"; then
    log "STAGE 0: Panel QC [chr${CHR}]"
    PANEL_FIXED="${PANEL_FIXED_DIR}/chr${CHR}_fixed.bcf"
    LOG="${LOG_DIR}/imputation/chr${CHR}_panel.log"

    bcftools norm -m -any "${PANEL_RAW}" 2>>"${LOG}" | \
        bcftools view -m2 -M2 -v snps -O b -o "${PANEL_FIXED}" >> "${LOG}" 2>&1
    bcftools index -f "${PANEL_FIXED}" >> "${LOG}" 2>&1

    checkpoint_done "panel_prep"
fi
PANEL_FIXED="${PANEL_FIXED_DIR}/chr${CHR}_fixed.bcf"

# =============================================================================
# STAGE 1: CHUNK GENOME
# =============================================================================
if ! checkpoint_exists "chunk"; then
    log "STAGE 1: Chunk genome [chr${CHR}]"
    CHUNKS_FILE="${PANEL_CHUNKS_DIR}/chunks_chr${CHR}.txt"
    LOG="${LOG_DIR}/imputation/chr${CHR}_chunk.log"

    "${G2_CHUNK}" --input "${PANEL_FIXED}" --region "chr${CHR}" --map "${MAP}" \
        --sequential --output "${CHUNKS_FILE}" >> "${LOG}" 2>&1

    checkpoint_done "chunk"
fi
CHUNKS_FILE="${PANEL_CHUNKS_DIR}/chunks_chr${CHR}.txt"
N_CHUNKS=$(wc -l < "${CHUNKS_FILE}")
log "Chunks: ${N_CHUNKS}"

# =============================================================================
# STAGE 2: SPLIT REFERENCE PANEL TO BINARY
# =============================================================================
if ! checkpoint_exists "split_reference"; then
    log "STAGE 2: Split reference panel [chr${CHR}]"
    LOG="${LOG_DIR}/imputation/chr${CHR}_split.log"

    while IFS="" read -r LINE || [ -n "$LINE" ]; do
        printf -v ID "%06d" "$(echo "$LINE" | cut -d' ' -f1)"
        IRG=$(echo "$LINE" | cut -d' ' -f3)
        ORG=$(echo "$LINE" | cut -d' ' -f4)
        OUT_BIN="${PANEL_BINARY_DIR}/1000GP_chr${CHR}_${ID}"
        [ -f "${OUT_BIN}.bin" ] && continue

        "${G2_SPLIT}" --reference "${PANEL_FIXED}" --map "${MAP}" \
            --input-region "${IRG}" --output-region "${ORG}" \
            --output "${OUT_BIN}" --threads "${IMPUTE_CPUS}" >> "${LOG}" 2>&1
    done < "${CHUNKS_FILE}"

    checkpoint_done "split_reference"
fi

# =============================================================================
# STAGE 3: JOINT-PHASE ALL SAMPLES, PER CHUNK
# =============================================================================
if ! checkpoint_exists "imputation"; then
    log "STAGE 3: GLIMPSE2 imputation [chr${CHR}]"
    LOG="${LOG_DIR}/imputation/chr${CHR}_phase.log"

    BAM_LIST="${LOG_DIR}/imputation/bam_list_chr${CHR}.txt"
    > "${BAM_LIST}"
    for SAMPLE in "${SAMPLES[@]}"; do
        BAM="${DECON_DIR}/${SAMPLE}.decontaminated.bam"
        [ -f "${BAM}" ] && echo "${BAM}" >> "${BAM_LIST}"
    done
    log "[chr${CHR}] BAMs in phasing list: $(wc -l < "${BAM_LIST}")"

    FAILED=0
    while IFS="" read -r LINE || [ -n "$LINE" ]; do
        printf -v ID "%06d" "$(echo "$LINE" | cut -d' ' -f1)"
        BIN_REF=$(ls "${PANEL_BINARY_DIR}/1000GP_chr${CHR}_${ID}_"*.bin 2>/dev/null | head -1)
        if [ -z "${BIN_REF}" ]; then
            log "[chr${CHR}] ERROR: binary not found for chunk ${ID}"
            FAILED=$((FAILED + 1)); continue
        fi
        OUT_BCF="${GLIMPSE_DIR}/chr${CHR}_${ID}_imputed.bcf"
        [ -s "${OUT_BCF}" ] && [ -f "${OUT_BCF}.csi" ] && continue
        rm -f "${OUT_BCF}" "${OUT_BCF}.csi"

        "${G2_PHASE}" --bam-list "${BAM_LIST}" --reference "${BIN_REF}" \
            --output "${OUT_BCF}" --threads "${IMPUTE_CPUS}" >> "${LOG}" 2>&1

        if [ $? -ne 0 ] || [ ! -s "${OUT_BCF}" ]; then
            log "[chr${CHR}] ERROR: chunk ${ID} failed"
            FAILED=$((FAILED + 1))
        else
            bcftools index -f "${OUT_BCF}" >> "${LOG}" 2>&1
        fi
    done < "${CHUNKS_FILE}"

    [ "${FAILED}" -gt 0 ] && die "${FAILED} chunk(s) failed for chr${CHR} — see ${LOG}"
    checkpoint_done "imputation"
fi

# =============================================================================
# STAGE 4: LIGATE, RAF-FILTER, RENAME
# =============================================================================
if ! checkpoint_exists "gl_merge"; then
    log "STAGE 4: Ligate + filter + rename [chr${CHR}]"
    LOG="${LOG_DIR}/gl_merge/chr${CHR}.log"

    LIGATED="${GLIMPSE_DIR}/chr${CHR}_all_ligated.bcf"
    MERGED_VCF="${IMPUTED_DIR}/chr${CHR}_imputed_filtered.vcf.gz"

    LIGATE_LIST="${LOG_DIR}/imputation/ligate_chr${CHR}.txt"
    > "${LIGATE_LIST}"
    while IFS="" read -r LINE || [ -n "$LINE" ]; do
        printf -v ID "%06d" "$(echo "$LINE" | cut -d' ' -f1)"
        echo "${GLIMPSE_DIR}/chr${CHR}_${ID}_imputed.bcf" >> "${LIGATE_LIST}"
    done < "${CHUNKS_FILE}"

    rm -f "${LIGATED}" "${LIGATED}.csi"
    "${G2_LIGATE}" --input "${LIGATE_LIST}" --output "${LIGATED}" \
        --threads "${IMPUTE_CPUS}" >> "${LOG}" 2>&1

    if [ $? -ne 0 ] || [ ! -s "${LIGATED}" ]; then
        die "Ligate failed for chr${CHR} — if this mentions 'Different number of samples', see docs/TROUBLESHOOTING.md (sample set changed mid-run)"
    fi
    bcftools index -f "${LIGATED}" >> "${LOG}" 2>&1

    SAMPLE_RENAME="${LOG_DIR}/imputation/sample_rename_chr${CHR}.txt"
    bcftools query -l "${LIGATED}" | awk -F'.' '{print $0"\t"$1}' > "${SAMPLE_RENAME}"

    bcftools view -i "INFO/RAF > ${RAF_FILTER}" "${LIGATED}" 2>>"${LOG}" | \
        bcftools reheader -s "${SAMPLE_RENAME}" | \
        bcftools view -O z -o "${MERGED_VCF}" >> "${LOG}" 2>&1

    rm -f "${LIGATED}" "${LIGATED}.csi"
    require_nonempty_file "${MERGED_VCF}" "Filter/rename failed for chr${CHR}"
    bcftools index --tbi "${MERGED_VCF}" >> "${LOG}" 2>&1

    N_SITES=$(bcftools view -H "${MERGED_VCF}" 2>/dev/null | wc -l)
    log "[chr${CHR}] Done — ${N_SITES} sites, ${N_SAMPLES} samples"

    checkpoint_done "gl_merge"
fi

echo "${CHR}" >> "${IMPUTED_DIR}/chr_complete.txt"
log "========================================"
log " 02_impute_chromosome.sh complete — chr${CHR}"
log "========================================"
