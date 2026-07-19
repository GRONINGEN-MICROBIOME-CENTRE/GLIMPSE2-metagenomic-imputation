#!/usr/bin/env bash
#SBATCH --job-name=merge_and_qc
#SBATCH --partition=regularmedium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/slurm/merge_%j.out

# =============================================================================
# 03_merge_and_qc.sh — Concatenate chromosomes, post-imputation QC, PLINK export
#
#   Stage 1: Concatenate all per-chromosome VCFs
#   Stage 2: VCF-level QC (imputation quality, MAF, missingness)
#            Auto-detects DR2 (Beagle-style) vs INFO (GLIMPSE2) quality field.
#   Stage 3: PLINK2 bed/bim/fam
#   Stage 4: PLINK2 pgen/psam/pvar
#   Stage 5: QC report
#
# USAGE:
#   sbatch scripts/03_merge_and_qc.sh
#
# Requires all chromosomes in CHROMOSOMES (config.sh) to have completed
# scripts/02_impute_chromosome.sh first.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

source "${REPO_ROOT}/config/config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

CHECKPOINT_KEY="merge"

set +u
source "${CONDA_SH}"
conda activate "${CONDA_ENV_IMPUTE}"
set -u
export BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS_IMPUTE}"

log "========================================"
log " 03_merge_and_qc.sh"
log " Host: $(hostname)"
log "========================================"

IMPUTED_DIR="${OUTDIR}/06_imputed"
QC_DIR="${OUTDIR}/07_postqc"
mkdir -p "${QC_DIR}"

# =============================================================================
# PRE-FLIGHT: verify all chromosome VCFs exist
# =============================================================================
log "Pre-flight: verifying imputed VCFs"
MISSING=0
for CHR in "${CHROMOSOMES[@]}"; do
    VCF="${IMPUTED_DIR}/chr${CHR}_imputed_filtered.vcf.gz"
    if [ ! -s "${VCF}" ] || [ ! -f "${VCF}.tbi" ]; then
        log "  MISSING: chr${CHR}"
        MISSING=$((MISSING + 1))
    fi
done
[ "${MISSING}" -eq 0 ] || die "${MISSING} chromosome VCF(s) missing — run scripts/02_impute_chromosome.sh for those first"
log "  All ${#CHROMOSOMES[@]} chromosomes present"

# =============================================================================
# STAGE 1: MERGE ALL CHROMOSOMES
# =============================================================================
if ! checkpoint_exists "merge_chromosomes"; then
    log "STAGE 1: Merge chromosomes"
    MERGED_VCF="${IMPUTED_DIR}/all_chr_imputed.vcf.gz"
    LOG="${LOG_DIR}/merge_chromosomes.log"

    CHR_VCFS=()
    for CHR in "${CHROMOSOMES[@]}"; do
        CHR_VCFS+=("${IMPUTED_DIR}/chr${CHR}_imputed_filtered.vcf.gz")
    done

    bcftools concat --naive-force --threads "${MERGE_CPUS}" \
        -O z -o "${MERGED_VCF}" "${CHR_VCFS[@]}" > "${LOG}" 2>&1
    require_nonempty_file "${MERGED_VCF}" "Chromosome merge failed"
    bcftools index --tbi "${MERGED_VCF}" >> "${LOG}" 2>&1

    N_TOTAL=$(bcftools view -H "${MERGED_VCF}" 2>/dev/null | wc -l)
    N_SAMPLES=$(bcftools query -l "${MERGED_VCF}" 2>/dev/null | wc -l)
    log "Merged: ${N_TOTAL} variants, ${N_SAMPLES} samples"

    checkpoint_done "merge_chromosomes"
fi
MERGED_VCF="${IMPUTED_DIR}/all_chr_imputed.vcf.gz"

# =============================================================================
# STAGE 2: VCF-LEVEL QC
# Auto-detects imputation quality field: DR2 (Beagle) or INFO (GLIMPSE2).
# =============================================================================
if ! checkpoint_exists "vcf_qc"; then
    log "STAGE 2: VCF-level QC"
    QC_VCF="${QC_DIR}/all_chr_imputed_QC.vcf.gz"
    LOG="${LOG_DIR}/postqc_filter.log"

    HAS_DR2=$(bcftools view -h "${MERGED_VCF}" | grep -c "##INFO.*ID=DR2" || true)
    HAS_INFO=$(bcftools view -h "${MERGED_VCF}" | grep -c "##INFO.*ID=INFO" || true)

    if [ "${HAS_DR2}" -gt 0 ]; then
        IMP_FILTER="INFO/DR2 >= ${IMPUTATION_QUALITY_THRESHOLD}"
        log "  Imputation quality field: DR2 (Beagle-style output)"
    elif [ "${HAS_INFO}" -gt 0 ]; then
        IMP_FILTER="INFO/INFO >= ${IMPUTATION_QUALITY_THRESHOLD}"
        log "  Imputation quality field: INFO (GLIMPSE2 output)"
    else
        IMP_FILTER=""
        log "  WARNING: no imputation quality field found — skipping quality filter"
    fi

    if [ -n "${IMP_FILTER}" ]; then
        bcftools filter --include "${IMP_FILTER}" "${MERGED_VCF}" 2>>"${LOG}" | \
            bcftools filter --include "MAF[0] >= ${MAF_FILTER}" 2>>"${LOG}" | \
            bcftools filter --include "F_MISSING < ${MAX_SITE_MISSING}" \
                --threads "${MERGE_CPUS}" -O z -o "${QC_VCF}" >> "${LOG}" 2>&1
    else
        bcftools filter --include "MAF[0] >= ${MAF_FILTER}" "${MERGED_VCF}" 2>>"${LOG}" | \
            bcftools filter --include "F_MISSING < ${MAX_SITE_MISSING}" \
                --threads "${MERGE_CPUS}" -O z -o "${QC_VCF}" >> "${LOG}" 2>&1
    fi

    require_nonempty_file "${QC_VCF}" "VCF QC filtering failed"
    bcftools index --tbi "${QC_VCF}" >> "${LOG}" 2>&1

    checkpoint_done "vcf_qc"
fi
QC_VCF="${QC_DIR}/all_chr_imputed_QC.vcf.gz"

# =============================================================================
# STAGE 3: PLINK BED
# =============================================================================
if ! checkpoint_exists "plink_bed"; then
    log "STAGE 3: PLINK2 bed/bim/fam"
    QC_PLINK="${QC_DIR}/genotypes_QC"
    LOG="${LOG_DIR}/postqc_plink.log"

    plink2 --vcf "${QC_VCF}" --make-bed --out "${QC_PLINK}" \
        --max-alleles 2 --mind "${MAX_SAMPLE_MISSING}" --geno "${MAX_SITE_MISSING}" \
        --maf "${MAF_FILTER}" --threads "${MERGE_CPUS}" >> "${LOG}" 2>&1

    checkpoint_done "plink_bed"
fi

# =============================================================================
# STAGE 4: PLINK PGEN
# =============================================================================
if ! checkpoint_exists "plink_pgen"; then
    log "STAGE 4: PLINK2 pgen"
    QC_PGEN="${QC_DIR}/genotypes_QC_pgen"
    LOG="${LOG_DIR}/postqc_pgen.log"

    plink2 --vcf "${QC_VCF}" --make-pgen --out "${QC_PGEN}" \
        --max-alleles 2 --mind "${MAX_SAMPLE_MISSING}" --geno "${MAX_SITE_MISSING}" \
        --maf "${MAF_FILTER}" --threads "${MERGE_CPUS}" >> "${LOG}" 2>&1 \
        || log "WARNING: pgen conversion failed — bed format is still available"

    checkpoint_done "plink_pgen"
fi

# =============================================================================
# STAGE 5: QC REPORT
# =============================================================================
if ! checkpoint_exists "qc_report"; then
    log "STAGE 5: QC report"
    QC_PLINK="${QC_DIR}/genotypes_QC"
    REPORT="${QC_DIR}/QC_report.txt"

    N_INPUT=$(bcftools view -H "${MERGED_VCF}" 2>/dev/null | wc -l)
    N_SAMPLES_IN=$(bcftools query -l "${MERGED_VCF}" 2>/dev/null | wc -l)
    N_QC=$(bcftools view -H "${QC_VCF}" 2>/dev/null | wc -l)
    N_BIM=$([ -f "${QC_PLINK}.bim" ] && wc -l < "${QC_PLINK}.bim" || echo 0)

    {
        echo "============================================================"
        echo " Post-Imputation QC Report"
        echo " Generated: $(date)"
        echo "============================================================"
        echo ""
        echo "PIPELINE: Trimmomatic -> Bowtie2 (human) -> Picard MarkDuplicates"
        echo "          -> Bowtie2 (bacterial decontamination) -> Joint GLIMPSE2"
        echo ""
        echo "COHORT: ${N_SAMPLES_IN} samples"
        echo ""
        echo "QC THRESHOLDS"
        echo "  Imputation quality >= ${IMPUTATION_QUALITY_THRESHOLD} (DR2 or INFO, auto-detected)"
        echo "  MAF >= ${MAF_FILTER}"
        echo "  Site missingness < ${MAX_SITE_MISSING}"
        echo "  Sample missingness < ${MAX_SAMPLE_MISSING}"
        echo ""
        echo "VARIANT COUNTS"
        printf "  %-30s %10d\n" "After imputation (merged):" "${N_INPUT}"
        printf "  %-30s %10d\n" "After VCF QC:" "${N_QC}"
        printf "  %-30s %10d\n" "After PLINK QC:" "${N_BIM}"
        echo ""
        echo "OUTPUT FILES"
        echo "  Merged VCF:  ${MERGED_VCF}"
        echo "  QC VCF:      ${QC_VCF}"
        echo "  PLINK bed:   ${QC_PLINK}.bed/.bim/.fam"
        echo "============================================================"
    } | tee "${REPORT}"

    checkpoint_done "qc_report"
fi

log "========================================"
log " 03_merge_and_qc.sh complete"
log "========================================"
