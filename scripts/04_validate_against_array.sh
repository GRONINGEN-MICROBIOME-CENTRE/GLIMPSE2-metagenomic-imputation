#!/usr/bin/env bash
# =============================================================================
# 04_validate_against_array.sh — Validate imputed genotypes against array data
#
# Only useful if you have genotype array data for some subset of your
# samples. Computes per-sample and per-MAF-bin concordance using
# bcftools gtcheck.
#
# USAGE:
#   bash scripts/04_validate_against_array.sh <imputed_qc_vcf> <output_dir> <run_name>
#
# EXAMPLE:
#   bash scripts/04_validate_against_array.sh \
#       output/07_postqc/all_chr_imputed_QC.vcf.gz \
#       output/validation/run1 \
#       run1
#
# Set ARRAY_VCF_HG38 in config/config.sh to your hg38-lifted array VCF.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

source "${REPO_ROOT}/config/config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

if [ $# -lt 3 ]; then
    die "Usage: $0 <imputed_qc_vcf> <output_dir> <run_name>"
fi
IMPUTED_VCF="$1"
VAL_OUTDIR="$2"
RUN_NAME="$3"

set +u
source "${CONDA_SH}"
conda activate "${CONDA_ENV_IMPUTE}"
set -u
export BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS_IMPUTE}"

mkdir -p "${VAL_OUTDIR}"
require_file "${IMPUTED_VCF}" "Imputed VCF not found"
require_file "${ARRAY_VCF_HG38}" "Array VCF not found (set ARRAY_VCF_HG38 in config.sh)"

log "========================================"
log " 04_validate_against_array.sh — ${RUN_NAME}"
log "========================================"

N_IMP=$(bcftools query -l "${IMPUTED_VCF}" | wc -l)
N_ARR=$(bcftools query -l "${ARRAY_VCF_HG38}" | wc -l)
log "Imputed VCF: ${N_IMP} samples"
log "Array VCF:   ${N_ARR} samples"

# =============================================================================
# STEP 1: OVERLAPPING SITES
# =============================================================================
ARRAY_OVERLAP="${VAL_OUTDIR}/array_overlap.vcf.gz"
ARRAY_OVERLAP_RENAMED="${VAL_OUTDIR}/array_overlap_renamed.vcf.gz"
IMPUTED_OVERLAP="${VAL_OUTDIR}/imputed_overlap.vcf.gz"

if [ ! -s "${IMPUTED_OVERLAP}" ]; then
    log "STEP 1: Finding overlapping sites"
    mkdir -p "${VAL_OUTDIR}/isec_tmp"

    bcftools isec -n=2 -w1,2 --threads 8 -O z \
        -p "${VAL_OUTDIR}/isec_tmp" "${ARRAY_VCF_HG38}" "${IMPUTED_VCF}"

    mv "${VAL_OUTDIR}/isec_tmp/0000.vcf.gz" "${ARRAY_OVERLAP}"
    mv "${VAL_OUTDIR}/isec_tmp/0001.vcf.gz" "${IMPUTED_OVERLAP}"
    bcftools index --tbi "${ARRAY_OVERLAP}"
    bcftools index --tbi "${IMPUTED_OVERLAP}"
    rm -rf "${VAL_OUTDIR}/isec_tmp"

    # PLINK-style VCFs often double sample names as FamilyID_IndividualID —
    # strip back to the individual ID, keeping everything after the first
    # underscore-delimited token, e.g. FAM01_SUBJ001 -> SUBJ001.
    bcftools query -l "${ARRAY_OVERLAP}" | awk -F'_' '{print $0, $1}' > "${VAL_OUTDIR}/rename_samples.txt"
    bcftools reheader -s "${VAL_OUTDIR}/rename_samples.txt" "${ARRAY_OVERLAP}" -o "${ARRAY_OVERLAP_RENAMED}"
    bcftools index --tbi "${ARRAY_OVERLAP_RENAMED}"
else
    log "Overlap files already exist — skipping"
fi

N_OVERLAP=$(bcftools view -H "${IMPUTED_OVERLAP}" 2>/dev/null | wc -l)
N_ARR_VARS=$(bcftools view -H "${ARRAY_VCF_HG38}" 2>/dev/null | wc -l)
log "Overlapping sites: ${N_OVERLAP} ($(awk "BEGIN{printf \"%.1f\", ${N_OVERLAP}/${N_ARR_VARS}*100}")% of array sites)"

# =============================================================================
# STEP 2: GENOTYPE CONCORDANCE
# =============================================================================
CONCORDANCE_FILE="${VAL_OUTDIR}/concordance_results.txt"
if [ ! -s "${CONCORDANCE_FILE}" ]; then
    log "STEP 2: Genotype concordance (bcftools gtcheck)"
    bcftools gtcheck -g "${ARRAY_OVERLAP_RENAMED}" "${IMPUTED_OVERLAP}" 2>/dev/null | \
        grep "^DCv2" | \
        awk '$2==$3 {printf "%s\t%d\t%d\t%.4f\n", $2, $6, $7, $7/$6}' | \
        sort -k4 -n > "${CONCORDANCE_FILE}"
fi
log "Concordance computed for $(wc -l < "${CONCORDANCE_FILE}") samples"

# =============================================================================
# STEP 3: CONCORDANCE BY MAF BIN
# =============================================================================
log "STEP 3: Concordance by MAF bin"
MAF_FILE="${VAL_OUTDIR}/concordance_by_maf.txt"
echo -e "MAF_bin\tN_sites\tConcordance" > "${MAF_FILE}"

MAF_BINS=(0.01:0.05 0.05:0.10 0.10:0.20 0.20:0.30 0.30:0.50)
for BIN in "${MAF_BINS[@]}"; do
    LOW="${BIN%%:*}"; HIGH="${BIN##*:}"
    BIN_IMP="${VAL_OUTDIR}/maf_${LOW}_imp.vcf.gz"
    BIN_ARR="${VAL_OUTDIR}/maf_${LOW}_arr.vcf.gz"

    bcftools view -i "MAF[0]>=${LOW} && MAF[0]<${HIGH}" -O z -o "${BIN_IMP}" "${IMPUTED_OVERLAP}" 2>/dev/null
    bcftools index --tbi "${BIN_IMP}" 2>/dev/null
    bcftools view -i "MAF[0]>=${LOW} && MAF[0]<${HIGH}" -O z -o "${BIN_ARR}" "${ARRAY_OVERLAP_RENAMED}" 2>/dev/null
    bcftools index --tbi "${BIN_ARR}" 2>/dev/null

    N_BIN=$(bcftools view -H "${BIN_IMP}" 2>/dev/null | wc -l)
    if [ "${N_BIN}" -gt 0 ]; then
        CONC=$(bcftools gtcheck -g "${BIN_ARR}" "${BIN_IMP}" 2>/dev/null | \
            grep "^DCv2" | awk '$2==$3 {s+=$7; n+=$6} END{printf "%.4f", s/n}')
    else
        CONC="NA"
    fi
    echo -e "${LOW}-${HIGH}\t${N_BIN}\t${CONC}" >> "${MAF_FILE}"
    rm -f "${BIN_IMP}" "${BIN_IMP}.tbi" "${BIN_ARR}" "${BIN_ARR}.tbi"
done

# =============================================================================
# STEP 4: SUMMARY REPORT
# =============================================================================
REPORT="${VAL_OUTDIR}/summary_report.txt"
OVERALL=$(awk '{s+=$3; n+=$2} END{printf "%.4f", s/n}' "${CONCORDANCE_FILE}")

{
    echo "============================================================"
    echo " Imputation Validation — ${RUN_NAME}"
    echo " Generated: $(date)"
    echo "============================================================"
    echo ""
    echo "OVERLAP: ${N_OVERLAP} / ${N_ARR_VARS} array sites"
    echo ""
    echo "PER-SAMPLE CONCORDANCE (sorted, worst first)"
    printf "  %-20s  %8s  %8s  %10s\n" "Sample" "N_sites" "N_conc" "Concordance"
    awk '{printf "  %-20s  %8d  %8d  %10s\n", $1, $2, $3, $4}' "${CONCORDANCE_FILE}"
    echo "  Overall: ${OVERALL}"
    echo ""
    echo "CONCORDANCE BY MAF BIN"
    tail -n +2 "${MAF_FILE}" | awk -F'\t' '{printf "  %-15s  %8d  %10s\n", $1, $2, $3}'
    echo "============================================================"
} | tee "${REPORT}"

log "Validation complete — ${RUN_NAME}"
