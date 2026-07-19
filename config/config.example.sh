#!/usr/bin/env bash
# =============================================================================
# config.example.sh — Central configuration for the pipeline
#
# Copy this file to config/config.sh and fill in the values for your cluster
# and dataset. config/config.sh is gitignored — never commit real paths,
# usernames, or credentials.
#
#   cp config/config.example.sh config/config.sh
#   # then edit config/config.sh
#
# Every script in scripts/ sources this file first. Changing a path here
# changes it everywhere — there should be no hardcoded paths inside the
# stage scripts themselves.
#
# WHY THIS INDIRECTION EXISTS (read this if you're new to the repo):
# You will NOT find the literal string "host_preprocessing" or
# "glimpse_imputation" anywhere inside scripts/*.sh — the scripts reference
# variables (CONDA_ENV_ALIGN, CONDA_ENV_IMPUTE, etc.) whose actual values
# live ONLY here, in config.sh. This means the same scripts work for anyone
# who clones this repo, no matter what they've named their own conda
# environments — they only ever edit this one file. To find what
# environment a script actually uses, grep THIS file, not the scripts:
#   grep CONDA_ENV_ config/config.sh
# Every script also prints its resolved environment name at startup
# (e.g. "Conda environment: host_preprocessing") so you can always confirm
# what actually ran from the log, without reading any source code.
# =============================================================================

# ---- Conda ---------------------------------------------------------------
CONDA_SH="/path/to/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_ALIGN="host_preprocessing"   # env with trimmomatic, bowtie2, picard, samtools
CONDA_ENV_IMPUTE="glimpse_imputation"  # env with GLIMPSE2, bcftools, plink2

# ---- Base directories ------------------------------------------------------
BASE_DIR="/path/to/project/genotype_imputation"
FASTQ_DIR="${BASE_DIR}/data/samples"          # one subdir per sample, or flat
OUTDIR="${BASE_DIR}/output"
LOG_DIR="${OUTDIR}/logs"
CHECKPOINT_DIR="${LOG_DIR}/checkpoints"

# ---- Reference data ---------------------------------------------------------
REF_GENOME="${BASE_DIR}/reference_genome/GRCh38_no_alt_plus_hs38d1.fa"
HUMAN_BT2_IDX="${BASE_DIR}/reference_genome/bowtie_index/GRCh38_no_alt_plus_hs38d1"
BACTERIAL_BT2_IDX="${BASE_DIR}/databases/bacterial_refseq/bacteria_index"
# Path to the Trimmomatic adapter FASTA shipped inside your CONDA_ENV_ALIGN
# environment, e.g.:
#   /path/to/miniconda3/envs/host_preprocessing/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
ADAPTERS_FASTA="/path/to/miniconda3/envs/${CONDA_ENV_ALIGN}/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa"

# 1000G / HGDP high-coverage GRCh38 panel (per-chromosome VCFs: chr{1..22}.vcf.gz)
PANEL_DIR="${BASE_DIR}/databases/HGDP_1000G_GRCh38"
GENETIC_MAP_DIR="${BASE_DIR}/databases/genetic_maps_glimpse2"

# GLIMPSE2 static binaries
GLIMPSE2_DIR="${BASE_DIR}/tools/glimpse2"
G2_CHUNK="${GLIMPSE2_DIR}/GLIMPSE2_chunk_static"
G2_SPLIT="${GLIMPSE2_DIR}/GLIMPSE2_split_reference_static"
G2_PHASE="${GLIMPSE2_DIR}/GLIMPSE2_phase_static"
G2_LIGATE="${GLIMPSE2_DIR}/GLIMPSE2_ligate_static"

BCFTOOLS_PLUGINS_ALIGN="/path/to/miniconda3/envs/${CONDA_ENV_ALIGN}/libexec/bcftools"
BCFTOOLS_PLUGINS_IMPUTE="/path/to/miniconda3/envs/${CONDA_ENV_IMPUTE}/libexec/bcftools"

# ---- Array genotype validation data (optional — only needed for scripts/04) ---
ARRAY_VCF_HG38="${BASE_DIR}/genotypes/final_QCed_genotypes/array_hg38.vcf.gz"

# ---- Resource defaults (override at sbatch submission time if needed) --------
# See submit/ for the sbatch wrappers that apply these.
PREPROCESS_CPUS=24
PREPROCESS_MEM="99G"
PREPROCESS_TIME="3-00:00:00"

IMPUTE_CPUS=60
IMPUTE_MEM="99G"
IMPUTE_TIME="10-00:00:00"

MERGE_CPUS=16
MERGE_MEM="64G"
MERGE_TIME="3-00:00:00"

# ---- Alignment / QC thresholds ------------------------------------------------
MIN_MAPQ=20
MIN_BASEQ=20

MAF_FILTER=0.01
MAX_SITE_MISSING=0.05
MAX_SAMPLE_MISSING=0.05
IMPUTATION_QUALITY_THRESHOLD=0.8   # applied to DR2 (Beagle) or INFO (GLIMPSE2), auto-detected
RAF_FILTER=0.001                   # minimum reference allele frequency kept after ligation

MIN_DEPTH_X=0.5                    # samples below this depth are flagged, not silently dropped

# ---- Chromosomes to process --------------------------------------------------
CHROMOSOMES=($(seq 1 22))
