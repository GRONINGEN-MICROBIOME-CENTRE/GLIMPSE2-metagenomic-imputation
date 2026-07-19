# Usage guide

## Setup

```bash
cp config/config.example.sh config/config.sh
```

Edit `config/config.sh` to point at your cluster's paths: conda environments,
reference genome, Bowtie2 indices, GLIMPSE2 binaries, reference panel, and
genetic maps. Every script sources this one file, so a path only needs to be
correct in one place.

Input FASTQs are expected as `<SAMPLE_ID>_R1.fastq.gz` /
`<SAMPLE_ID>_R2.fastq.gz`, either directly under `FASTQ_DIR` or in a
per-sample subdirectory `FASTQ_DIR/<SAMPLE_ID>/`. If your sequencing
provider delivers a different naming convention (e.g. `_1`/`_2`), see
[TROUBLESHOOTING.md](TROUBLESHOOTING.md#fastq-naming-mismatch).

## Stage 1 — Preprocessing (`scripts/01_preprocess_sample.sh`)

Runs per sample: Trimmomatic → Bowtie2 (human) → Picard MarkDuplicates →
Bowtie2 (bacterial decontamination). Output:
`output/03_decontaminated/<SAMPLE>.decontaminated.bam`.

```bash
sbatch scripts/01_preprocess_sample.sh SAMPLE001
# or for a whole list:
bash submit/submit_all_samples.sh my_sample_list.txt
```

### Depth QC

After decontamination, mean depth is written to
`output/logs/decon_summary.tsv`. Samples below `MIN_DEPTH_X` (default 0.5x
— set in `config.sh`) are **flagged**, not silently dropped, in
`output/logs/depth_flags.tsv`. Below ~0.25x, concordance drops sharply in
our validation; between 0.25–0.5x it's usable but degraded (~97% vs ~99%
concordance). Decide per-project whether to exclude flagged samples before
Stage 2 — see [TROUBLESHOOTING.md](TROUBLESHOOTING.md#excluding-a-sample).

## Stage 2 — Joint imputation (`scripts/02_impute_chromosome.sh`)

One SLURM task per chromosome, submitted as an array job:

```bash
sbatch --array=1-22 scripts/02_impute_chromosome.sh
```

All samples currently in `output/03_decontaminated/` are imputed **jointly**
per chunk (GLIMPSE2's joint-calling mode), which is what gives it an edge
over per-sample imputation at low coverage — samples with more data help
inform genotype calls at low-coverage sites in other samples.

**Do not run this interactively** (`ssh` into a node and run the script by
hand) for anything beyond a quick single-chunk test. Interactive sessions
can be dropped silently (network blip, terminal closed) with no SLURM
accounting record and no error in the log — just a job that stops mid-chunk
with nothing to tell you why. Always submit via `sbatch`/`--array` so you
get a real exit code and can check `sacct`.

## Stage 3 — Merge & QC (`scripts/03_merge_and_qc.sh`)

```bash
sbatch scripts/03_merge_and_qc.sh
```

Concatenates all 22 chromosomes, applies quality/MAF/missingness filters
(auto-detecting whether the VCF has a Beagle-style `DR2` field or a
GLIMPSE2-style `INFO` field), and exports PLINK bed/pgen. Produces
`output/07_postqc/QC_report.txt`.

## Stage 4 — Validation against array data (optional)

If you have genotype array data (hg38) for a subset of samples:

```bash
bash scripts/04_validate_against_array.sh \
    output/07_postqc/all_chr_imputed_QC.vcf.gz \
    output/validation/run1 run1
```

Produces per-sample and per-MAF-bin concordance against `ARRAY_VCF_HG38`
(set in `config.sh`).

## Resuming

Every stage checkpoints its progress under `output/logs/checkpoints/`:

- Stage 1: `<SAMPLE>.<stage_name>.done`
- Stage 2: `chr<N>.<stage_name>.done`
- Stage 3: `merge.<stage_name>.done`

Re-running any script skips steps whose checkpoint already exists. This
makes it safe to resubmit after a failed or interrupted job — **as long as
your inputs haven't changed**. If you add or remove samples, see
[TROUBLESHOOTING.md](TROUBLESHOOTING.md#sample-set-changed-between-runs)
before resubmitting Stage 2 for any chromosome that was partway through.

## Adjusting resources

Default SLURM resource requests are in the `#SBATCH` headers of each script
and mirrored in `config/config.sh` (`PREPROCESS_*`, `IMPUTE_*`, `MERGE_*`).
Override at submission time without editing the script, e.g.:

```bash
sbatch --mem=128G --array=1-22 scripts/02_impute_chromosome.sh
```

Large/dense chromosomes (chr1–chr3 especially) are the most likely to need
more memory as your sample count grows — see
[TROUBLESHOOTING.md](TROUBLESHOOTING.md#out-of-memory-on-large-chromosomes).
