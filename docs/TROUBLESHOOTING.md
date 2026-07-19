# Troubleshooting

Real issues encountered running this pipeline, kept here so they don't have
to be re-diagnosed from scratch.

## FASTQ naming mismatch

**Symptom:** `sbatch scripts/01_preprocess_sample.sh SAMPLE` finishes in a
few seconds, with `ERROR: R1/R2 not found` in the log.

**Cause:** the script looks for `<SAMPLE>_R1.fastq.gz` / `_R2.fastq.gz`.
Some sequencing providers deliver `<SAMPLE>_1.fastq.gz` / `_2.fastq.gz`
instead.

**Fix:** check what's actually there, then rename or symlink to match:
```bash
ls FASTQ_DIR/<SAMPLE>/
# if you see _1.fastq.gz / _2.fastq.gz instead of _R1/_R2:
mv <SAMPLE>_1.fastq.gz <SAMPLE>_R1.fastq.gz
mv <SAMPLE>_2.fastq.gz <SAMPLE>_R2.fastq.gz
```
A loop for a whole batch:
```bash
for d in FASTQ_DIR/*/; do
    SAMPLE=$(basename "$d")
    F1="${d}${SAMPLE}_1.fastq.gz"; F2="${d}${SAMPLE}_2.fastq.gz"
    if [ -f "$F1" ] && [ ! -f "${d}${SAMPLE}_R1.fastq.gz" ]; then
        mv -v "$F1" "${d}${SAMPLE}_R1.fastq.gz"
        mv -v "$F2" "${d}${SAMPLE}_R2.fastq.gz"
    fi
done
```

## Truncated trimmed FASTQ

**Symptom:** Stage 3 (Bowtie2 alignment) dies partway through with:
```
Error, fewer reads in file specified with -1 than in file specified with -2
terminate called after throwing an instance of 'int'
(ERR): bowtie2-align died with signal 6 (ABRT)
```
even though Stage 2 (Trimmomatic) reported "Already trimmed — skipping".

**Cause:** an earlier, interrupted run of Stage 2 left a partially-written
(but non-empty) trimmed FASTQ behind. The checkpoint's existence check only
tests `[ -f ] && [ -s ]` — file exists and is non-empty — not that it's
complete. A truncated gzip stream passes that check.

**Fix (built into this version of the pipeline):**
`scripts/01_preprocess_sample.sh` now calls
`require_paired_fastq_counts` after trimming, which compares R1/R2 read
counts and fails loudly if they differ, instead of letting the mismatch
surface hours later mid-alignment. If you hit this on an older checkpoint,
clear it and rerun:
```bash
rm -f output/logs/checkpoints/<SAMPLE>.stage2_trimmomatic.done
rm -f output/01_trimmed/<SAMPLE>_R1_paired.fastq.gz output/01_trimmed/<SAMPLE>_R2_paired.fastq.gz
sbatch scripts/01_preprocess_sample.sh <SAMPLE>
```

## Sample set changed between runs

**This is the most important gotcha in the pipeline — read this before
adding or removing samples mid-project.**

**Symptom:** `scripts/02_impute_chromosome.sh` fails at the ligate step
with:
```
ERROR: Different number of samples in chr<N>_<chunk>_imputed.bcf.
```

**Cause:** `bam_samples.txt` is regenerated from whatever BAMs are
currently in `output/03_decontaminated/` every time Stage 2 runs. Stage 3
(the chunk-imputation loop) skips any chunk whose output BCF already
exists — it only checks file presence, not which sample set produced it.
If a chromosome's chunk loop is interrupted partway through, and you then
add or remove a sample from `03_decontaminated/` before resubmitting, the
chunks completed *before* the change and the chunks completed *after* will
have different sample counts. `GLIMPSE2_ligate` correctly refuses to stitch
them together.

This is exactly what happened during development: a sample with
effectively zero usable reads (0.0003x depth) was excluded partway through
a chromosome's imputation run, and the chunks already done before the
exclusion silently carried the old sample count forward.

**Diagnosis:** check chunk BCF timestamps for a clean before/after split:
```bash
ls -la output/06_glimpse_imputed/chr<N>_*_imputed.bcf | sort -k6,7
```
A jump in timestamp and a jump in file size at the same point is the
signature of this issue.

**Fix — safest option, always correct:** wipe and redo Stage 3 + 4 for the
affected chromosome(s):
```bash
CHR=<N>
rm -f output/logs/checkpoints/chr${CHR}.imputation.done
rm -f output/logs/checkpoints/chr${CHR}.gl_merge.done
rm -f output/06_glimpse_imputed/chr${CHR}_*_imputed.bcf*
rm -f output/06_glimpse_imputed/chr${CHR}_all_ligated.bcf*
sbatch scripts/02_impute_chromosome.sh ${CHR}
```

**Fix — surgical option (faster, only if you're confident about the
split):** remove only the chunks computed before the sample-set change,
identified by timestamp, then rerun — the loop will redo just those and
reuse the rest.

**Built-in guard:** this version of the pipeline records the sample count
alongside each chromosome's checkpoints
(`output/logs/checkpoints/chr<N>.sample_count.txt`) and refuses to proceed
past an already-completed imputation checkpoint if the current sample
count doesn't match, printing a pointer back to this section instead of
letting it fail later at ligate with a more cryptic message.

## Excluding a sample

To exclude a sample from imputation entirely (e.g. it's below your depth
threshold, or turns out to be unmapped/unidentifiable), don't just filter
`bam_samples.txt` — it's regenerated from the directory contents every run
and your edit will be overwritten. Move the BAM out of the input directory
instead:
```bash
mkdir -p output/03_decontaminated/excluded
mv output/03_decontaminated/<SAMPLE>.decontaminated.bam* output/03_decontaminated/excluded/
```
Then follow the "sample set changed" fix above for any chromosome that was
partway through imputation when you did this.

## Out-of-memory on large chromosomes

**Symptom:** `sacct` shows `OUT_OF_MEMORY` (or a chromosome's phase log
just stops mid-chunk with no error at all — see next section) for large
chromosomes (chr1–chr3 especially), while small chromosomes complete fine.

**Cause:** GLIMPSE2's joint-phasing memory footprint scales with sample
count and chunk density. Default resource requests
(`IMPUTE_MEM=99G` in `config.sh`) were tuned for one cohort size; a
different project with more (or differently-distributed) samples can
exceed that on the largest chromosomes even if the total sample count is
smaller than a previous successful run.

**Fix:** bump memory for the affected chromosomes:
```bash
sbatch --mem=128G --array=1,2,3 scripts/02_impute_chromosome.sh
```

## Job died with no error and no SLURM record

**Symptom:** a chromosome's phase log stops mid-sentence, no error, no
"Killed" message — and `sacct`/`squeue` show no record of the job at all.

**Cause:** almost always means the script was run **interactively**
(`ssh`'d into a node, run in a terminal/tmux session directly) rather than
via `sbatch`, and the session was dropped (network blip, terminal closed,
node recycled). Interactive runs get no SLURM accounting.

**Fix:** always submit via `sbatch --array`, even for a single chromosome
you're rerunning:
```bash
sbatch --array=<N> scripts/02_impute_chromosome.sh
```
This guarantees a real exit code and `MaxRSS`/state in `sacct` if it fails
again, instead of another silent stall.

## Stale/leftover checkpoint from an earlier, unrelated run

**Symptom:** a chromosome's checkpoint (e.g. `chr4.imputation.done`) exists,
but there's no corresponding chunk data or final VCF, and the run behaves
as if that stage were skipped incorrectly.

**Cause:** checkpoints from an earlier pipeline run (different dataset,
earlier test, or a since-cleaned-up run) can be left behind in
`output/logs/checkpoints/` if the working directory is reused across
projects.

**Fix:** before starting a new project or sample batch, confirm the
checkpoint directory doesn't have stale entries from a previous run:
```bash
ls output/logs/checkpoints/ | grep "^chr<N>\."
```
If a checkpoint has no matching output, just remove it and rerun that
stage. Consider using a separate `OUTDIR` per project in `config.sh` to
avoid this entirely.
