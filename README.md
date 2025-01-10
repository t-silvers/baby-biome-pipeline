# 

## Usage

```bash
module load snakemake/8.26.0

echo "$SLURM_MAX_SUBMIT_JOB_LIMIT" # 225

ALLDATA_WORKDIR="/raven/ptmp/thosi/baby-biome/20250109-all-data"
mkdir -p $ALLDATA_WORKDIR

cd /raven/u/thosi/dev/projects/wide-variant/workflows/orchestration

snakemake \
  --rerun-triggers mtime \
  --resources njobs=$SLURM_MAX_SUBMIT_JOB_LIMIT \
  --directory "$ALLDATA_WORKDIR" \
  --workflow-profile "$GROUP_HOME/config/snakemake/profiles/widevariant" \
  --retries 2 \
  -n

```