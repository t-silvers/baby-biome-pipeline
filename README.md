# Baby biome workflow

## Usage

```bash
# Requires smk>=8.26
module load snakemake/8.26.0

cd wide-variant

snakemake \
  --snakefile workflows/orchestration/Snakefile \
  --rerun-triggers mtime \
  --resources njobs=100 \
  --directory /path/to/workdir \
  --workflow-profile config/example-workflow-profile \
  --retries 5
```