## Usage

from [here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources):

> In general, resources are just names to the Snakemake scheduler, i.e., Snakemake does not check on the resource consumption of jobs in real time. Instead, resources are used to determine which jobs can be executed at the same time without exceeding the limits specified at the command line. Apart from making Snakemake aware of hybrid-computing architectures (e.g. with a limited number of additional devices like GPUs) this allows us to control scheduling in various ways, e.g. to limit IO-heavy jobs by assigning an artificial IO-resource to them and limiting it via the --resources flag. If no limits are given, the resources are ignored in local execution.

```bash
ml task
labtasks smk PROJ="bifido-multilib" EXTRA="--rerun-triggers mtime --resources njobs=200"
```