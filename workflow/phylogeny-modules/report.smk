rule plot_trees:
    input:
        workflow.source_path('results/trees/{cohort}/msa-filtered.raxml.bestTree'),
        workflow.source_path('results/trees/{cohort}/metadata.txt'),
    output:
        workflow.source_path('results/trees/{cohort}/phylo.png'),
    params:
        script=workflow.source_path('../scripts/plot_tree.R'),
    resources:
        cpus_per_task=2,
        njobs=50
    envmodules:
        'R/4.4'
    shell:
        'Rscript {params.script} -t {input[0]} -m {input[1]} -o {output}'


def aggregate_tree_plots(wildcards):
    return expand(
        workflow.source_path('results/trees/{cohort}/phylo.png'),
        cohort=list(cohort_to_glob.keys())
    )


rule all_tree_plots:
    input:
        aggregate_tree_plots
