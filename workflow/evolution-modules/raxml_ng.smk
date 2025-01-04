def get_tree_cpus(wildcards, attempt):
    return attempt * 0 + 8


def get_tree_mem(wildcards, attempt):
    return attempt * 15_000


def get_tree_time(wildcards, attempt):
    return attempt * 60


rule raxml_ng:
    input:
        'results/trees/{cohort}/msa-filtered.fas',
    output:
        multiext(
            'results/trees/{cohort}/msa-filtered.raxml',
            '.reduced.phy',
            '.rba',
            '.bestTreeCollapsed',
            '.bestTree',
            '.mlTrees',
            '.support',
            '.bestModel',
            '.bootstraps',
            '.log'
        )
    params:
        extra='--all --model GTR+G --bs-trees 10',
        outgroup='--outgroup reference',
        prefix='results/trees/{cohort}/msa-filtered',
    resources:
        cpus_per_task=get_tree_cpus,
        mem_mb=get_tree_mem,
        runtime=get_tree_time,
        njobs=1
    envmodules:
        'raxml-ng/1.2.2_MPI'
    shell:
        '''
        export OMP_PLACES=threads

        raxml-ng \
          {params.extra} \
          {params.outgroup} \
          --msa {input} \
          --threads {resources.cpus_per_task} \
          --prefix {params.prefix} \
          --blopt nr_safe \
          --redo
        
        touch {output}
        '''
