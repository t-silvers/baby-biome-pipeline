checkpoint check_msas:
    input:
        collect_msas
    output:
        directory('results/trees')
    localrule: True
    resources:
        njobs=50
    shell:
        '''
        mkdir -p results/trees
        for i in results/pseudogenomes/species=*_family=*_msa.fas; do
          if [ -s "$i" ]; then
            cp "$i" "results/trees/"
          fi
        done
        '''


rule filter_invariant_sites:
    input:
        'results/trees/{cohort}_msa.fas'
    output:
        'results/trees/{cohort}_msa_filtered.fas',
    resources:
        njobs=1
    envmodules:
        'snp-sites/2.5.1'
    shell:
        'snp-sites -o {output} {input}'


def get_tree_cpus(wildcards, attempt):
    return attempt * 2 + 2


def get_tree_mem(wildcards, attempt):
    return attempt * 16_000


def get_tree_time(wildcards, attempt):
    return attempt * 30


rule raxml_ng:
    input:
        'results/trees/{cohort}_msa_filtered.fas',
    params:
        extra='--all --model GTR+G --bs-trees 1000',
        outgroup='--outgroup reference',
        prefix='results/trees/{cohort}',
    output:
        multiext(
            'results/trees/{cohort}.raxml',
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


def collect_trees(wildcards):
    import os

    tree_dir = checkpoints.check_msas.get(**wildcards).output[0]

    wc = glob_wildcards(
        os.path.join(tree_dir, '{cohort}_msa.fas')
    )

    return expand(
        'results/trees/{cohort}.raxml.bestTree',
        cohort=wc.cohort
    )


rule:
    input:
        collect_trees
    output:
        'results/trees.done'
    localrule: True
    resources:
        njobs=50
    shell:
        'touch {output}'