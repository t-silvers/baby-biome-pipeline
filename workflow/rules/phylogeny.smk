rule filter_to_pseudogenome:
    input:
        'data/variants/species={species}/family={family}/snvs.duckdb',
    output:
        'data/pseudogenomes/species={species}/family={family}/depth={dp}/filtered_calls.parquet',
    resources:
        cpus_per_task=8,
        mem_mb=31_000,
        runtime=20,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {input} -c ".read ' + 
            workflow.source_path('../scripts/filter_calls.sql') + 
            '" > {output}'
        )


rule pseudogenome_to_msa:
    input:
        'data/pseudogenomes/species={species}/family={family}/depth={dp}/filtered_calls.parquet',
    output:
        'data/pseudogenomes/species={species}/family={family}/depth={dp}/msa.parquet',
        'results/pseudogenomes/species={species}/family={family}/depth={dp}/msa.fas',
    resources:
        cpus_per_task=8,
        mem_mb=16_000,
        runtime=5,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export FILTERED_CALLS={input};' +
            'export PSEUDOGENOMES={output[0]};' +

            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' -c ".read ' + 
            workflow.source_path('../models/msa.sql') + 
            '" > {output[0]}'

            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' -c ".read ' + 
            workflow.source_path('../scripts/calls_to_msa.sql') + 
            '" > {output[1]}'
        )


rule filter_invariant_sites:
    input:
        'results/pseudogenomes/species={species}/family={family}/depth={dp}/msa.fas',
    output:
        'results/pseudogenomes/species={species}/family={family}/depth={dp}/msa-filtered.fas',
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
        'results/pseudogenomes/species={species}/family={family}/depth={dp}/msa-filtered.fas',
    params:
        extra='--all --model GTR+G',
        outgroup='--outgroup reference',
        prefix='results/trees/species={species}/family={family}/depth={dp}',
    output:
        multiext(
            'results/trees/species={species}/family={family}/depth={dp}.raxml',
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


rule aggregate_nondepth:
    input:
        expand(
            'results/trees/species={species}/family={family}/depth={dp}.raxml.bestTree',
            dp = [5, 8, 10, 20]
        )
    output:
        'logs/smk/trees_species={species}_family={family}.done'
    shell:
        'touch {output}'


def collect_msas(wildcards):
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[1]
    )

    identification = pd.read_csv(
        checkpoints.reference_identification
        .get(**wildcards)
        .output[0]
    )

    identification = identification[
        identification['species'].isin(
            config['wildcards']['species'].split('|')
        )
    ]

    return list(
        samplesheet
        .filter(['family', 'sample'])
        .merge(identification, on='sample')
        .filter(['species', 'family'])
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: 'logs/smk/trees_species={species}_family={family}.done'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule aggregate_trees:
    input:
        collect_msas