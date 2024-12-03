rule:
    input:
        'results/variants/species={species}/family={family}/pseudogenomes.csv',
        'results/samplesheet.duckdb',
    output:
        'results/variants/species={species}/family={family}/msa.fas',
    resources:
        cpus_per_task=12,
        mem_mb=16_000,
        runtime=5,
    envmodules:
        'duckdb/nightly'
    # TODO: Add reference genome to MSA
    shell:        
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export FILTERED_CALLS={input[0]};' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {input[1]} -c ".read ' + 
            workflow.source_path('../scripts/calls_to_msa.sql') + 
            '" > {output}'
        )


rule filter_invariant_sites:
    input:
        'results/variants/species={species}/family={family}/msa.fas',
    output:
        'results/variants/species={species}/family={family}/msa_filtered.fas',
    envmodules:
        'snp-sites/2.5.1'
    shell:
        'snp-sites -o {output} {input}'


rule raxml_ng:
    input:
        'results/variants/species={species}/family={family}/msa_filtered.fas',
    params:
        extra='--all --model GTR+G --bs-trees 1000',
        # outgroup=lambda wildcards: '--outgroup ' + config['outgroup'][wildcards.donor][wildcards.species]['ID'],
        outgroup='',
        prefix='results/variants/species={species}/family={family}/msa_filtered',
    output:
        multiext(
            'results/variants/species={species}/family={family}/msa_filtered.raxml',
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
        cpus_per_task=48,
        mem_mb=64_000,
        runtime=120,
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
          --redo
        
        touch {output}
        '''


def collect_trees(wildcards):
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

    return list(
        identification
        [identification['species'].isin(config['wildcards']['species'].split('|'))]
        .merge(
            samplesheet.filter(['family', 'sample']),
            on='sample'
        )
        .filter(['species', 'family'])
        .drop_duplicates()
        .transpose()
        .apply(lambda df: 'results/variants/species={species}/family={family}/msa_filtered.raxml.bestTree'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule:
    input:
        collect_trees
    output:
        temp('results/trees.done')
    localrule: True
    shell:
        'touch {output}'