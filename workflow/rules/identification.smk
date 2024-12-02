def paired_fastqs(wildcards):
    import pandas as pd

    return (
        pd.read_csv(
            checkpoints.fastqs.get(**wildcards).output[0]
        )
        .dropna()
        .query(
            f"sample == {wildcards['sample']}"
        )
        .filter(like='fastq', axis=1)
        .values
        .flatten()
    )


rule kraken2:
    input:
        paired_fastqs
    output:
        'results/kraken2/{sample}.kraken2',
        'results/kraken2/{sample}.k2report',
    params:
        extra='--paired --report-minimizer-data --minimum-hit-groups 3',
        k2db=config['public_data']['kraken_db'],
    resources:
        cpus_per_task=4,
        mem_mb=128_000,
        runtime=5,
    envmodules:
        'intel/2025.0',
        'impi/2021.14',
        'kraken2/2.1.3'
    shell:
        '''
        export OMP_PLACES=threads
        
        kraken2 {params.extra} --db {params.k2db} --threads {resources.cpus_per_task} --output {output[0]} --report {output[1]} {input}
        '''


rule bracken:
    input:
        'results/kraken2/{sample}.kraken2',
        'results/kraken2/{sample}.k2report',
    output:
        'results/bracken/{sample}.bracken',
        'results/bracken/{sample}.breport',
    params:
        k2db=config['public_data']['kraken_db']
    resources:
        runtime=5,
    envmodules:
        'bracken/2.9'
    shell:
        'bracken -d {params.k2db} -r 100 -i {input[1]} -o {output[0]} -w {output[1]}'


def abundance_output(wildcards):
    import pandas as pd

    sample_ids = (
        pd.read_csv(
            checkpoints.fastqs.get(**wildcards).output[0]
        )
        # NOTE: Not all FASTQs can be resolved to IDs
        .dropna()
        ['sample'].astype(str)
    )

    return expand(
        [
            'results/bracken/{sample}.bracken',
            'results/bracken/{sample}.breport',
        ],
        sample=sample_ids
    )


rule collect_bracken:
    input:
        abundance_output
    output:
        'data/bracken.duckdb',
    params:
        glob="'results/bracken/*.bracken'",
    resources:
        cpus_per_task=8,
        mem_mb=4_000,
        runtime=15
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export GLOB={params.glob};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output} -c ".read ' + 
            workflow.source_path('../scripts/models/bracken.sql') + 
            '"'
        )


# rule reference_identification:
#     input:
#         'data/bracken.duckdb',
#         'data/sample_info.duckdb',
#     output:
#         'results/identification.csv'
#     params:
#         f=config['identification']['minimum_fraction_reference'],
#         p=config['identification']['minimum_readspow_reference'],
#     localrule: True
#     envmodules:
#         'duckdb/nightly'
#     shell:
#         '''
#         export READ_FRAC={params.f} \
#                READ_POW={params.p}

#         # TODO: Must update
#         duckdb -readonly -init config/duckdbrc-local {input} -c '.read workflow/scripts/match_reference_genome.sql' > {output}
#         '''