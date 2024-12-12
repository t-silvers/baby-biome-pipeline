rule allele_depth_plot_data:
    input:
        'data/variants/species={species}/family={family}/snvs.duckdb',
        'results/pseudogenomes/species={species}/family={family}/filtered_calls.parquet',
    output:
        'reports/variants/species={species}/family={family}/allele_depth_plot_data.csv',
    resources:
        cpus_per_task=8,
        mem_mb=48_000,
        runtime=10,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export FILTERED_CALLS={input[1]};' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {input[0]} -c ".read ' + 
            workflow.source_path('../scripts/data_for_dp_plots.sql') + 
            '" > {output}'
        )


def collect_allele_depth_plot_data(wildcards):
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
        .apply(lambda df: 'reports/variants/species={species}/family={family}/allele_depth_plot_data.csv'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule:
    input:
        collect_allele_depth_plot_data
    output:
        'results/report.done'
    localrule: True
    resources:
        njobs=50
    shell:
        'touch {output}'