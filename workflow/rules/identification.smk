checkpoint taxprofiler_samplesheet:
    input:
        'results/samplesheet.csv'
    output:
        'results/samplesheets/taxprofiler.csv'
    localrule: True
    run:
        import pandas as pd

        (
            pd.read_csv(input[0])
            
            # NOTE: Temporary
            .query("species_stdized == 'Bifidobacterium_spp'")
            .assign(
                run_accession=lambda df: df['sample'],
                
                # NOTE: Must conform to EBI ENA controlled vocabulary
                instrument_platform='ILLUMINA'
            )
            .filter(['sample', 'run_accession', 'instrument_platform', 'fastq_1', 'fastq_2'])
            .dropna()
            .sort_values('sample')
            
            # TODO: Handle better
            .drop_duplicates(subset=['fastq_1'])
            .to_csv(output[0], index=False)
        )


rule taxprofiler:
    input:
        input='results/samplesheets/taxprofiler.csv',
    output:
        'results/taxprofiler/multiqc/multiqc_report.html'
    params:
        pipeline='taxprofiler',
        profile='singularity',
        nxf='-work-dir results/taxprofiler/work -config ' + taxprofiler_cfg,
        # nxf='-work-dir results/taxprofiler/work --custom_config_base ' + taxprofiler_cfg,
        databases=workflow.source_path('../../config/databases.csv'),
        outdir='results/taxprofiler/',
    handover: True
    localrule: True
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.04.4',
        'jdk/17.0.10'
    wrapper:
        'https://raw.githubusercontent.com/fm-key-lab/snakemake-wrappers/nf-core/bio/nf-core'


# rule collect_taxprofiler:
#     input:
#         abundance_output
#     output:
#         'data/bracken.duckdb',
#     params:
#         glob="'results/bracken/*.bracken'",
#     resources:
#         cpus_per_task=8,
#         mem_mb=4_000,
#         runtime=15
#     envmodules:
#         'duckdb/nightly'
#     shell:
#         (
#             'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
#             'export GLOB={params.glob};' +
#             'duckdb -init ' + 
#             workflow.source_path('../../config/duckdbrc-slurm') +
#             ' {output} -c ".read ' + 
#             workflow.source_path('../scripts/models/bracken.sql') + 
#             '"'
#         )


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