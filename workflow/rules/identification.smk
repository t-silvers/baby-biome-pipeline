checkpoint taxprofiler_samplesheet:
    input:
        'results/samplesheet.csv'
    output:
        'results/samplesheets/taxprofiler.csv'
    localrule: True
    resources:
        njobs=50
    run:
        import pandas as pd

        (
            pd.read_csv(input[0])            
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
        extra='--run_kraken2 --run_bracken --run_krona --run_profile_standardisation',
        databases=workflow.source_path('../../config/databases.csv'),
        outdir='results/taxprofiler/',
    handover: True
    localrule: True
    resources:
        njobs=200
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.04.4',
        'jdk/17.0.10'
    container:
        'docker://nfcore/taxprofiler'
    wrapper:
        'https://raw.githubusercontent.com/fm-key-lab/snakemake-wrappers/nf-core/bio/nf-core'


rule collect_species_identification:
    input:
        'results/taxprofiler/multiqc/multiqc_report.html',
    output:
        multiext('results/bracken', '.duckdb', '.csv')
    params:
        glob="'results/taxprofiler/bracken/*/*.bracken.tsv'",
    resources:
        cpus_per_task=8,
        mem_mb=4_000,
        runtime=15,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export TAXPROFILER_BRACKEN={params.glob};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output[0]} -c ".read ' + 
            workflow.source_path('../scripts/models/bracken.sql') + 
            '" > {output[1]}'
        )


checkpoint reference_identification:
    input:
        'results/bracken.csv',
        'results/bracken.duckdb',
        'results/samplesheet.duckdb',
    output:
        'results/identification.csv'
    params:
        f=config['identification']['minimum_fraction_reference'],
        p=config['identification']['minimum_readspow_reference'],
    localrule: True
    resources:
        njobs=50
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export BRACKEN={input[0]};' +
            'export READ_FRAC={params.f};' +
            'export READ_POW={params.p};' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {input[2]} -c ".read ' + 
            workflow.source_path('../scripts/match_reference_genome.sql') + 
            '" > {output}'
        )