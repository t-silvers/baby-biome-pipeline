checkpoint taxprofiler_samplesheet:
    input:
        'resources/samplesheets/main.csv'
    output:
        'resources/samplesheets/taxprofiler.csv'
    log:
        'logs/smk/identification/taxprofiler_samplesheet.log'
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
        input='resources/samplesheets/taxprofiler.csv',
    output:
        'results/taxprofiler/multiqc/multiqc_report.html'
    params:
        databases=config['identification']['taxprofiler']['databases'],
        extra=config['identification']['taxprofiler']['args'],
        nxf=config['identification']['taxprofiler']['nxf_args'],
        nxf_log='logs/nxf/taxprofiler.log',
        outdir='results/taxprofiler/',
        pipeline='taxprofiler',
        profile=config['identification']['taxprofiler']['profiles']
    log:
        'logs/smk/identification/taxprofiler.log'
    handover: True
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


rule taxprofiler_bracken:
    input:
        ancient('results/taxprofiler/multiqc/multiqc_report.html'),
    output:
        'data/bracken.duckdb',
        'results/taxprofiler/bracken/summary.csv',
    params:
        glob="'results/taxprofiler/bracken/*/*.bracken.tsv'",
        pat="'results/taxprofiler/bracken/\S*/(\d+)_\d+_\S*.bracken.tsv$'",
        model=config['models']['bracken']
    log:
        'logs/smk/identification/taxprofiler_bracken.log'
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
            'export TAXPROFILER_PAT={params.pat};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output[0]} -c ".read {params.model}" > {output[1]}'
        )


# TODO: Should log exact reference genome used (in addition to species)
checkpoint reference_identification:
    input:
        'results/taxprofiler/bracken/summary.csv',
        'data/samplesheet.duckdb',
    output:
        'resources/reference_genomes.csv'
    params:
        f=config['identification']['minimum_fraction_reference'],
        p=config['identification']['minimum_readspow_reference'],
    log:
        'logs/smk/identification/reference_identification.log'
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
            ' {input[1]} -c ".read ' + 
            workflow.source_path('../scripts/match_reference_genome.sql') + 
            '" > {output}'
        )