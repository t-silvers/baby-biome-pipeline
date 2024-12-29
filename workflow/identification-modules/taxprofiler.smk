rule taxprofiler_samplesheet:
    input:
        data_dir / 'data/samplesheet.csv',
    output:
        data_dir / 'data/samplesheets/taxprofiler.csv',
    log:
        log_dir / 'smk/identification/taxprofiler_samplesheet.log'
    resources:
        njobs=1,
    run:
        import pandas as pd

        TAXPROFILER_COLS = ['sample', 'run_accession', 'instrument_platform', 'fastq_1', 'fastq_2']

        (
            pd.read_csv(input[0])
            .assign(
                run_accession=lambda df: df['sample'],

                # NOTE: Must conform to EBI ENA controlled vocabulary
                instrument_platform='ILLUMINA'
            )
            .filter(TAXPROFILER_COLS)
            .dropna()
            .sort_values('sample')
            
            # TODO: Handle better
            .drop_duplicates(subset=['fastq_1'])
            
            .to_csv(output[0], index=False)
        )


rule taxprofiler:
    input:
        data_dir / 'data/samplesheets/taxprofiler.csv',
    output:
        data_dir / 'results/taxprofiler/multiqc/multiqc_report.html',
    params:
        # Dirs
        outdir=data_dir / 'results/taxprofiler/',
        workdir=log_dir / 'nxf/taxprofiler/work',
        
        # Generic params
        config=config['identification']['taxprofiler']['config'],
        profile=config['identification']['taxprofiler']['profiles'],
        
        # Pipeline params
        databases=config['identification']['taxprofiler']['databases'],
        extra=config['identification']['taxprofiler']['extra'],
    log:
        log_dir / 'smk/identification/taxprofiler.log'
    handover: True
    resources:
        njobs=295,
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run nf-core/taxprofiler \
          -config {params.config} \
          -profile {params.profile} \
          -resume \
          -work-dir {params.workdir} \
          --databases {params.databases} \
          --input {input} \
          --outdir {params.outdir} \
          {params.extra}
        '''


rule taxprofiler_bracken:
    """Collect kraken2-bracken output from taxprofiler."""
    input:
        ancient(data_dir / 'results/taxprofiler/multiqc/multiqc_report.html'),
    output:
        touch(data_dir / 'results/taxprofiler/bracken/bracken/{sample}_{sample}_bracken.bracken.tsv'),
    resources:
        njobs=1
    localrule: True


rule bracken_to_parquet:
    input:
        data_dir / 'results/taxprofiler/bracken/bracken/{sample}_{sample}_bracken.bracken.tsv',
    output:
        data_dir / 'data/identification/tool=taxprofiler/family={family}/id={id}/library={library}/{sample}.bracken.parquet',
    params:
        model=workflow.source_path(models['bracken']['staging']),
        pat='results/taxprofiler/bracken/\\S*/(\\d+)_\\d+_\\S*.bracken.tsv$',
    resources:
        cpus_per_task=2,
        runtime=5,
        njobs=1,
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)