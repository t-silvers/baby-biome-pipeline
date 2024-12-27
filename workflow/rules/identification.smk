rule taxprofiler_samplesheet:
    input:
        data_dir / 'data/samplesheets/main.csv',
    output:
        data_dir / 'data/samplesheets/taxprofiler.csv',
    log:
        log_dir / 'smk/identification/taxprofiler_samplesheet.log'
    resources:
        njobs=1
    run:
        import pandas as pd

        TAXPROFILER_COLS = ['sample', 'run_accession', 'instrument_platform', 'fastq_1', 'fastq_2']

        (
            pd.read_csv(input[0])

            # ----------
            # TODO: TEMP
            .head(2)
            # ----------

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
        njobs=295
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


rule bracken_to_parquet:
    input:
        data_dir / 'results/taxprofiler/bracken/bracken/{sample}_{sample}_bracken.bracken.tsv',
    output:
        data_dir / 'data/identification/family={family}/id={id}/library={library}/{sample}.bracken.parquet',
    params:
        model=workflow.source_path(models['bracken']['staging']),
        pat='results/taxprofiler/bracken/\\S*/(\\d+)_\\d+_\\S*.bracken.tsv$',
    resources:
        cpus_per_task=2,
        runtime=5,
        njobs=1
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)


def aggregate_bracken(wildcards):
    import pandas as pd

    def covars_to_partitions(df, covars=None):
        """Form hive-partitioned path substring from df columns."""
        EXCLUDE = ['species', 'sample']

        if 'partitions' in df.columns:
            raise ValueError

        if covars is None:
            covars = [col for col in df.columns if col not in EXCLUDE]
        
        def row_to_partitions(row):
            return '/'.join([f'{k}={v}' for k, v in row.items()])

        return df.filter(covars).apply(row_to_partitions, axis=1)

    def bracken_pq(df):
        path = 'data/identification/family={family}/id={id}/library={library}/{sample}.bracken.parquet'.format(**df.to_dict())
        return (data_dir / path).as_posix()

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[0]
    )

    return (
        samplesheet

        # ----------
        # TODO: TEMP
        .head(2)
        # ----------

        .filter(['sample', 'family', 'id', 'library'])
        .transpose()
        .apply(lambda df: bracken_pq(df))
        .values
        .flatten()
    )


# TODO: Should log exact reference genome used (in addition to species)
checkpoint reference_identification:
    input:
        data_dir / 'data/samplesheets/main.csv',
        aggregate_bracken
    output:
        data_dir / 'data/identification/reference_genomes.csv',
    params:
        bracken_glob=data_dir / 'data/identification/*/*/*/*.bracken.parquet',
        model=workflow.source_path(models['reference_genome']),
        read_frac=config['identification']['minimum_fraction_reference'],
        read_pow=config['identification']['minimum_readspow_reference'],
    log:
        log_dir / 'smk/identification/reference_identification.log'
    resources:
        cpus_per_task=8,
        mem_mb=4_000,
        runtime=15,
        njobs=1
    run:
        params.update({'samplesheet': input[0], 'output': output[0]})
        transform(params['model'], params, log=log[0])