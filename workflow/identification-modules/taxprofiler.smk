import pathlib


def get_bracken_db_info() -> dict:
    import pandas as pd

    bracken_dbs = (
        pd.read_csv(config['identification']['taxprofiler']['databases'])
        .query('tool == "bracken"')
    )

    if bracken_dbs.shape[0] != 1:
        raise ValueError('Expected exactly one bracken database.')

    return bracken_dbs.to_dict('records')[0]


bracken_db_info = get_bracken_db_info()


rule taxprofiler_samplesheet:
    input:
        results / 'samplesheets/samplesheet.csv',
    output:
        results / 'samplesheets/taxprofiler.csv',
    log:
        logdir / 'smk/identification/taxprofiler_samplesheet.log'
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

        # TODO: Remove samples that have already been analyzed; otherwise,
        #       relies on nxf dependency tracking. Use run db (sample,tool,...).


rule taxprofiler:
    input:
        results / 'samplesheets/taxprofiler.csv',
    output:
        results / 'taxprofiler/multiqc/multiqc_report.html',
        results / 'taxprofiler/bracken/bracken_key_k2db_combined_reports.txt',
    params:
        # Dirs
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent.parent,
        workdir=logdir / 'nxf/taxprofiler_work',

        # Generic params
        config=config['identification']['taxprofiler']['config'],
        profile=config['identification']['taxprofiler']['profiles'],
        version=config['identification']['taxprofiler']['version'],

        # Pipeline params
        databases=config['identification']['taxprofiler']['databases'],
        extra=config['identification']['taxprofiler']['extra'],
    log:
        logdir / 'smk/identification/taxprofiler.log'
    handover: True
    resources:
        njobs=max_submit,
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run nf-core/taxprofiler \
          -config {params.config} \
          -profile {params.profile} \
          -r {params.version} \
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
        ancient(results / 'taxprofiler/multiqc/multiqc_report.html'),
    output:
        touch(results / 'taxprofiler/bracken/{db_name}/{sample}_{sample}_{db_name}.bracken.tsv'),
    resources:
        njobs=1
    localrule: True


# TODO: Could relax requirement for 1 db and collect later (while enforcing that reference
#       genome identification is based on 1 specific db).
def aggregate_bracken_dbs(wildcards):
    """Taxprofiler bracken will be nested by database name."""
    return results / 'taxprofiler/bracken/{db_name}/{{sample}}_{{sample}}_{db_name}.bracken.tsv'.format(**bracken_db_info)


rule bracken_to_parquet:
    input:
        aggregate_bracken_dbs,
    output:
        data / 'identification/tool=taxprofiler/family={family}/id={id}/library={library}/{sample}.bracken.parquet',
    params:
        model=workflow.source_path(models['bracken']['staging']),
        pat=results / 'taxprofiler/bracken/{db_name}/(\\d+)_\\d+_{db_name}.bracken.tsv$'.format(**bracken_db_info),
    resources:
        cpus_per_task=2,
        runtime=5,
        njobs=1,
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)