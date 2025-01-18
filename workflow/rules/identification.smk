taxprofiler_params = config['params']['taxprofiler']


checkpoint taxprofiler_samplesheet:
    input:
        sample_info='resources/samplesheet.csv',
    output:
        samplesheet='resources/samplesheets/taxprofiler.csv',
    params:
        # Files
        db=data / 'identification_cache.duckdb',
        input=input['sample_info'],
        output=output['samplesheet'],
        
        # Params
        # NOTE: Must conform to EBI ENA controlled vocabulary
        instrument_platform='ILLUMINA',
    log:
        'logs/smk/taxprofiler_samplesheet.log'
    run:
        transform(
            models['nfcore_inputs']['taxprofiler'],
            params, db=params['db'], log=log[0], readonly=True
        )


rule taxprofiler:
    input:
        'resources/samplesheets/taxprofiler.csv',
    output:
        'results/taxprofiler/multiqc/multiqc_report.html',
        'results/taxprofiler/bracken/bracken_key_k2db_combined_reports.txt',
    params:
        pipeline=taxprofiler_params['pipeline'],

        # Dirs
        outdir=subpath(output[0], ancestor=2),
        workdir='logs/nxf/taxprofiler_work',

        # Generic params
        config=taxprofiler_params['config'],
        profile=taxprofiler_params['profiles'],
        version=taxprofiler_params['version'],

        # Pipeline params
        databases=taxprofiler_params['databases'],
        extra=taxprofiler_params['extra'],
    log:
        'logs/smk/identification/taxprofiler.log'
    handover: True
    envmodules:
        'apptainer',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run {params.pipeline} \
          -config {params.config} \
          -profile {params.profile} \
          -resume \
          -work-dir {params.workdir} \
          --databases {params.databases} \
          --input {input} \
          --outdir {params.outdir} \
          {params.extra}
        '''


rule clean_bracken:
    input:
        'results/taxprofiler/multiqc/multiqc_report.html',
    output:
        data / 'identification/tool=taxprofiler/db={db_name}/family={family}/id={id}/library={library}/{sample}.bracken.parquet',
    params:
        input=lambda wildcards: f'results/taxprofiler/bracken/{wildcards.db_name}/{wildcards.sample}_{wildcards.sample}_{wildcards.db_name}.bracken.tsv',
        output=output[0],
        sample=lambda wildcards: wildcards.sample,
    run:
        transform(models['bracken']['staging'], params)


# TODO: Break dependency on a specific mapper
# TODO: Provide separate trimming rule for SRST2
#       (avoids re-running expensive nf-core
#       pipelines just to get trimmed reads,
#       without needing to store them).
# TODO: Use db to track progress
checkpoint srst2_samplesheet:
    input:
        rules.all_bactmap.input,
        sample_info='resources/samplesheet-with-ref.csv',
    output:
        samplesheet='resources/samplesheets/srst2_{species}.csv',
    params:
        glob=lambda wildcards: f'results/bactmap/{wildcards.species}/fastp/*.trim.fastq.gz',
        pat='fastp/(\\d+)_(1|2).trim.fastq.gz$',
        output=output['samplesheet']
    log:
        'logs/smk/srst2_samplesheet_{species}.log'
    run:
        transform(models['srst2']['samplesheet'], params, log=log[0])


def trimmed_fastqs(wildcards):
    import pandas as pd

    sample = int(wildcards.sample)
    
    return (
        pd.read_csv(
            checkpoints.srst2_samplesheet
            .get(species=wildcards.species)
            .output['samplesheet']
        )
        .query('sample == @sample')
        .filter(['trim_fastq_1', 'trim_fastq_2'])
        .values
        .flatten()
    )


rule srst2:
    """Profile sequence types using srst2.

    Notes:
      - `srst2` allows `--use_existing_bowtie2_sam` and `--use_existing_pileup`, but
          1. uses restrictive pattern matching
          2. the mapping reference differs from the ST reference, which is only
              the core genes; the mapping reference will be of an 'unknown'
              (findable) sequence type
      - input may be zipped or not; however, `srst2` will infer this based on
          extension and fail if the inferred type differs from actual type.
      - `--forward` and `--reverse` flags refer to the basename before extension
          (i.e., '_R1')
      - `--use_existing_pileup` is a boolean flag,
          `parser.add_argument('--use_existing_pileup', action="store_true" ...`
      - whatever is passed to `--output` will not be the output; rather, the
          output for `--output sample` will be `sample__mlst__<db>__results.txt`
    """
    input:
        trimmed_fastqs,
    output:
        'results/srst2/{species}/{sample}.txt',
    params:
        # Dirs
        prefix=lambda wildcards: f'results/srst2/{wildcards.species}/{wildcards.sample}',

        # srst2 params
        extra=config['params']['srst2']['extra'],
        mlst_db=lambda wildcards: config['data']['public']['mlst_schema'][wildcards.species]['db'],
        mlst_definitions=lambda wildcards: config['data']['public']['mlst_schema'][wildcards.species]['definitions'],
        species_alias=lambda wildcards: config['data']['public']['mlst_schema'][wildcards.species]['alias'],
    envmodules:
        'srst2/0.2.0'
    container:
        'docker://staphb/srst2'
    shell:
        """
        srst2 \
          --input_pe {input} \
          --output {params.prefix} \
          {params.extra} \
          --mlst_db '{params.mlst_db}' \
          --mlst_definitions '{params.mlst_definitions}' \
          --threads {resources.cpus_per_task}

        rm -f '{params.prefix}__{wildcards.sample}.{params.species_alias}.pileup'
        rm -f '{params.prefix}__{wildcards.sample}.{params.species_alias}.scores'
        rm -f '{params.prefix}__{wildcards.sample}.{params.species_alias}.sorted.bam'

        mv '{params.prefix}__mlst__{params.species_alias}__results.txt' {output}
        """


rule clean_srst2:
    input:
        'results/srst2/{species}/{sample}.txt',
    output:
        data / 'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}.txt',
    params:
        input=input[0],
        output=output[0],
    run:
        transform(models['srst2']['staging'], params)