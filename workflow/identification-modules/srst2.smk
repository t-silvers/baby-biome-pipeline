rule bactmap_trimmed_fastqs:
    """Collect fastp output from bactmap."""
    input:
        results / 'bactmap/{species}/pipeline_info/pipeline_report.html',
    output:
        multiext(
            (results / 'bactmap/{species}/fastp/{sample}').as_posix(),
            '_1.trim.fastq.gz', '_2.trim.fastq.gz',
        ),
    resources:
        njobs=1,
    localrule: True


rule srst2:
    """Profile sequence types using srst2.

    Notes
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
        multiext(
            (results / 'bactmap/{species}/fastp/{sample}').as_posix(),
            '_1.trim.fastq.gz', '_2.trim.fastq.gz',
        ),
    output:
        results / 'srst2/{species}/{sample}_results.txt',
    params:
        # Dirs
        prefix=lambda wildcards: results / 'srst2' / wildcards.species / wildcards.sample,
        
        # srst2 params
        extra=config['identification']['srst2']['extra'],
        mlst_db=lambda wildcards: config['public_data']['mlst'][wildcards.species]['db'],
        mlst_definitions=lambda wildcards: config['public_data']['mlst'][wildcards.species]['definitions'],
        species_alias=lambda wildcards: config['public_data']['mlst'][wildcards.species]['alias'],
    resources:
        cpus_per_task=16,
        mem_mb=4_000,
        runtime=15,
        njobs=1,
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
        results / 'srst2/{species}/{sample}_results.txt',
    output:
        data / 'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}_results.txt',
    params:
        model=workflow.source_path(models['srst2']),
    resources:
        cpus_per_task=2,
        runtime=5,
        njobs=1,
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)


def aggregate_srst2(wildcards):
    import pandas as pd

    # TODO: Could replace with `mlst --list` check
    available_schemas = list(config['public_data']['mlst'].keys())

    def srst2_results(df):
        SRST2_TEMPLATE = 'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}_results.txt'

        return data_path_from_template(SRST2_TEMPLATE, df.to_dict())

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .rename(columns={'reference_genome': 'species'})
        .query('species in @available_schemas')
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: srst2_results(df))
        .values
        .flatten()
    )


# PHONY
rule all_srst2:
    input:
        aggregate_srst2
