include: '../identification-modules/taxprofiler.smk'


def aggregate_bracken(wildcards):
    import pandas as pd

    def bracken_pq(df):
        path = 'data/identification/tool=taxprofiler/family={family}/id={id}/library={library}/{sample}.bracken.parquet'.format(**df.to_dict())
        return (data_dir / path).as_posix()

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[0]
    )

    return (
        samplesheet
        .transpose()
        .apply(lambda df: bracken_pq(df))
        .values
        .flatten()
    )


# TODO: Should log exact reference genome used (in addition to species)
checkpoint reference_identification:
    input:
        data_dir / 'data/samplesheet.csv',
        aggregate_bracken
    output:
        data_dir / 'data/identification/reference_genomes.csv',
    params:
        # TODO: 'bracken_glob' must be dynamic
        bracken_glob=data_dir / 'data/identification/*/*/*/*/*.bracken.parquet',
        
        model=workflow.source_path(models['bracken']['reference_genome']),
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


include: '../identification-modules/srst2.smk'


# PHONY
rule all_identification:
    input:
        aggregate_srst2