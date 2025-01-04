include: '../identification-modules/taxprofiler.smk'
# include: '../identification-modules/mlst.smk'
include: '../identification-modules/srst2.smk'

BRACKEN_TEMPLATE = 'identification/tool=taxprofiler/family={family}/id={id}/library={library}/{sample}.bracken.parquet'

IDVARS = ['family', 'id', 'library', 'sample']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: Adapt this logic to non-multi-library data (and non-module?) use case

rule identification_db:
    output:
        data / 'identification.duckdb',
    params:
        model=workflow.source_path(models['workflow']['identification']['create']),
        idtools_seed=workflow.source_path(seeds['types']['idtools']),
    log:
        logdir / 'smk/identification/identification_db.log'
    resources:
        njobs=1,
    run:
        transform(params['model'], params, db=output[0], log=log[0])


rule update_identification_db:
    input:
        ancient(data / 'identification.duckdb'),
    output:
        results / 'samplesheets/identification_progress.csv',
    params:
        model=workflow.source_path(models['workflow']['identification']['insert']),
    log:
        logdir / 'smk/identification/update_identification_db.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        runtime=5,
        njobs=1
    run:        
        params.update({'output': output[0]})
        transform(params['model'], params, db=input[0], log=log[0])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def aggregate_bracken(wildcards):
    import pandas as pd

    def bracken_pq(df):
        return data_path_from_template(BRACKEN_TEMPLATE, df.to_dict())

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
        results / 'samplesheets/samplesheet.csv',
        aggregate_bracken
    output:
        results / 'samplesheets/reference_genomes.csv',
    params:        
        model=workflow.source_path(models['bracken']['reference_genome']),

        # Params
        bracken_glob=data / BRACKEN_TEMPLATE.format(**{v: "*" for v in IDVARS}),
        read_frac=config['identification']['minimum_fraction_reference'],
        read_pow=config['identification']['minimum_readspow_reference'],
    log:
        logdir / 'smk/identification/reference_identification.log'
    resources:
        cpus_per_task=8,
        mem_mb=4_000,
        runtime=15,
        njobs=1
    run:
        import pandas as pd

        params.update({'samplesheet': input[0], 'output': output[0]})
        transform(params['model'], params, log=log[0])

        if pd.read_csv(output[0]).empty:
            raise ValueError('No reference genomes found.')


# PHONY
rule all_identification:
    input:
        aggregate_srst2