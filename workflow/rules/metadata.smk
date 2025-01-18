"""Create and log metadata about pipeline.

If no databases (`data / '{{ DB }}.duckdb'`) provided, creates one. If no data
infrastructure library provided, uses default scripts for staging, cleaning,
and inserting into databases.

Usage example:

  snakemake log_metadata
"""
localrules:
    create_metadata_db,
    update_progress_cache,


# PHONY
rule metadata:
    input:
        ''


rule create_metadata_db:
    output:
        metadata_db=data / 'metadata.duckdb',
    params:
        enum_seeds=seeds['types'],
    log:
        'logs/smk/create_metadata_db.log',
    run:
        transform(
            models['pipeline_metadata']['create_db'],
            params['enum_seeds'],
            db=output['metadata_db'],
            log=log[0],
        )


rule update_progress_cache:
    input:
        data / ('.{db_name}_duckdb'),
    output:
        update(data / (wildcards.db_name + '.duckdb')),
    params:
        db=lambda wildcards: data / (wildcards.db_name + '.duckdb'),
        glob=lambda wildcards: progress_cache_glob(wildcards),
        pat=lambda wildcards: progress_cache_pat(wildcards),
    log:
        'logs/smk/update_{db_name}.log'
    run:
        transform(models[wildcards.db_name]['insert_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


rule update_progress_cache:
    input:
        data / ('.{db_name}_duckdb'),
    output:
        update(data / (wildcards.db_name + '.duckdb')),
    params:
        db=lambda wildcards: data / (wildcards.db_name + '.duckdb'),
        glob=lambda wildcards: progress_cache_glob(wildcards),
        pat=lambda wildcards: progress_cache_pat(wildcards),
    log:
        'logs/smk/update_{db_name}.log'
    run:
        transform(models[wildcards.db_name]['insert_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


def progress_cache_glob(wildcards):
    if wildcards.db_name == 'identification_cache':
        return data / 'identification/tool=taxprofiler/**/*.bracken.parquet'
    elif wildcards.db_name == 'mapping_cache':
        return data / 'variants/tool=*/**/*.parquet'


def progress_cache_pat(wildcards):
    if wildcards.db_name == 'identification_cache':
        return 'identification/tool=(taxprofiler)/.*/(\\d+).bracken.parquet$'
    elif wildcards.db_name == 'mapping_cache':
        return f'variants/tool=({config["tools"]["mapping"]})/.*/(\\d+).parquet$'