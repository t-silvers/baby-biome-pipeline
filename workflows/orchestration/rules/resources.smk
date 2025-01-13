# TODO: Find donor metadata for B002_B001_Lib pool, ie, in samplesheet:
#       `select * from samplesheet where relationship is null and family ilike 'B%';`

localrules:
    finalize_resource,
    fix_index_swap,
    get_fastqs,
    get_resource,
    prepare_resource,
    prepare_samplesheet,
    samplesheet,
    update_samplesheet_db,

ruleorder:
    get_fastqs > get_resource

ruleorder:
    fix_index_swap > finalize_resource


rule create_db:
    output:
        data / ('.{db_name}_duckdb'),
    params:
        db=lambda wildcards: data / (wildcards.db_name + '.duckdb'),
        idtools=seeds['types']['idtools'],
        genomes=seeds['types']['reference_genomes'],
        maptools=seeds['types']['maptools'],
        relationship=seeds['types']['relationship'],
        species=seeds['types']['species'],
        timepoint=seeds['types']['timepoint'],
    log:
        'logs/smk/create_{db_name}_db.log',
    localrule: True
    run:
        transform(models[wildcards.db_name]['create_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


def progress_cache_glob(wildcards):
    if wildcards.db_name == 'identification_cache':
        return data / 'identification/tool=taxprofiler/**/*.bracken.parquet'
    elif wildcards.db_name == 'mapping_cache':
        return data / 'variants/tool=*/**/*.cleaned.parquet'


def progress_cache_pat(wildcards):
    if wildcards.db_name == 'identification_cache':
        return 'identification/tool=(taxprofiler)/.*/(\\d+).bracken.parquet$'
    elif wildcards.db_name == 'mapping_cache':
        return f'variants/tool=({config["tools"]["mapping"]})/.*/(\\d+).cleaned.parquet$'


rule update_progress_cache:
    input:
        data / ('.{db_name}_duckdb'),
    output:
        temp('results/.{db_name}'),
    params:
        db=lambda wildcards: data / (wildcards.db_name + '.duckdb'),
        glob=lambda wildcards: progress_cache_glob(wildcards),
        pat=lambda wildcards: progress_cache_pat(wildcards),
    log:
        'logs/smk/update_{db_name}.log'
    localrule: True
    run:
        transform(models[wildcards.db_name]['insert_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


rule get_fastqs:
    output:
        resources / 'library={library}/fastqs-raw.ext',
    params:
        glob=lambda wildcards: config['data'][wildcards.library]['fastqs'] + '*.fastq.gz',
        output=output[0],
    log:
        'logs/smk/get_{library}_fastqs.log'
    run:
        transform(models['fastqs']['staging'], params, log=log[0])


rule get_resource:
    output:
        resources / 'library={library}/{resource}-raw.ext',
    params:
        path=lambda wildcards: config['data'][wildcards.library][wildcards.resource],
    log:
        'logs/smk/get_{library}_{resource}.log'
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule prepare_resource:
    input:
        resources / 'library={library}/{resource}-raw.ext',
    output:
        resources / 'library={library}/{resource}-pp.csv',
    params:
        pat=lambda wildcards: config['data'][wildcards.library]['pat'],
        input=input[0],
        output=output[0],
    log:
        'logs/smk/prepare_{library}_{resource}.log'
    run:
        transform(
            models[wildcards.resource]['clean'][wildcards.library],
            params, log=log[0]
        )


rule finalize_resource:
    input:
        resources / 'library={library}/{resource}-pp.csv',
    output:
        resources / 'library={library}/{resource}.csv',
    shell:
        'cp "{input}" {output}'


rule fix_index_swap:
    input:
        resources / 'library=230119_B001_Lib/fastqs-pp.csv',
        resources / 'library=230119_B001_Lib/sequencing-pp.csv'
    output:
        resources / 'library=230119_B001_Lib/fastqs.csv',
        resources / 'library=230119_B001_Lib/sequencing.csv'
    params:
        fastqs=input[0],
        sequencing=input[1],
        f_output=output[0],
        s_output=output[1],
    log:
        'logs/smk/resources/230119_B001_Lib/fix_index_swap.log'
    run:
        transform(
            workflow.source_path('../scripts/fix_index_swap_230119_B001_Lib.sql'),
            params, log=log[0]
        )


rule prepare_samplesheet:
    input:
        expand(
            resources / 'library={library}/{resource}.csv',
            library=config['wildcards']['library'].split('|'),
            resource=config['wildcards']['resource'].split('|'),
        )
    output:
        'resources/samplesheet-raw.csv',
    params:
        fastqs_glob=resources / 'library=*/fastqs.csv',
        samples_glob=resources / 'library=*/samples.csv',
        sequencing_glob=resources / 'library=*/sequencing.csv',
        output=output[0],
    log:
        'logs/smk/prepare_samplesheet.log'
    run:
        transform(models['samplesheet']['clean'], params, log=log[0])


rule update_samplesheet_db:
    input:
        'resources/samplesheet-raw.csv',
        data / '.samplesheet_duckdb',
    output:
        'resources/samplesheet-pp.csv',
    params:
        db=data / 'samplesheet.duckdb',
        input=input[0],
        output=output[0],
    log:
        'logs/smk/update_samplesheet_db.log'
    run:
        transform(models['samplesheet']['update_db'], params, db=params['db'], log=log[0])


rule samplesheet:
    input:
        'resources/samplesheet-pp.csv',
    output:
        'resources/samplesheet.csv',
    log:
        'logs/smk/samplesheet.log'
    run:
        import pandas as pd

        samplesheet = pd.read_csv(input[0]).dropna(subset=['fastq_1', 'fastq_2'])
        
        for var in config['wildcards']:
            if var in samplesheet.columns:
                # TODO: Replace `isin` with regex (as with `wildcard_constraints`).
                #       Deprecate `wc_params` in common.smk.
                var_mask = samplesheet[var].isin(config['wildcards'][var].split('|'))
                samplesheet = samplesheet[var_mask]

            if samplesheet.empty:
                raise ValueError('No samplesheet rows matched wildcards.')

        samplesheet.to_csv(output[0], index=False)