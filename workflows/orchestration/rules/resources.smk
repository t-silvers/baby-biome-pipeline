# TODO: Find donor metadata for B002_B001_Lib pool, ie, in samplesheet:
#       `select * from samplesheet where relationship is null and family ilike 'B%';`

localrules:
    create_identification_db,
    create_samplesheet_db,
    create_variants_db,
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


rule create_identification_db:
    output:
        data / '.identification_duckdb',
    params:
        db=data / 'identification.duckdb',
    log:
        'logs/smk/create_identification_db.log'
    run:
        transform(models['create_identification_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


rule create_samplesheet_db:
    output:
        data / '.samplesheet_duckdb',
    params:
        db=data / 'samplesheet.duckdb',
        relationship_seed=seeds['relationship_type'],
        species_seed=seeds['species_type'],
        timepoint_seed=seeds['timepoint_type'],
    log:
        'logs/smk/create_samplesheet_db.log'
    run:
        transform(models['create_samplesheet_db'], params, db=params['db'], log=log[0])
        shell('touch ' + output[0])


def variants_db_path(wildcards):
    return data / 'variants' / f'species={wildcards.species}' / f'family={wildcards.family}' / f'tool={wildcards.mapping_tool}' / 'all.duckdb'


rule create_variants_db:
    output:
        data / '.variants_{species}_{family}_{mapping_tool}_duckdb',
    params:
        db=lambda wildcards: variants_db_path(wildcards),
        contigs_seed=seeds['contigs_type'],
        species=lambda wildcards: wildcards.species,
    log:
        'logs/smk/create_variants_db_{species}_{family}_{mapping_tool}.log'
    run:
        transform(models['create_variants_db'], params, db=params['db'], log=log[0])
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
        transform(models['stg_fastqs'], params, log=log[0])


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
        transform(models[f'clean_{wildcards.resource}_{wildcards.library}'], params, log=log[0])


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
        transform(workflow.source_path('../scripts/fix_index_swap_230119_B001_Lib.sql'), params)


# NOTE: Global -> Local resources


rule prepare_samplesheet:
    input:
        expand(
            resources / 'library={library}/{resource}.csv',
            library=wc_params['library'],
            resource=wc_params['resource'],
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
        transform(models['prepare_samplesheet'], params, log=log[0])


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
        transform(models['update_samplesheet_db'], params, db=params['db'], log=log[0])


checkpoint samplesheet:
    input:
        'resources/samplesheet-pp.csv',
    output:
        'resources/samplesheet.csv',
    log:
        'logs/smk/samplesheet.log'
    run:
        import pandas as pd

        # TODO: Can either filter by requested wildcards
        #       OR
        #       Determine wildcards values

        samplesheet = pd.read_csv(input[0]).dropna(subset=['fastq_1', 'fastq_2'])

        # TODO: Replace `isin` with regex (as with `wildcard_constraints`).
        #       Deprecate `wc_params` in common.smk.
        samplesheet = samplesheet[
            samplesheet['family'].isin(wc_params['family']) &
            samplesheet['library'].isin(wc_params['library']) &
            samplesheet['species_stdized'].isin(wc_params['species_plate'])
        ]

        if samplesheet.empty:
            raise ValueError('No samplesheet rows matched wildcards.')

        # TODO: TEMP
        # samplesheet = samplesheet.sample(5)

        samplesheet.to_csv(output[0], index=False)