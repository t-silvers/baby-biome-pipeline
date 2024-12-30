ruleorder: get_fastqs > get_resource

RESOURCES = ['fastqs', 'samples', 'sequencing']

wildcard_constraints:
    resource='|'.join(RESOURCES)


def resource_loc(wildcards):
    return config['data'][wildcards.library][wildcards.resource]


def resource_tform_model(wildcards):
    libs = list(models[wildcards.resource].keys())
    if wildcards.library in libs:
        path = models[wildcards.resource][wildcards.library]
    else:
        path = models[wildcards.resource]["preparing"]
    return workflow.source_path(path)


rule get_fastqs:
    output:
        resources / 'library={library}/fastqs-raw.ext',
    params:
        glob=lambda wildcards: config['data'][wildcards.library]['fastqs'] + '*.fastq.gz',
        model=workflow.source_path(models['fastqs']['staging']),
    log:
        logdir / 'smk/resources/{library}/get_fastqs.log'
    resources:
        njobs=1,
    run:
        params.update({'output': output[0]})
        transform(params['model'], params)


rule get_resource:
    output:
        resources / 'library={library}/{resource}-raw.ext',
    params:
        path=resource_loc,
    log:
        logdir / 'smk/resources/{library}/get_{resource}.log'
    resources:
        njobs=1,
        slurm_partition='datatransfer',
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
        model=resource_tform_model,
        pat=lambda wildcards: config['data'][wildcards.library]['pat'],
    log:
        logdir / 'smk/resources/{library}/prepare_{resource}.log'
    resources:
        njobs=1,
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)


rule finalize_resource:
    input:
        resources / 'library={library}/{resource}-pp.csv',
    output:
        resources / 'library={library}/{resource}.csv',
    log:
        logdir / 'smk/resources/{library}/finalize_{resource}.log'
    resources:
        njobs=1,
    shell:
        'cp "{input}" {output}'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: Adapt this logic to non-multi-library data, non-module use case

rule samplesheet_db:
    output:
        data / 'samplesheet.duckdb',
    params:
        model=workflow.source_path(models['samplesheet']['create']),
        relationship_seed=workflow.source_path(seeds['types']['relationship']),
        species_seed=workflow.source_path(seeds['types']['species']),
        timepoint_seed=workflow.source_path(seeds['types']['timepoint']),
    log:
        logdir / 'smk/resources/samplesheet_db.log'
    resources:
        njobs=1,
    run:
        transform(params['model'], params, db=output[0], log=log[0])


checkpoint samplesheet:
    input:
        ancient(data / 'samplesheet.duckdb'),
        expand(
            resources / 'library={library}/{resource}.csv',
            library=wc_params['libraries'],
            resource=RESOURCES
        )
    output:
        results / 'samplesheets/samplesheet.csv',
    params:
        model=workflow.source_path(models['samplesheet']['insert']),
        fastqs_glob=resources / 'library=*/fastqs.csv',
        samples_glob=resources / 'library=*/samples.csv',
        sequencing_glob=resources / 'library=*/sequencing.csv',
    log:
        logdir / 'smk/resources/samplesheet.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        runtime=5,
        njobs=1
    run:
        params.update({'output': output[0]})
        transform(params['model'], params, log=log[0])

        import pandas as pd

        samplesheet = pd.read_csv(output[0])
        samplesheet = samplesheet[
            samplesheet['family'].isin(wc_params['families']) &
            samplesheet['library'].isin(wc_params['libraries']) &
            samplesheet['species_stdized'].isin(wc_params['species_plate'])
        ]

        if samplesheet.empty:
            raise ValueError('No samplesheet rows matched wildcards.')

        samplesheet.to_csv(output[0], index=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~