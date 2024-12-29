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
        data_dir / 'resources/library={library}/fastqs-raw.ext',
    params:
        glob=lambda wildcards: config['data'][wildcards.library]['fastqs'] + '*.fastq.gz',
        model=workflow.source_path(models['fastqs']['staging']),
    log:
        log_dir / 'smk/resources/{library}/get_fastqs.log'
    resources:
        njobs=1,
    run:
        params.update({'output': output[0]})
        transform(params['model'], params)


rule get_resource:
    output:
        data_dir / 'resources/library={library}/{resource}-raw.ext',
    params:
        path=resource_loc,
    log:
        log_dir / 'smk/resources/{library}/get_{resource}.log'
    resources:
        njobs=1,
        slurm_partition='datatransfer',
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule prepare_resource:
    input:
        data_dir / 'resources/library={library}/{resource}-raw.ext',
    output:
        data_dir / 'resources/library={library}/{resource}-pp.csv',
    params:
        model=resource_tform_model,
        pat=lambda wildcards: config['data'][wildcards.library]['pat'],
    log:
        log_dir / 'smk/resources/{library}/prepare_{resource}.log'
    resources:
        njobs=1,
    run:
        params.update( {'input': input[0], 'output': output[0]})
        transform(params['model'], params)


rule finalize_resource:
    input:
        data_dir / 'resources/library={library}/{resource}-pp.csv',
    output:
        data_dir / 'resources/library={library}/{resource}.csv',
    log:
        log_dir / 'smk/resources/{library}/finalize_{resource}.log'
    resources:
        njobs=1,
    shell:
        'cp "{input}" {output}'


checkpoint samplesheet:
    input:
        expand(
            data_dir / 'resources/library={library}/{resource}.csv',
            library=wc_params['libraries'],
            resource=RESOURCES
        )
    output:
        data_dir / 'data/samplesheet.csv',
    params:
        fastqs_glob=data_dir / 'resources/library=*/fastqs.csv',
        model=workflow.source_path(models['samplesheet']),
        samples_glob=data_dir / 'resources/library=*/samples.csv',
        sequencing_glob=data_dir / 'resources/library=*/sequencing.csv',
    log:
        log_dir / 'smk/resources/samplesheet.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        runtime=5,
        njobs=1
    run:
        params.update({'output': output[0]})
        transform(params['model'], params, log=log[0])