# NOTE: Snakemake wildcards and Python scripting could simplify this rule
#       but at the cost of less portability as a module.

ruleorder: get_fastqs > get_resource

RESOURCES = ['samples', 'sequencing']


rule get_fastqs:
    output:
        'resources/library={library}/raw-fastqs.ext'
    params:
        glob=lambda wildcards: f"'{config['data'][wildcards.library]['fastqs']}*.fastq.gz'",
        pat=lambda wildcards: config['data'][wildcards.library]['pat'],
        model=config['models']['fastqs']
    log:
        'logs/smk/resources/{library}/get_fastqs.log'
    resources:
        njobs=20
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export GLOB={params.glob} PAT={params.pat};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' -c ".read {params.model}" > {output}'
        )


rule get_resource:
    output:
        'resources/library={library}/raw-{resource}.ext'
    params:
        path=lambda wildcards: config['data'][wildcards.library][wildcards.resource]
    log:
        'logs/smk/resources/{library}/get_{resource}.log'
    resources:
        njobs=20,
        slurm_partition='datatransfer'
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule raw_resources:
    input:
        expand(
            'resources/library={library}/raw-{resource}.ext',
            library=config['wildcards']['libraries'].split('|'),
            resource=['fastqs'] + RESOURCES
        )


rule clean_resource:
    input:
        ancient('resources/library={library}/raw-{resource}.ext')
    output:
        'resources/library={library}/{resource}.csv'
    params:
        model=lambda wildcards: config['models'][wildcards.resource][wildcards.library]
    log:
        'logs/smk/resources/{library}/clean_{resource}.log'
    resources:
        njobs=50
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export FN="{input}";' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' -c ".read {params.model}" > {output}'
        )


rule cleaned_resources:
    input:
        expand(
            'resources/library={library}/{resource}.csv',
            library=config['wildcards']['libraries'].split('|'),
            resource=['fastqs'] + ['samples', 'sequencing']
        )


# TODO: Fix full join on controls
# TODO: Figure out these: `select * from samplesheet where family = 'B001' and relationship is null;`
checkpoint samplesheet:
    input:
        expand(
            'resources/library={library}/{resource}.csv',
            library=config['wildcards']['libraries'].split('|'),
            resource=['fastqs'] + RESOURCES
        )
    output:
        'data/samplesheet.duckdb',
        'resources/samplesheets/main.csv',
    params:
        model=config['models']['samplesheet']
    log:
        'logs/smk/resources/samplesheet.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        runtime=5,
        njobs=1
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output[0]} -c ".read {params.model}" > {output[1]}'
        )