# NOTE: Snakemake wildcards and Python scripting could simplify this rule
#       but at the cost of less portability as a module.

RESOURCES = ['samples', 'sequencing']

ruleorder: get_fastqs > get_resource


rule get_fastqs:
    output:
        'resources/library={library}/raw-fastqs.ext'
    params:
        glob=lambda wildcards: f"'{config[wildcards.library]['fastqs']}*.fastq.gz'",
        pat=lambda wildcards: config[wildcards.library]['pat'],
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
        path=lambda wildcards: config[wildcards.library][wildcards.resource]
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
    output:
        temp(touch('logs/smk/resources/raw_resources.done'))


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
    output:
        temp('logs/cleaned_resources.done')


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
        njobs=50
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {output[0]} -c ".read {params.model}" > {output[1]}'
        )