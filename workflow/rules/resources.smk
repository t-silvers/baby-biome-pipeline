# NOTE: Snakemake wildcards and Python scripting could simplify this rule
#       but at the cost of less portability?

RESOURCES = ['fastqs', 'samples', 'sequencing']


rule get_fastqs:
    output:
        'resources/library={library}/raw-fastqs.ext'
    params:
        glob=lambda wildcards: f"'{config[wildcards.library]['fastqs']}*.fastq.gz'",
        pat=lambda wildcards: config[wildcards.library]['pat']
    resources:
        njobs=20
        slurm_partition='datatransfer'
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export GLOB={params.glob} PAT={params.pat};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' -c ".read ' +
            workflow.source_path('../scripts/models/fastqs.sql') +
            '" > {output}'
        )


rule:
    name:
        'get_{resource}'
    output:
        'resources/library={library}/raw-{resource}.ext'
    params:
        path=lambda wildcards: config[wildcards.library][wildcards.resource]
    resources:
        njobs=20
        slurm_partition='datatransfer'
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule:
    name:
        'clean_{resource}'
    input:
        ancient('resources/library={library}/raw-{resource}.ext')
    output:
        'resources/library={library}/{resource}.csv'
    params:
        model=lambda wildcards: f'workflow/models/{wildcards.resource}-{wildcards.library}.sql'
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


checkpoint samplesheet:
    input:
        expand(
            'resources/library={library}/{resource}.csv',
            library=config['wildcards']['libraries'].split('|'),
            resource=RESOURCES
        )
    output:
        'data/samplesheet.duckdb',
        'resources/samplesheets/main.csv',
    resources:
        njobs=50
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {output[0]} -c ".read ' + 
            workflow.source_path('../scripts/create_samplesheet.sql') + 
            '" > {output[1]}'
        )