rule sample_info:
    output:
        'results/raw_sample_info.xlsx'
    params:
        path='/path/to/file'
    resources:
        slurm_partition='datatransfer'
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule seq_info:
    output:
        'results/raw_seq_info.csv'
    params:
        path='/path/to/file'
    resources:
        slurm_partition='datatransfer'
    envmodules:
        'rclone/1.67.0'
    shell:
        'rclone copyto "nextcloud:{params.path}" {output}'


rule log_samples:
    input:
        'results/raw_sample_info.xlsx'
    output:
        'data/sample_info.duckdb'
    params:
        model='/path/to/model.sql'
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export FN="{input}";' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {output} -c ".read {params.model}"'
        )


rule log_seq:
    input:
        'results/raw_seq_info.csv'
    output:
        'data/seq_info.duckdb'
    params:
        model='/path/to/model'
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export FN="{input}";' + 
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {output} -c ".read {params.model}"'
        )


rule glob_fastqs:
    output:
        'results/fastqs.csv'
    params:
        glob="'/path/to/*.fastq.gz'",
        pat=''
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


checkpoint fastqs:
    input:
        db='data/sample_info.duckdb',
        fastqs='results/fastqs.csv',
    output:
        'results/samplesheet.csv'
    localrule: True
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export FASTQS="{input.fastqs}";' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' {input.db} -c ".read ' + 
            workflow.source_path('../scripts/create_samplesheet.sql') + 
            '" > {output}'
        )