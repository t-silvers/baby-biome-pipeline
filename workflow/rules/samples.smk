rule fastqs:
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


rule clean_samples:
    input:
        'results/raw_sample_info.xlsx'
    output:
        temp('results/sample_info.csv')
    params:
        model='/path/to/model.sql'
    envmodules:
        'duckdb/1.0'
    shell:
        (
            'export FN="{input}";' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-local') +
            ' -c ".read {params.model}" > {output}'
        )


rule clean_seq:
    input:
        'results/raw_seq_info.csv'
    output:
        temp('results/seq_info.csv')
    params:
        model='/path/to/model'
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
        'results/sample_info.csv',
        'results/seq_info.csv',
        'results/fastqs.csv',
    output:
        multiext('results/samplesheet', '.duckdb', '.csv')
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