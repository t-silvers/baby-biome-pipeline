rule:
    input:
        'results/variants/species={species}/family={family}/variants.duckdb',
    output:
        'results/variants/species={species}/family={family}/pseudogenomes.csv',
    params:
        dp=config['pseudogenome']['total_read_depth'],
        maf=config['pseudogenome']['maf'],
        qual=config['pseudogenome']['quality'],
        sdp=config['pseudogenome']['strand_read_depth'],
    resources:
        cpus_per_task=12,
        mem_mb=16_000,
        runtime=5,
    envmodules:
        'duckdb/nightly'
    # TODO: Must revise `filter_genotype_calls.sql`
    #       Ensure reference is still in table (for MSA)
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export DP={params.dp};' +
            'export MAF={params.maf};' +
            'export QUAL={params.qual};' +
            'export STRAND_DP={params.sdp};' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {input} -c ".read ' + 
            workflow.source_path('../scripts/filter_genotype_calls.sql') + 
            '" > {output}'
        )