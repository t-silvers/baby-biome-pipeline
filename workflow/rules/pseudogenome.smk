rule:
    input:
        'results/data/variants/{species}.duckdb',
    output:
        'results/{species}/pseudogenomes/annot_filtered_calls.csv',
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
    shell:
        '''
        export  MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
                DP={params.dp} \
                MAF={params.maf} \
                QUAL={params.qual} \
                STRAND_DP={params.sdp}

        duckdb -readonly -init config/duckdbrc-slurm {input} \
          -c ".read workflow/scripts/filter_genotype_calls.sql" > {output}
        '''
