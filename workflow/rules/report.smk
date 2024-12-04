rule allele_depth_plot_data:
    input:
        'results/variants/{species}.duckdb',
        'results/pseudogenomes/species={species}/family={family}/filtered_calls.parquet',
    output:
        'reports/variants/species={species}/family={family}/allele_depth_plot_data.csv',
    resources:
        cpus_per_task=8,
        mem_mb=48_000,
        runtime=10,
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export FILTERED_CALLS={input[1]};' +
            'duckdb -readonly -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {input[0]} -c ".read ' + 
            workflow.source_path('../scripts/data_for_dp_plots.sql') + 
            '" > {output}'
        )
