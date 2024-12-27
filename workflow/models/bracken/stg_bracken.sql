copy (
    with
        raw_bracken_results as (
            select * from read_csv(
                '{{ input }}',
                auto_detect = false,
                header = true,
                delim = '\t',
                columns = {
                    'name': 'varchar',
                    'taxonomy_id': 'varchar',
                    'taxonomy_lvl': 'varchar',
                    'kraken_assigned_reads': 'ubigint',
                    'added_reads': 'ubigint',
                    'new_est_reads': 'ubigint',
                    'fraction_total_reads': 'float4'
                },
                filename = true
            )
        ),

        cleaned as (
            select * exclude("filename")
                , regexp_extract("filename", '{{ pat }}', 1) as "sample"
            from raw_bracken_results
        )

    select * from cleaned
) to '{{ output }}' (format parquet);