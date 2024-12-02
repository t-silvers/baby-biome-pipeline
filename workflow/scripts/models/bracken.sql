create table bracken as

with raw_bracken_results as (

    select
        *
    
    from
        read_csv(
            getenv('GLOB'),
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

    select
        regexp_extract(
            "filename",
            'bracken/(\d+).bracken',
            1
        ) as "sample"
        , * exclude("filename")
    
    from
        raw_bracken_results

)

select * from cleaned;