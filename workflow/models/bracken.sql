create table bracken as

with raw_bracken_results as (

    select
        *
    
    from
        read_csv(
            getenv('TAXPROFILER_BRACKEN'),
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
            'results/taxprofiler/bracken/\S*/(\d+)_\d+_\S*.bracken.tsv$',
            1
        ) as "sample"
        , * exclude("filename")
    
    from
        raw_bracken_results

)

select * from cleaned;

copy (
    
    select * from bracken

) to '/dev/stdout' (format csv);