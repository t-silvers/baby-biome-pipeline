set preserve_insertion_order = false;

create table annotated_vcfs as

with vcfs as (
    select
        *
    
    from
        read_csv(
            getenv('FILTERED_VCFS'),
            auto_detect = False,
            columns = {
                'CHROM': 'varchar',
                'POS': 'uinteger',
                'ID': 'varchar',
                'REF': 'varchar',
                'ALT': 'varchar',
                'QUAL': 'decimal(4, 1)',
                'FILTER': 'varchar',
                'INFO': 'varchar',
                'FORMAT': 'varchar',
                'FORMAT_values': 'varchar',
            },
            compression = 'gzip',
            delim = '\t',
            filename = True,
            header = True,
            ignore_errors = True,
            new_line = '\n', 
            nullstr = '.'
        )

),

final as (

    select
        regexp_extract("filename", 'variants/(\d+).filtered', 1) as "sample"
        , CHROM as chromosome
        , POS as position
        , REF as reference
        , ALT as alternate
        , QUAL as quality
        , array_transform(string_split(regexp_extract(INFO, 'ADF=([0-9,]+);', 1), ','), x -> cast(x as usmallint)) as allelic_read_depth_forward
        , array_transform(string_split(regexp_extract(INFO, 'ADR=([0-9,]+);', 1), ','), x -> cast(x as usmallint)) as allelic_read_depth_reverse
        , array_transform(string_split(regexp_extract(INFO, 'AD=([0-9,]+);', 1), ','), x -> cast(x as usmallint)) as allelic_read_depth
        , case when INFO ilike '%INDEL%' then true else false end as indel

    from
        vcfs

)

select * from final;