create temp table vcf_clean as
with
    raw_vcf_pq as (
        select * from read_parquet('{{ input }}', hive_partitioning = false)
    ),

    select_fields as (
        select row_number() over () as "index"
            , chromosome
            , position
            , reference
            , alternate
            , quality
            , "filter"
            , info_INDEL
            , info_DP
            , info_AD
            , info_ADF
            , info_ADR
            , columns('format_\d+_SP') as format_SP
        from raw_vcf_pq
    ),

    optim_dtypes as (
        select "index"
            , chromosome
            , position
            , reference
            , array_transform(alternate, x -> nullif(x, '')) as alternate
            , cast(quality as decimal(4, 1)) as quality
            , array_transform("filter", x -> x = 'PASS') as "filter"
            , info_INDEL
            , cast(info_DP as usmallint) as info_DP
            , array_transform(info_AD, x -> cast(x as usmallint)) as info_AD
            , array_transform(info_ADF, x -> cast(x as usmallint)) as info_ADF
            , array_transform(info_ADR, x -> cast(x as usmallint)) as info_ADR
            , cast(format_SP as usmallint) as format_SP
        from select_fields
    )

select * from optim_dtypes;


create temp table alt_density as
with
    alt_density_data as (
        select "index"
            , chromosome
            , position
            , cast(alternate[1] is not null and not info_INDEL as int) as info_snv
            , cast(info_INDEL as int) as info_indel
        from vcf_clean
    ),

    density_calc as (
        select "index"
            , cast(sum(info_snv) over adw as usmallint) as SNV_density
            , cast(sum(info_indel) over adw as usmallint) as INDEL_density
        from alt_density_data
        window adw as (
            partition by chromosome
            order by position
            range between cast('{{ alt_density_window }}' as int) preceding and 
                          cast('{{ alt_density_window }}' as int) following
        )
    )

select * from density_calc;


copy (
    select v.* exclude("index")
        , d.SNV_density
        , d.INDEL_density
    from vcf_clean v
    join alt_density d on (v.index = d.index)
) to '{{ output }}' (format parquet);