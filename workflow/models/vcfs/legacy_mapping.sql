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
            , info_RPBZ
            , info_MQBZ
            , info_BQBZ
            , info_SCBZ
            , info_DP4
            , columns('format_.*_PL') as format_PL
        from raw_vcf_pq
    ),

    optim_dtypes as (
        select "index"
            , chromosome
            , position
            , reference
            , alternate
            , cast(quality as decimal(5, 1)) as quality
            , array_transform("filter", x -> x = '') as "filter"
            , cast(info_DP as usmallint) as info_DP
            , cast(info_RPBZ as decimal(2, 1)) as info_RPBZ
            , cast(info_MQBZ as decimal(2, 1)) as info_MQBZ
            , cast(info_BQBZ as decimal(2, 1)) as info_BQBZ
            , cast(info_SCBZ as decimal(2, 1)) as info_SCBZ
            , array_transform(info_DP4, x -> cast(x as usmallint)) as info_DP4
            , array_transform(format_PL, x -> cast(x as usmallint)) as format_PL
        from select_fields
    )

select * from optim_dtypes;


-- TODO: Change `bcftools view` args to estimate INDEL density
create temp table alt_density as
with
    alt_density_data as (
        select "index"
            , chromosome
            , position
            , cast(alternate[1] is not null as int) as info_snv
        from vcf_clean
    ),

    density_calc as (
        select "index"
            , cast(sum(info_snv) over adw as usmallint) as SNV_density
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
    from vcf_clean v
    join alt_density d on (v.index = d.index)
) to '{{ output }}' (format parquet);