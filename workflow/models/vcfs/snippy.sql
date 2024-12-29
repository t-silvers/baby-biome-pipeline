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
            , info_DP
            , info_RO
            , info_AO
            , info_TYPE
            , columns('format_\d+_GL') as format_GL
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
            , cast(info_RO as usmallint) as info_RO
            , array_transform(info_AO, x -> cast(x as usmallint)) as info_AO
            , info_TYPE
            , array_transform(format_GL, x -> cast(x as decimal(4, 1))) as format_GL
        from select_fields
    )

select * from optim_dtypes;


-- TODO: `snps.subs` output would be more accurate for (SNV) density calculation.
--       Could instead sum over allele lengths rather than occurrences.
create temp table alt_density as
with
    alt_density_data as (
        select "index"
            , chromosome
            , position
            , cast(alternate[1] is not null and (info_TYPE[1] != 'del' or info_TYPE[1] != 'ins') as int) as info_snv
            , cast((info_TYPE[1] != 'del' or info_TYPE[1] != 'ins') as int) as info_indel
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