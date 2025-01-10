copy (
    with
        raw_vcf_pq as (
            select
                try_cast(
                    regexp_extract("filename", '/(\d+).parquet$', 1)
                    as uinteger
                ) as "sample"
                , chromosome
                , position
                , reference
                , alternate
                , quality
                , "filter"
                , info_AD
                , info_ADF
                , info_ADR
                , info_MQ
                , info_DP
                , info_INDEL
                , columns('format_.*_GT') as format_GT
                , columns('format_.*_PL') as format_PL
                , columns('format_.*_SP') as format_SP
            from read_parquet(
                '{{ input }}',
                filename = true,
                hive_partitioning = false
            )
        ),

        filtered as (
            select * exclude("filter")
            from raw_vcf_pq
            where quality >= cast('{{ quality }}' as float)
                and info_DP >= cast('{{ dp }}' as float)
                and (info_AD[1] / array_reduce(info_AD, (x, y) -> x + y)) < (1 - cast('{{ maf }}' as float))
                and array_reduce(info_AD, (x, y) -> x + y) >= cast('{{ ad }}' as float)
                and array_reduce(info_ADF, (x, y) -> x + y) >= cast('{{ ad_strand }}' as float)
                and array_reduce(info_ADR, (x, y) -> x + y) >= cast('{{ ad_strand }}' as float)
                and info_MQ >= cast('{{ mq }}' as float)
                and format_SP < cast('{{ sp }}' as float)
        ),

        cleaned as (
            select chromosome as contig
                , position as start_pos
                , "sample"
                , position + length(reference) - 1 as end_pos
                , cast(quality as decimal(4, 1)) as qual
                , array_concat([reference], alternate) as alleles

                -- info fields
                , array_transform(info_AD, x -> cast(x as usmallint)) as info_AD
                , array_transform(info_ADF, x -> cast(x as usmallint)) as info_ADF
                , array_transform(info_ADR, x -> cast(x as usmallint)) as info_ADR
                , cast(info_DP as usmallint) as info_DP
                , cast(info_MQ as decimal(4, 1)) as info_MQ
                , case
                    when info_INDEL is true then ['indel']
                    when info_INDEL is false then ['snp']
                    else [null]
                end
                as info_TYPE
                
                -- format fields
                , cast(format_GT as varchar) as format_GT
                , array_transform(format_PL, x -> cast(x as usmallint)) as format_PL
                , cast(format_SP as float) as format_SP

            from filtered
        ),

        final as (select * from cleaned)
    select * from final
) to '{{ output }}' (format parquet);