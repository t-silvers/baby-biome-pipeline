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
                , columns('format_.*_AD') as info_AD
                , info_DP
                , info_MQM
                , info_TYPE
                , columns('format_.*_GT') as format_GT
                , columns('format_.*_GL') as format_GL
            from read_parquet(
                '{{ input }}',
                filename = true,
                hive_partitioning = false
            )
        ),

        filtered as (
            select * exclude("filter")
            from raw_vcf_pq
            where quality >= 25
                and info_DP >= 10
                and (info_AD[1] / array_reduce(info_AD, (x, y) -> x + y)) < .05
                and array_reduce(info_AD, (x, y) -> x + y) >= 6
        ),

        cleaned as (
            select chromosome as contig
                , position as start_pos
                , "sample"
                , position + length(reference) - 1 as end_pos
                , cast(quality as decimal(5, 1)) as qual
                , array_concat([reference], alternate) as alleles

                -- info fields
                , array_transform(info_AD, x -> cast(x as usmallint)) as info_AD
                , cast(info_DP as usmallint) as info_DP
                , cast(info_MQM[1] as decimal(4, 1)) as info_MQ
                , array_transform(
                    info_TYPE,
                    x -> case
                        when x = 'complex' then 'indel'
                        when x = 'del' then 'indel'
                        when x = 'ins' then 'indel'
                        when x = 'mnp' then 'indel'
                        when x = 'snp' then x
                        else null
                    end
                ) as info_TYPE

                -- format fields
                , cast(format_GT as varchar) as format_GT
                , array_transform(format_GL, x -> cast(x as float)) as format_PL

            from filtered
        ),

        final as (select * from cleaned)
    select * from final
) to '{{ output }}' (format parquet);