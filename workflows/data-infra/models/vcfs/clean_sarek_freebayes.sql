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
                , info_EPP
                , info_EPPR
                , info_MQM
                , info_MQMR
                , info_SAP
                , info_SRP
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
            where quality >= cast('{{ quality }}' as float)
                and info_DP >= cast('{{ dp }}' as float)
                and (info_AD[1] / array_reduce(info_AD, (x, y) -> x + y)) < (1 - cast('{{ maf }}' as float))
                and array_reduce(info_AD, (x, y) -> x + y) >= cast('{{ ad }}' as float)
                and list_min(info_EPP) < cast('{{ epp }}' as float)
                and info_EPPR < cast('{{ epp }}' as float)
                and list_max(info_MQM) >= cast('{{ mq }}' as float)
                and info_MQMR >= cast('{{ mq }}' as float)
                and (
                    abs(list_max(info_MQM) - info_MQMR) < 10
                    and (list_max(info_MQM) / (info_MQMR + .1)) >= .4
                    and (list_max(info_MQM) / info_MQMR) <= 2.5
                )
                and list_min(info_SAP) < cast('{{ sp }}' as float)
                and info_SRP < cast('{{ sp }}' as float)
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
                , array_transform(
                    array_concat(
                        [info_EPPR], info_EPP),
                        x -> cast(x as decimal(4, 1)
                    )
                ) as info_EPP
                , array_transform(
                    array_concat(
                        [info_MQMR], info_MQM),
                        x -> cast(x as decimal(4, 1)
                    )
                ) as info_MQM
                , array_transform(
                    array_concat(
                        [info_SRP], info_SAP),
                        x -> cast(x as decimal(4, 1)
                    )
                ) as info_SP
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