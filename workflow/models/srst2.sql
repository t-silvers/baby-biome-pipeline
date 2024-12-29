copy (
    with
        raw_srst2_results as (
            select * from read_csv(
                '{{ input }}',
                delim = '\t',
                header = true,
                columns = {
                    'Sample': 'varchar',
                    'ST': 'varchar',
                    'adk': 'varchar',
                    'fumC': 'varchar',
                    'gyrB': 'varchar',
                    'icd': 'varchar',
                    'mdh': 'varchar',
                    'purA': 'varchar',
                    'recA': 'varchar',
                    'mismatches': 'varchar',
                    'uncertainty': 'varchar',
                    'depth': 'double',
                    'maxMAF': 'double',
                },
                nullstr = '-',
                auto_detect = false
            )
        ),

        final as (
            select "Sample" as "sample"
                , try_cast(regexp_extract(ST, '\d+') as bigint) as st
                , case when contains(ST, '?') then true else false end as low_depth
                , case when contains(ST, 'NF') then true else false end as not_found
                , case when contains(ST, 'ND') then true else false end as not_done
                , case when ST = 'failed' then true else false end as failed
                , depth
                , maxMAF as max_maf
            from raw_srst2_results
        )
    select * from final
) to '{{ output }}' (delimiter '\t');