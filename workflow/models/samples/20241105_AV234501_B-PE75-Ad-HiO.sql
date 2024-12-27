load spatial;

copy (
    with
        raw_sample_info as (
            select * from st_read(
                '{{ input }}',
                open_options = ['HEADERS=FORCE']
            )
        ),

        cleaned as (
            select ID as id
                , Family as family
                , Family || '_' || "Subject" as donor
                , "Subject" as relationship
                , Timepoint as timepoint
                , 'Escherichia_coli' as species
            from raw_sample_info
        ),

        final as (
            select id
                , family
                , donor
                , relationship
                , timepoint
                , species
            from  cleaned
            where id is not null
        )

    select * from final
) to '{{ output }}' (format csv);