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
                , regexp_extract(ID, '([BP]\d+|Ctr\d+|Control\d+)_\d+', 1) as family
                , family || '_' || "subject" as donor
                , "subject" as relationship
                , timepoint
                , species
            from raw_sample_info
        ),

        final as (
            select  id
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