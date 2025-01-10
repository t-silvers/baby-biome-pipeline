load spatial;

copy (
    with
        raw_sample_info_sheet1 as (
            select * from st_read(
                '{{ input }}',
                layer = 'LibraryPlate1-3',
                open_options = ['HEADERS=FORCE']
            )
        ),

        raw_sample_info_sheet2 as (
            select * from st_read(
                '{{ input }}',
                layer = 'LibraryPlate4-8'
            )
        ),

        cleaned_sheet1 as (
            select "gDNA from 96-Well" as id
                , Field2 as species
                , regexp_extract(Field3, '^(\D+)[- ]', 1) as relationship
                , regexp_extract(Field3, '[- ](\w+)$', 1) as timepoint
            from raw_sample_info_sheet1
        ),

        cleaned_sheet2 as (
            select Field1 as id
                , Field2 as gDNA
                , case
                    when Field4 is null then Field3
                    else concat(Field3, ' (', Field4, ')')
                end as species
                , Field5
                , Field6 as relationship
                , Field7 as timepoint
            from raw_sample_info_sheet2
        ),

        cleaned_combined as (
            select * exclude(species)
                , regexp_replace(trim(species), ' ', '_', 'g') as species
            from (select * from cleaned_sheet1)
            union by name
            select * from cleaned_sheet2
        ),

        final as (
            select id
                , regexp_extract(id, '([BP]\d+|Ctr\d+|Control\d+)_\d+', 1) as family
                , family || '_' || relationship as donor
                , relationship
                , timepoint
                , species
            from cleaned_combined
            where id is not null
        )

    select * from final
) to '{{ output }}' (format csv);