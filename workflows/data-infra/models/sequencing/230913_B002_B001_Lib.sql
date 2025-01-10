copy (
    with
        raw_seq_info as (
            select * from read_csv(
                '{{ input }}',
                auto_detect = false,
                columns = {
                    'Original sample plate (or first arrayed plate, if source is tubes)': varchar,
                    'Well': varchar,
                    'Library plate name': varchar,
                    'Library plate well': varchar,
                    'Barcode 1 (Plate or Row)': varchar,
                    'Barcode 2 (Well or Column)': varchar,
                    'SampleName (optional)': varchar,
                    'Library prep method (optional)': varchar,
                    'Notes (optional)': varchar
                },
                skip = 13
            )
        ),

        cleaned as (
            select
                case when "SampleName (optional)" ilike 'control%' then 'control'
                    else "SampleName (optional)"
                end as id
                , cast(
                    regexp_extract(
                        "Library plate name",
                        '^Library Plate (\d+)$',
                        1
                    )
                    as usmallint
                ) as plate
                , "Library plate well" as well
                , "Barcode 1 (Plate or Row)" as barcode_1
                , "Barcode 2 (Well or Column)" as barcode_2
                , "Notes (optional)" as notes
            from raw_seq_info
        ),

        final as (
            select * exclude(id)
                , regexp_replace(id, '-', '_') as id
            from cleaned
            where id is not null
        )

    select * from final
) to '{{ output }}' (format csv);