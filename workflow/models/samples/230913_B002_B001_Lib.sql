load spatial;

copy (
    with
        raw_sample_info as (
            select * from st_read(
                -- '{{ input }}',
                '/nexus/posix0/MPIIB-keylab/lab_data/infant_mbiome/resources/library=240704_B002_B001_Lib_AVITI_reseq/samples-raw.ext',
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