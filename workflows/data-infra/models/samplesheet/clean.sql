copy (
    with
        fastqs as (
            select * from read_csv('{{ fastqs_glob }}', hive_partitioning = true)
        ),

        samples as (
            select * from read_csv('{{ samples_glob }}', hive_partitioning = true)
        ),

        sequencing as (
            select * from read_csv('{{ sequencing_glob }}', hive_partitioning = true)
        ),

        junction_table as (
            select distinct on(id, library) id, library
            from
                (
                    select id, library from fastqs
                    union
                    select id, library from samples
                    union
                    select id, library from sequencing
                )
        ),

        combined as (
            select j.*
                , fastqs.* exclude(id, library)
                , samples.* exclude(id, library)
                , sequencing.* exclude(id, library)
            from junction_table j
            full join fastqs
            on
                j.id = fastqs.id
                and j.library = fastqs.library
            full join samples
            on
                j.id = samples.id
                and j.library = samples.library
            full join sequencing
            on
                j.id = sequencing.id
                and j.library = sequencing.library
        ),

        combined_w_unique_sample as (
            select row_number() over () as "sample", * from combined
        ),

        cleaned_relationship as (
            select "sample"
                , case
                    -- baby (focal)
                    when relationship ilike 'baby%' then 'baby'
                    -- sibling
                    when relationship ilike 'kind%' then 'sibling'
                    when relationship ilike 'sibling%' then 'sibling'
                    -- mother
                    when relationship ilike 'mutter' then 'mother'
                    when relationship ilike 'mother' then 'mother'
                    -- father
                    when relationship ilike 'vater' then 'father'
                    when relationship ilike 'father' then 'father'
                    else null
                end as relationship
            from combined_w_unique_sample
        ),

        cleaned_timepoint_unit as (
            select "sample"
                , timepoint
                , case
                    -- before
                    when timepoint_unit ilike 'vor' then 'before'
                    when timepoint_unit ilike 'before' then 'before'
                    -- months
                    when timepoint_unit ilike 'm' then 'months'
                    when timepoint_unit ilike 'monate' then 'months'
                    when timepoint_unit ilike 'months' then 'months'
                    -- weeks
                    when timepoint_unit ilike 'w' then 'weeks'
                    when timepoint_unit ilike 'wochen' then 'weeks'
                    when timepoint_unit ilike 'weeks' then 'weeks'
                    else null
                end
                as collection_interval_unit
            from (
                    select "sample"
                        , timepoint
                        , regexp_extract(timepoint, '(\D+)$', 1) as timepoint_unit
                    from combined_w_unique_sample
                )
        ),

        cleaned_timepoint as (
            select "sample"
                , timepoint
                , concat_ws(
                    ' ', cast(timepoint_value as varchar), collection_interval_unit
                ) as collection_interval_category
                , cast(
                    round(
                        case
                            -- before
                            when collection_interval_unit = 'before' then -1
                            -- months, 1 day ~= (365 * 4 + 1) / (12 * 4)
                            when collection_interval_unit = 'months' then timepoint_value * ((365 * 4 + 1) / (12 * 4))
                            -- weeks
                            when collection_interval_unit = 'weeks' then timepoint_value * 7
                            else null
                        end
                    ) || ' days'
                    as interval
                ) as collection_interval
            from (
                select "sample"
                    , timepoint
                    , try_cast(
                        regexp_extract(timepoint, '^(\d+)', 1)
                        as usmallint
                    ) as timepoint_value
                    , collection_interval_unit
                from cleaned_timepoint_unit
            )
        ),

        cleaned_species as (
            select "sample"
                , species
                , case
                    -- Bacteroides
                    when species ilike 'bacteroides%' then 'Bacteroides_xylanisolvens'
                    -- Bifidobacterium spp
                    when species ilike 'bifidobacterium%' then 'Bifidobacterium_spp'
                    -- Enterococcus faecalis
                    when species ilike '%enterococcus%' then 'Enterococcus_faecalis'
                    -- Escherichia coli
                    when species ilike 'escherichia%' then 'Escherichia_coli'
                    -- Klebsiella oxytoca
                    when species ilike 'klebsiella%' then 'Klebsiella_oxytoca'
                    -- Staphylococcus aureus
                    when species ilike 'staphylococcus%' then 'Staphylococcus_aureus'
                    else null
                end as species_stdized
            from combined_w_unique_sample
        ),

        final as (
            select t1.family
                , t2.relationship
                , t1.donor
                , t1.id
                , t3.timepoint
                , t3.collection_interval_category
                , t3.collection_interval
                , t4.species
                , t4.species_stdized
                , t1.library
                , t1.plate
                , t1.well
                , t1.barcode_1
                , t1.barcode_2
                , t1.notes
                , t1.fastq_1
                , t1.fastq_2
            from combined_w_unique_sample t1
            join cleaned_relationship t2 on (t1.sample = t2.sample)
            join cleaned_timepoint t3 on (t1.sample = t3.sample)
            join cleaned_species t4 on (t1.sample = t4.sample)
            where fastq_1 is not null
            and fastq_2 is not null
        )
    select * from final
) to '{{ output }}' (format csv);