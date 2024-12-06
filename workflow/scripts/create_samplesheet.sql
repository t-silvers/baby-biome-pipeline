create table samplesheet as

with samples as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^resources/library=(.*)/samples.csv$', 1) as library

    from
        read_csv(
            'resources/library=*/samples.csv',
            filename=true
        )

), sequencing as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^resources/library=(.*)/sequencing.csv$', 1) as library

    from
        read_csv(
            'resources/library=*/sequencing.csv',
            filename=true
        )

), fastqs as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^resources/library=(.*)/fastqs.csv$', 1) as library

    from
        read_csv(
            'resources/library=*/fastqs.csv',
            filename=true
        )

), junction_table as (

    select
        distinct on(id, library) id, library 
    
    from
        (

            select id, library from samples
            
            union
            select id, library from sequencing
            
            union
            select id, library from fastqs

        )

), combined as (

    select
        j.*
        , samples.* exclude(id, library)
        , sequencing.* exclude(id, library)
        , fastqs.* exclude(id, library)

    from junction_table j
    
    full join samples
    
    on
        j.id = samples.id
        and j.library = samples.library

    full join sequencing
    
    on
        j.id = sequencing.id
        and j.library = sequencing.library

    full join fastqs
    
    on
        j.id = fastqs.id
        and j.library = fastqs.library

), combined_w_unique_sample as (

    select row_number() over () as "sample" from combined

), cleaned_relationship as (

    select
        "sample"
        , case
            
            -- baby (focal)
            when relationship ilike 'baby%' then 'B'
            
            -- sibling
            when relationship ilike 'kind%' then 'S'
            when relationship ilike 'sibling%' then 'S'
            
            -- mother
            when relationship ilike 'mutter' then 'M'
            when relationship ilike 'mother' then 'M'
            
            -- father
            when relationship ilike 'vater' then 'F'
            when relationship ilike 'father' then 'F'
            
            else null
        end as relationship
    
    from combined_w_unique_sample

), cleaned_timepoint as (

    select
        "sample"
        , timepoint
        , cast(
            round(
                case
                    -- before
                    when timepoint_unit ilike 'vor' then -1
                    when timepoint_unit ilike 'before' then -1

                    -- months
                    when timepoint_unit ilike 'm' then timepoint_value * 30.4
                    when timepoint_unit ilike 'monate' then timepoint_value * 30.4
                    when timepoint_unit ilike 'months' then timepoint_value * 30.4
                    
                    -- weeks
                    when timepoint_unit ilike 'w' then timepoint_value * 7
                    when timepoint_unit ilike 'wochen' then timepoint_value * 7
                    when timepoint_unit ilike 'weeks' then timepoint_value * 7
                    
                    else null
                end
            ) || ' days'
            as interval
        ) as specimen_collection_interval
    
    from
        (
            
            select
                "sample"
                , timepoint
                , try_cast(
                    regexp_extract(
                        timepoint,
                        '^(\d+)',
                        1
                    )
                    as usmallint
                ) as timepoint_value
                , regexp_extract(
                    timepoint,
                    '(\D+)$',
                    1
                ) as timepoint_unit

            from combined_w_unique_sample

        )

), cleaned_species as (

    select
        "sample"
        , species
        , case
            
            -- Bacteroides
            when species ilike 'bacteroides%' then 'Bacteroides_xylanisolvens'
            
            -- Bifidobacterium spp
            when species ilike 'bifidobacterium%' then 'Bifidobacterium_spp'
            
            -- Escherichia coli
            when species ilike 'escherichia%' then 'Escherichia_coli'
            
            -- Enterococcus faecalis
            when species ilike '%enterococcus%' then 'Enterococcus_faecalis'
            
            -- Klebsiella oxytoca
            when species ilike 'klebsiella%' then 'Klebsiella_oxytoca'
            
            -- Staphylococcus aureus
            when species ilike 'staphylococcus%' then 'Staphylococcus_aureus'
            
            else null
        end as species_stdized

    from combined_w_unique_sample

), final as (

    select 
        t1.sample
        , t1.family
        , t2.relationship
        , t1.donor
        , t1.id
        , t3.timepoint
        , t3.specimen_collection_interval
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

    from  with_sample_and_donor_id t1

    join cleaned_relationship t2 on (t1.sample = t2.sample)

    join cleaned_timepoint t3 on (t1.sample = t3.sample)

    join cleaned_species t4 on (t1.sample = t4.sample)

)

select * from final;


copy (

    select * from samplesheet

) to '/dev/stdout' (format csv);