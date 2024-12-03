-- TODO: Rough script. Should refactor cleans to select subset of table.

create table samplesheet as

with sample_info as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^results/sample_info_(.*).csv$', 1) as library

    from
        read_csv(
            'results/sample_info_*.csv',
            filename=true
        )

), seq_info as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^results/seq_info_(.*).csv$', 1) as library

    from
        read_csv(
            'results/seq_info_*.csv',
            filename=true
        )

), fastqs as (

    select
        * exclude("filename")
        , regexp_extract("filename", '^results/fastqs_(.*).csv$', 1) as library

    from
        read_csv(
            'results/fastqs_*.csv',
            filename=true
        )

), junction_table as (

    select
        distinct on(id, library) id, library 
    
    from
        (

            select id, library from sample_info
            
            union
            select id, library from seq_info
            
            union
            select id, library from fastqs

        )

), combined as (

    select
        j.*
        , sample_info.* exclude(id, library)
        , seq_info.* exclude(id, library)
        , fastqs.* exclude(id, library)

    from
        junction_table j
    
    full join
        sample_info
    
    on
        j.id = sample_info.id
        and j.library = sample_info.library

    full join
        seq_info
    
    on
        j.id = seq_info.id
        and j.library = seq_info.library

    full join
        fastqs
    
    on
        j.id = fastqs.id
        and j.library = fastqs.library

), cleaned_relationship as (

    select
        * exclude (relationship)
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
    
    from
        combined

), with_sample_and_donor_id as (

    select
        row_number() over () as "sample"
        , * exclude(donor)
        , family || '_' || relationship as donor
    
    from
        cleaned_relationship

), cleaned_timepoint as (

    select
        * exclude(timepoint_num, timepoint_unit)
        , cast(
            case
                -- before
                when timepoint_unit ilike 'vor' then -1
                when timepoint_unit ilike 'before' then -1

                -- months
                when timepoint_unit ilike 'm' then timepoint_num * 30.4
                when timepoint_unit ilike 'monate' then timepoint_num * 30.4
                when timepoint_unit ilike 'months' then timepoint_num * 30.4
                
                -- weeks
                when timepoint_unit ilike 'w' then timepoint_num * 7
                when timepoint_unit ilike 'wochen' then timepoint_num * 7
                when timepoint_unit ilike 'weeks' then timepoint_num * 7
                
                else null
            end
            as smallint
        ) as timepoint_day
    
    from
        (
            
            select
                *
                , try_cast(
                    regexp_extract(
                        timepoint,
                        '^(\d+)',
                        1
                    )
                    as usmallint
                ) as timepoint_num
                , regexp_extract(
                    timepoint,
                    '(\D+)$',
                    1
                ) as timepoint_unit

            from
                with_sample_and_donor_id

        )

), cleaned_species as (

    select
        *
        , case
            
            -- Escherichia coli
            when species ilike 'escherichia%' then 'Escherichia_coli'
            
            -- Bacteroides
            when species ilike 'bacteroides%' then 'Bacteroides_xylanisolvens'
            
            -- Enterococcus faecalis
            when species ilike '%enterococcus%' then 'Enterococcus_faecalis'
            
            -- Staphylococcus aureus
            when species ilike 'staphylococcus%' then 'Staphylococcus_aureus'
            
            -- Bifidobacterium spp
            when species ilike 'bifidobacterium%' then 'Bifidobacterium_spp'
            
            -- Klebsiella oxytoca
            when species ilike 'klebsiella%' then 'Klebsiella_oxytoca'
            
            else null
        end as species_stdized

    from cleaned_timepoint

), final as (

    select * from cleaned_species

)

select * from final;


copy (

    select * from samplesheet

) to '/dev/stdout' (format csv);