load spatial;

create table sample_info as 

with raw_sample_info as (

    select * from st_read(getenv('FN'), open_options = ['HEADERS=FORCE'])
),

cleaned as (

    select
        ID as id
        , Family as family
        , "Subject" as relationship
        , Timepoint as timepoint

    from raw_sample_info

    where id != 'control' and id != 'empty_well'

    -- To have sample field ~match~ id field
    order by id

),

with_sample_and_donor_id as (

    select
        row_number() over () as "sample"
        , *
        , family || '_' || relationship as donor_id
    
    from cleaned

),

cleaned_relationship as (

    select
        "sample"
        , case
            when relationship ilike 'Baby%' then 'B'
            when relationship ilike 'Kind%' then 'K'
            when relationship ilike 'Mutter' then 'M'
            when relationship ilike 'Vater' then 'F'
            else null
        end as relationship
    
    from with_sample_and_donor_id

),

cleaned_timepoint as (

    select
        "sample"
        , timepoint_cat
        , cast(
            case
                when timepoint_unit = 'vor' then -1
                when timepoint_unit = 'M' then timepoint_num * 30.4
                when timepoint_unit = 'W' then timepoint_num * 7
                else null
            end
            as smallint
        ) as timepoint_day
    
    from
        (
            
            select
                "sample"
                , timepoint as timepoint_cat
                , try_cast(
                    regexp_extract(
                        timepoint,
                        '^(\d+)(M|W)$',
                        1
                    )
                    as usmallint
                ) as timepoint_num
                , regexp_extract(
                    timepoint,
                    '(\D+)$',
                    1
                ) as timepoint_unit

            from with_sample_and_donor_id

        )

),

final as (

    select 
        t1.sample
        , t1.id
        , t1.donor_id
        , t1.family
        , t2.relationship
        , t3.timepoint_cat
        , t3.timepoint_day

    from  with_sample_and_donor_id t1

    inner join 
        cleaned_relationship t2 on t1.sample = t2.sample

    inner join 
        cleaned_timepoint t3 on t2.sample = t3.sample

)

select * from final order by "sample";