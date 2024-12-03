set preserve_insertion_order = false;

create table annotated_vcfs as

with raw_vcfs as (
    
    select
        "sample"
        , family
        , species
        , chromosome
        , position
        , reference
        , alternate
        , quality
        , info_ADF
        , info_ADR
        , info_AD
        -- TODO: Which column for `case when INFO ilike '%INDEL%' then true else false end as indel`
        -- , info
    
    from
        read_parquet(
            getenv('VCFS'), 
            hive_partitioning = true
            -- TODO: Not sure whether hive paritioning recognizes non-asterisk partitions
            -- hive_types = {
            --     'species': varchar,
            --     'family': varchar,
            --     'sample': uinteger
            -- }
        )

), cleaned as (

    select
        * exclude(alternate, info_ADF, info_ADR)
        , array_transform(
            alternate, x -> nullif(x, '')
        ) as alternate
        , array_transform(
            info_ADF, x -> cast(x as usmallint)
        ) as info_ADF
        , array_transform(
            info_ADR, x -> cast(x as usmallint)
        ) as info_ADR

    from
        raw_vcfs

), final as (

    select * from cleaned

)

select * from final;