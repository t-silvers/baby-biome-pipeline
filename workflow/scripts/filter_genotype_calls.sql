create temp table filtered_calls as

with family_filtered as (

    select
        *
    
    from
        annotated_vcfs
    
    where
        family = getenv('FAMILY')

), 

snvs_only as (

    select 
        * exclude(info_INDEL)

    from
        family_filtered

    where
        not info_INDEL

),

-- simple filters

quality_depth_filtered as (

    select 
        * exclude(quality, info_SNPGAP, info_INDELGAP_3, indel_INDELGAP_10, format_SP)

    from
        snvs_only

    where
        quality >= cast(getenv('QUAL') as float)
        and list_reduce(
            info_AD, (x, y) -> x + y
        ) >= cast(getenv('DP') as int)
        and info_SNPGAP > cast(getenv('SNPGAP') as int)
        and ((indel_INDELGAP_10 != true) or (indel_INDELGAP_10 is null))
        and format_SP < list_reduce(
            info_AD, (x, y) -> x + y
        ) / 2

), 

coalesced_alternates as (

    select
        * exclude(alternate)
        , alternate[1] as alternate
    
    from
        quality_depth_filtered

),

maf_readbias_filtered as (

    select 
        species
        , family
        , relationship
        , donor
        , id
        , timepoint
        , timepoint_day
        , chromosome
        , position
        , reference
        , case 
            when alternate is null then reference 
            when info_ADF[2] >= cast(getenv('STRAND_DP') as int)
                and info_ADR[2] >= cast(getenv('STRAND_DP') as int)
                and list_transform(info_AD, x -> x / list_aggregate(info_AD, 'sum'))[2] > cast(getenv('MAF') as float)
                then alternate
            when info_ADF[1] >= cast(getenv('STRAND_DP') as int)
                and info_ADR[1] >= cast(getenv('STRAND_DP') as int)
                and list_transform(info_AD, x -> x / list_aggregate(info_AD, 'sum'))[1] > cast(getenv('MAF') as float)
                then reference
            else null
        end as allele

    from
        coalesced_alternates

    where
        allele is not null

),

variable_positions as (

    select
        chromosome, position

    from
        maf_readbias_filtered

    group by
        chromosome, position

    having
        count(distinct allele) > 1
        and count(id) > .95 * (
            select 
                count(distinct id) 
            from 
                maf_readbias_filtered
        )

),

variable_filtered as (

    select
        species
        , family
        , relationship
        , donor
        , id
        , timepoint
        , timepoint_day
        , chromosome
        , position
        , reference
        , allele

    from
        maf_readbias_filtered

    where
        (chromosome, position) in (
            select
                (chromosome, position)
            from
                variable_positions
        )

    order by
        species
        , family
        , id
        , chromosome
        , position

),

filter_samples as (
    
    select
        id
    
    from
        variable_filtered

    group by
        id

    having
        count(*) > .95 * (
            select 
                count(distinct position) 
            from 
                variable_filtered
        )

),

final as (
    
    select
        species
        , family
        , relationship
        , donor
        , id
        , timepoint
        , timepoint_day
        , chromosome
        , position
        , reference
        , allele
    
    from
        variable_filtered

    where
        id in (
            
            select id from filter_samples
        )

    order by
        species
        , family
        , id
        , chromosome
        , position

)

select * from final;

copy (

    select
        *

    from
        filtered_calls

) to '/dev/stdout' (format parquet);