create temp table filtered_calls as

with snvs_only as (

    select 
        * exclude(indel)

    from
        annotated_vcfs

    where
        indel

),

quality_depth_filtered as (

    select 
        * exclude(quality)

    from
        snvs_only

    where
        quality >= cast(getenv('QUAL') as float)
        and list_reduce(
            allelic_read_depth, (x, y) -> x + y
        ) >= cast(getenv('DP') as int)

),

maf_readbias_filtered as (

    select 
        "sample"
        , chromosome
        , position
        , case 
            when alternate is null then reference 
            when allelic_read_depth_forward[2] >= cast(getenv('STRAND_DP') as int)
                and allelic_read_depth_reverse[2] >= cast(getenv('STRAND_DP') as int)
                and list_transform(allelic_read_depth, x -> x / list_aggregate(allelic_read_depth, 'sum'))[2] >= cast(getenv('MAF') as float)
                then alternate
            when allelic_read_depth_forward[1] >= cast(getenv('STRAND_DP') as int)
                and allelic_read_depth_reverse[1] >= cast(getenv('STRAND_DP') as int)
                and list_transform(allelic_read_depth, x -> x / list_aggregate(allelic_read_depth, 'sum'))[1] >= cast(getenv('MAF') as float)
                then reference
            else null
        end as allele

    from
        quality_depth_filtered

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
        and count("sample") > .95 * (
            select 
                count(distinct "sample") 
            from 
                maf_readbias_filtered
        )

),

final as (

    select
        *

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
        "sample"
        , chromosome
        , position

)

select * from final;

copy (

    select
        *

    from
        filtered_calls

) to '/dev/stdout' (format csv);