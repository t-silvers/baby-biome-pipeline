create temp table samplesheet as

with fastqs as (

    select
        id
        , columns('fastq_[1,2]')

    from
        read_csv('results/fastqs.csv')

),

final as (

    select
        sample_info.sample
        , columns('fastq_[1,2]')
    
    from
        fastqs

    left join 
        sample_info on fastqs.id = sample_info.id

    where
        sample_info.sample is not null
    
    order by
        sample_info.sample

)

select * from final;

copy (

    select 
        *

    from
        samplesheet

) to '/dev/stdout' (format csv);