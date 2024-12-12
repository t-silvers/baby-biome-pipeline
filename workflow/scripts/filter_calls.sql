create temp table filtered_calls as

with dp_filtered as (

    select
        row_number() over () as idx
        , * exclude(is_forward, dp)
        , sum(dp) as dp 
    
    from prioritized 
    
    group by species, family, id, library, chromosome, position, reference, allele 
    
    having sum(dp) >= cast(getenv('DP') as int)

),

ref_or_alt as (

    select a.idx, a.allele from
        
        (

            select idx, allele, dp, row_number() over (
                    
                partition by idx order by dp desc
            
            ) rnkd
            
            from dp_filtered
        
        ) a

    where a.rnkd = 1 

),

final as (

    select d.* exclude(idx, allele, dp), r.allele
    
    from ref_or_alt r

    join dp_filtered d on (r.idx = d.idx and r.allele = d.allele)

)

select * from final;

copy (

    select * from filtered_calls

) to '/dev/stdout' (format parquet);