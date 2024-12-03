create temp table matched_reference_genome as

with bracken_parsed as (

    select
        "sample"
        , "name"
        , new_est_reads
        , fraction_total_reads

    from
        read_csv(getenv('BRACKEN'))

),

sample_info_parsed as (

    select
        "sample"
        , id
    
    from
        samplesheet

),

joined as (

    select
        s.*
        , b.* exclude("sample")
    
    from
        bracken_parsed b
    
    left join
        sample_info_parsed s on b.sample = s.sample

),

top_species_id as (

    select
        "sample"
        , max(fraction_total_reads) as max_frac
    
    from
        joined
    
    where
        new_est_reads > pow(
            10,
            cast(getenv('READ_POW') as int)
        )
        and fraction_total_reads > cast(getenv('READ_FRAC') as float)
    
    group by "sample"

),

final as (

    select
        "sample"
        , id
        , regexp_replace(
            trim("name"), ' ', '_', 'g'
        ) as species
        , new_est_reads
        , fraction_total_reads

    from
        joined
    
    where
        ("sample", fraction_total_reads) in (

            select
                ("sample", max_frac)
            
            from
                top_species_id
        
        )

)

select * from final;

copy (

    select 
        "sample"
        , species

    from
        matched_reference_genome
    
    order by
        cast("sample" as int)

) to '/dev/stdout' (format csv);