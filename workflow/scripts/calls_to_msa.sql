set enable_progress_bar = false;

create temp table msa as

with annot_filtered_calls as (

    select
        *

    from
        read_csv(getenv('FILTERED_CALLS'))

),

min_gap_samples as (
    
    select
        "sample"
    
    from
        annot_filtered_calls
    
    group by
        "sample"
    
    having
        count(*) > 0.8 * (
            
            select
                count(distinct(chromosome, position))
            
            from
                annot_filtered_calls
        )

),

samples_passing as (

    select
        *

    from
        annot_filtered_calls

    where
        "sample" in (

            select
                "sample"
            
            from
                min_gap_samples

        )

),

cartesian_template as (

    select
        s.sample,
        c.chromosome,
        c.position

    from
        (

            select
                distinct chromosome, position 
            
            from
                samples_passing

        ) c
    
    cross join
        (
            
            select
                distinct "sample"
            
            from
                samples_passing
        
        ) s

),

imputed as (

    select
        template.sample
        , template.chromosome
        , template.position
        , calls.allele
    
    from
        cartesian_template template
    
    left join
        samples_passing calls

    on
        template.sample = calls.sample
        and template.chromosome = calls.chromosome
        and template.position = calls.position

),

final as (

    select
        * exclude(allele)
        , ifnull(allele, 'N') as allele
    
    from
        imputed
    
    where
        "sample" is not null
    
    order by
        "sample"
        , chromosome
        , position

)

select * from final;


-- Add annotations
create temp table annot_msa as

with sample_annots as (
    
    select
        "sample"
        , id
    
    from
        sample_info

),
final as (
    
    select
        sample_annots.id
        , msa.chromosome
        , msa.position
        , msa.allele
    
    from
        msa
    
    left join
        sample_annots
    
    on
        msa.sample = sample_annots.sample
    
    order by
        sample_annots.id
        , msa.chromosome
        , msa.position

)

select * from final;


-- Output as fasta
copy (

    select
        '>' || id || chr(10) || string_agg(allele, '')
    
    from
        annot_msa
    
    group by
        id
    
    order by
        id

) to '/dev/stdout' (
    delimiter '',
    header false,
    quote ''
);