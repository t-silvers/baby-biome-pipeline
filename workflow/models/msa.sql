set enable_progress_bar = false;


create temp table filtered_calls as

select * from read_parquet(getenv('FILTERED_CALLS'));


create temp table with_outgroup as

with min_gap_samples as (
    
    select id from filtered_calls
    
    group by id having count(*) > 0.8 * (
        
        select count(distinct(chromosome, position)) from filtered_calls

    )

),

samples_passing as (

    select * from filtered_calls where id in (select id from min_gap_samples)

),

cartesian_template as (

    select s.id, c.chromosome, c.position from
        (
            
            select distinct chromosome, position from samples_passing
            
        ) c
    
    cross join (select distinct id from samples_passing) s

),

reference_outgroup as (

    select
        'reference' as id
        , c.chromosome
        , c.position
        , v.allele
    
    from
        (
            select
                distinct on(chromosome, position)
                chromosome
                , position
            
            from cartesian_template

        ) c
    
    left join
        (
            
            select
                distinct on(chromosome, position)
                chromosome
                , position
                , reference as allele
            
            from samples_passing

        ) v
    
    on
        c.chromosome = v.chromosome
        and c.position = v.position

),

imputed as (

    select
        template.id
        , template.chromosome
        , template.position
        , calls.allele
    
    from cartesian_template template
    
    left join samples_passing calls

    on
        template.id = calls.id
        and template.chromosome = calls.chromosome
        and template.position = calls.position

),

add_outgroup as (

    select
        * exclude(allele)
        , ifnull(allele, 'N') as allele
    
    from imputed

    where id is not null
    
    union by name
        (
            select * from reference_outgroup
        )
    
    order by id, chromosome, position

) select * from add_outgroup;


create temp table annot_msa as

select s.species, s.family, s.id, s.library, v.*

from with_outgroup v

left join
    (

        select distinct on(id)
            * exclude(chromosome, position, reference, allele)
        
        from filtered_calls

    ) s

on cast(v.id as varchar) = cast(s.id as varchar);


select * from final;

copy (

    select * from annot_msa

) to '/dev/stdout' (format parquet);