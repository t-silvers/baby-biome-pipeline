set enable_progress_bar = false;

copy (

    select
        '>' || id || chr(10) || string_agg(allele, '')
    
    from
        (
            select
                id, allele
            
            from
                annot_msa
            
            order by
                id
                , chromosome
                , position
        )
    
    group by
        id
    
    order by
        id

) to '/dev/stdout' (
    delimiter '',
    header false,
    quote ''
);