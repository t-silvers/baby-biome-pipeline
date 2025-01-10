create temp table parsed_sequencing_notes as 
select * exclude(id)
    , id as original_id
    , case
        when notes ilike '%fastqID%' then regexp_extract(
            notes,
            'fastqID_of_sample:(plate\d+_[A-L]\d+|(B\d+)_\d+)$',
            1
        )
        else original_id 
    end as id
from read_csv('{{ sequencing }}');


copy (
    with incorrect_fastqs as (
        select * from read_csv('{{ fastqs }}')
    ),

    final as (
        select coalesce(t2.id, t1.id) as id
            , columns('fastq_[1,2]')
        from incorrect_fastqs t1
        left join parsed_sequencing_notes t2
        on t1.id = t2.original_id
    )

    select * from final
) to '{{ f_output }}' (format csv);


copy (
    with incorrect_sequencing as (
        select * from read_csv('{{ sequencing }}')
    ),

    final as (
        select coalesce(t2.id, t1.id) as id
            , t1.* exclude(id)
        from incorrect_sequencing t1
        left join parsed_sequencing_notes t2
        on t1.id = t2.original_id
    )

    select plate
        , well
        , barcode_1
        , barcode_2
        , notes
        , id
        , library
    from final
) to '{{ s_output }}' (format csv);