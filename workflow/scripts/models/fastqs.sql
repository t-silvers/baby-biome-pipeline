create temp table fastqs as

with fastqs_glob as (

    select
        "file"
    
    from
        glob(getenv('GLOB'))

),

parsed_filename as (

    select
        "file"
        , extracted.* exclude ("read")
        , 'fastq_' || extracted.read as "read"
    
    from
        (
            select
                "file"
                , regexp_extract(
                    "file",
                    getenv('PAT'),
                    ['library', 'id', 'family', 'read']
                ) as extracted
            
            from
                fastqs_glob

        )

),

pivot_on_reads as (

    pivot
        parsed_filename
        
    on
        "read" 
    
    using
        first("file") 
    
    group by
        library
        , id

),

final as (

    select
        * exclude(id)
        , replace(id, '-', '_') as id
    
    from
        pivot_on_reads
    
    order by
        id

)

select * from final;

copy (

    select
        id
        , columns('fastq_[1,2]')

    from
        fastqs

) to '/dev/stdout' (format csv);