copy (
    with
        fastqs_glob as (
            select * from read_csv('{{ input }}')
        ),

        parsed_filename as (
            select "file"
                , extracted.* exclude ("read")
                , 'fastq_' || extracted.read as "read"
            from (
                select "file"
                    , regexp_extract(
                        "file",
                        '{{ pat }}',
                        ['library', 'id', 'family', 'read']
                    ) as extracted
                from fastqs_glob
            )
        ),

        pivot_on_reads as (
            pivot parsed_filename
            on "read" 
            using first("file") 
            group by library, family, id
        ),

        final as (
            select * exclude(id)
                , replace(id, '-', '_') as id
            from pivot_on_reads
        )

    select id, columns('fastq_[1,2]') from final
) to '{{ output }}' (format csv);