copy (
    with
        sample_info as (
            select "sample", family, id, library, species_stdized
            from read_csv('{{ samplesheet }}')
        ),

        bracken_output as (
            select "sample"
                , family
                , id
                , library
                , "name"
                , new_est_reads
                , fraction_total_reads
            from read_parquet('{{ bracken_glob }}')
        ),

        joined as (
            select s.*, b.* exclude("sample", family, id, library)
            from bracken_output b
            left join sample_info s
            on  b.sample = s.sample
            and b.family = s.family
            and b.id = s.id
            and b.library = s.library
        ),

        top_species_id as (
            select "sample"
                , max(fraction_total_reads) as max_frac
            from joined
            where new_est_reads > pow(10, cast('{{ read_pow }}' as int))
              and fraction_total_reads > cast('{{ read_frac }}' as float)
            group by "sample"
        ),

        final as (
            select "sample"
                , family
                , id
                , library
                , regexp_replace(
                    trim("name"), ' ', '_', 'g'
                ) as reference_genome
                , new_est_reads
                , fraction_total_reads
            from joined
            where ("sample", fraction_total_reads) in (
                    select ("sample", max_frac) from top_species_id
                  )
              and split_part("name", ' ', 1) = split_part(species_stdized, '_', 1)
            order by
                cast("sample" as int)
        )

    select * from final
) to '{{ output }}' (format csv);