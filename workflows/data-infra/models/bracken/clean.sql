copy (
    with
        bracken as (
            select "sample"
                , "name"
                , new_est_reads
                , fraction_total_reads
            from read_parquet('{{ bracken_glob }}', hive_partitioning = false)
        ),

        top_species_id as (
            select "sample", max(fraction_total_reads) as max_frac
            from bracken
            group by "sample"
        ),

        final as (
            select * exclude("name")
                , regexp_replace(
                    trim("name"), ' ', '_', 'g'
                ) as reference_genome
            from bracken
            where ("sample", fraction_total_reads) in (
                select ("sample", max_frac) from top_species_id
            )
        )
    select * from final
) to '{{ output }}' (format csv);