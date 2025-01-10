copy (
    with
        samplesheet as (
            select * from read_csv('{{ input }}')
        ),

        completed_samples as (
            select "sample" from reference_genomes
            where "reference_genome" is not null
        ),

        taxprofiler_samplesheet as (
            select 
                "sample"
                , "sample" as run_accession
                , '{{ instrument_platform }}' as instrument_platform
                , fastq_1
                , fastq_2
            from samplesheet
        ),

        final as (
            select * from taxprofiler_samplesheet
            where "sample" not in (
                select "sample" from completed_samples
            )
        )
    select * from final
) to '{{ output }}' (format csv);