set preserve_insertion_order = false;

create table samtools_stats as

with raw_multiqc_fastp as (

    select
        "sample"
        , reads_mapped_and_paired_percent
        , reads_properly_paired_percent

    from
        read_csv(
            '/dev/stdin',
            columns = {
                'sample': int,
                'reads_mapped_and_paired_percent': float,
                'reads_properly_paired_percent': float
            }
        )

)

select * from raw_multiqc_bcftools_vqs;