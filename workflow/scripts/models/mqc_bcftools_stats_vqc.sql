set preserve_insertion_order = false;

-- See https://github.com/MultiQC/MultiQC/blob/a70cdd4e60eab452b9fabdd14f035e7483fe3a0b/multiqc/modules/bcftools/stats.py#L192-L202

create table variant_quality_scores as

with raw_multiqc_bcftools_vqs as (

    select
        "sample"
        , quality
        , count

    from
        read_csv(
            '/dev/stdin',
            auto_detect = false,
            columns = {
                'sample': 'varchar',
                'quality': 'bigint',
                'count': 'bigint',
            },
            header = false
        )

)

select * from raw_multiqc_bcftools_vqs;