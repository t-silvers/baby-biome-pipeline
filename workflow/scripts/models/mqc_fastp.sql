set preserve_insertion_order = false;

create table fastp_total_reads as

with raw_multiqc_fastp as (

    select
        "sample"
        , before_filtering_total_reads
        , after_filtering_total_reads

    from
        read_csv(
            '/dev/stdin',
            auto_detect = false,
            columns = {
                'sample': 'varchar',
                'before_filtering_total_reads': 'bigint',
                'after_filtering_total_reads': 'bigint',
            },
            header = false
        )

)

select * from raw_multiqc_fastp;