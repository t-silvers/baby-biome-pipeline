set preserve_insertion_order = false;
set max_temp_directory_size = '{{ max_temp_directory_size }}';

create temp table new_variants as
select * from read_parquet(
    '{{ vcf_pq_glob }}',
    union_by_name = true,
    hive_partitioning = false
);

-- NOTE: Relies on automatic type conversion for "contig" field
insert or ignore into variants
    by name
    select * from new_variants;

copy (
    with
        snvs as (
            select "sample", contig, start_pos, alleles
            from variants 
            where not list_contains(info_TYPE, 'indel')
        ),
        
        shared_positions as (
            select contig, start_pos
            from snvs
            group by contig, start_pos
            having count(*) >= cast('{{ frac_cohort_share_pos }}' as float) * (
                select count(distinct "sample")
                from variants
            )
        ),

        shared_snvs as (
            select * from snvs
            where (contig, start_pos) in (
                select (contig, start_pos)
                from shared_positions
            )
        ),

        final as (
            select "sample"
                , contig
                , start_pos as position
                , alleles[1] as allele
            from shared_snvs
        )
    select * from final
) to '{{ output }}' (format parquet);