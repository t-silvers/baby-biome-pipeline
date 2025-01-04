create temp table new_progress as
with
    raw_progress as (
        select info.* from (
            select
                regexp_extract(
                    "file",
                    'variants/tool=(bactmap|legacy_mapping|snippy|sarek_bcftools|sarek_deepvariant|sarek_freebayes|sarek_haplotypecaller)/.*/(\d+).cleaned.parquet$',
                    ['tool', 'sample']
                ) as info
            from glob('variants/tool=*/**/*.cleaned.parquet')
        )
        where info.sample != ''
    ),

    final as (
        select "sample"
            , array_sort(array_agg(distinct tool)) as tool
        from raw_progress
        group by "sample"
    )
select * from final;

insert or replace into mapping_progress
    by name
    select * from new_progress;

copy (
    select "sample", unnest(tool) as tool
    from mapping_progress
) to '{{ output }}' (format csv);