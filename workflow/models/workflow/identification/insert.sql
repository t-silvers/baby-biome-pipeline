create temp table new_progress as
with
    raw_progress as (
        select info.* from (
            select
                regexp_extract(
                    "file",
                    'identification/tool=(taxprofiler|srst2)/.*/(\d+).(bracken.parquet|results.txt)$',
                    ['tool', 'sample', 'ext']
                ) as info
            from glob('identification/tool=*/**/*')
        )
    ),

    final as (
        select "sample"
            , array_sort(array_agg(distinct tool)) as tool
        from raw_progress
        group by "sample"
    )
select * from final;

insert or replace into identification_progress
    by name
    select * from new_progress;

copy (
    select "sample", unnest(tool) as tool
    from identification_progress
) to '{{ output }}' (format csv);