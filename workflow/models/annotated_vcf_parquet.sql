set preserve_insertion_order = false;

create table annotated_vcfs as

with family_sample_info as (

    select
        family
        , relationship
        , donor
        , id
        , timepoint
        , timepoint_day
        , "sample"
    
    from
        read_csv(
            'results/samplesheet.csv'
        )

), 

raw_vcfs as (
    
    select
        "sample"
        , species
        , chromosome
        , position
        , reference
        , alternate
        , quality
        , "filter"
        , info_ADF
        , info_ADR
        , info_AD
        , info_INDEL
    
    from
        read_parquet(
            getenv('VCFS_PQ'),
            hive_partitioning = true,
            hive_types = {
                'species': varchar,
                'sample': uinteger
            }
        )
    
    where
        "sample" in (

            select "sample" from family_sample_info

        )

), 

cleaned_vcfs as (

    select
        * exclude(alternate, info_ADF, info_ADR)
        , array_transform(
            alternate, x -> nullif(x, '')
        ) as alternate
        , array_transform(
            info_ADF, x -> cast(x as usmallint)
        ) as info_ADF
        , array_transform(
            info_ADR, x -> cast(x as usmallint)
        ) as info_ADR

    from
        raw_vcfs

), 

gap_metrics as (

    select 
        "sample"
        , chromosome
        , position
        , case
            when info_INDEL then null
            else nearest_snv
        end as info_SNPGAP
        , case
            when info_INDEL then null
            else indel_within_3
        end as info_INDELGAP_3
        , case
            when info_INDEL then null
            else indel_within_10
        end as indel_INDELGAP_10

    from
        (

            select
                *
                , position - lag(position) over (
                    
                    partition by
                        "sample", chromosome, info_INDEL
                    
                    order by
                        position
                
                ) as lag_dist
                , lead(position) over (
                    
                    partition by
                        "sample", chromosome, info_INDEL
                    
                    order by
                        position
                
                ) - position as lead_dist
                , case
                    when lag_dist is null then lead_dist
                    when lead_dist is null then lag_dist
                    when lag_dist > lead_dist then lead_dist
                    when lead_dist >= lag_dist then lag_dist
                    else null
                end as nearest_snv
                , max(info_INDEL) over (
                    partition by "sample", chromosome 
                    order by position
                    range between 3. preceding 
                        and 3. following
                )
                as indel_within_3
                , max(info_INDEL) over (
                    partition by "sample", chromosome 
                    order by position
                    range between 10. preceding 
                        and 10. following
                )
                as indel_within_10

            from
                cleaned_vcfs

            where
                alternate[1] is not null

        )

-- Add strand bias metric

),

strand_bias_metric as (
    
    select
        cast(
            regexp_extract(
                "filename",
                'variants/(\d+).filtered',
                1
            ) as uinteger
        ) as "sample"
        , CHROM as chromosome
        , POS as position
        , try_cast(
            list_extract(
                string_split(FORMAT_values, ':'),
                list_position(
                    string_split(FORMAT, ':'),
                    'SP'
                )
            ) as usmallint
        ) as format_SP
    
    from
        read_csv(
            getenv('VCFS'),
            auto_detect = False,
            columns = {
                'CHROM': 'varchar',
                'POS': 'uinteger',
                'ID': 'varchar',
                'REF': 'varchar',
                'ALT': 'varchar',
                'QUAL': 'decimal(4, 1)',
                'FILTER': 'varchar',
                'INFO': 'varchar',
                'FORMAT': 'varchar',
                'FORMAT_values': 'varchar',
            },
            compression = 'gzip',
            delim = '\t',
            filename = True,
            header = True,
            ignore_errors = True,
            new_line = '\n', 
            nullstr = '.'
        )

), 

vcfs_with_metrics as (

    select
        v.*
        , gap.info_SNPGAP
        , gap.info_INDELGAP_3
        , gap.indel_INDELGAP_10
        , sp.format_SP
    
    from
        cleaned_vcfs v
    
    left join
        gap_metrics gap
    
    on
        v.sample = gap.sample
        and v.chromosome = gap.chromosome
        and v.position = gap.position

    left join
        strand_bias_metric sp
    
    on
        v.sample = sp.sample
        and v.chromosome = sp.chromosome
        and v.position = sp.position

), 

combined as (

    select
        s.* exclude("sample")
        , v.* exclude("sample")
    
    from
        family_sample_info s
        , vcfs_with_metrics v
    
    where
        s.sample = v.sample

), 

final as (

    select
        species
        , family
        , relationship
        , donor
        , id
        , timepoint
        , timepoint_day
        , chromosome
        , position
        , reference
        , alternate
        , quality
        , "filter"
        , info_INDEL
        , info_AD
        , info_ADF
        , info_ADR
        , info_SNPGAP
        , info_INDELGAP_3
        , indel_INDELGAP_10
        , format_SP

    from
        combined

    order by
        species
        , id
        , chromosome
        , position

)

select * from final;

copy
    annotated_vcfs
    
to 'results/variants' (
        format parquet,
        partition_by (species, family, relationship, donor, id, chromosome),
        overwrite_or_ignore 1
);