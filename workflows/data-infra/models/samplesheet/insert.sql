create temp table new_samples as
with
    prepared_samplesheet as (
        select * from read_csv('{{ input }}', hive_partitioning = false)
    ),

    -- For duplicates records on primary key, combine data into arrays
    final as (
        select library, id
            , array_agg(distinct family) as family
            , array_agg(distinct relationship) as relationship
            , array_agg(distinct donor) as donor
            , array_agg(distinct timepoint) as timepoint
            , array_agg(distinct collection_interval_category) as collection_interval_category
            , array_agg(distinct collection_interval) as collection_interval
            , array_agg(distinct species) as species
            , array_agg(distinct species_stdized) as species_stdized
            , array_agg(distinct plate) as plate
            , array_agg(distinct well) as well
            , array_agg(distinct barcode_1) as barcode_1
            , array_agg(distinct barcode_2) as barcode_2
            , array_agg(distinct notes) as notes
            , array_agg(distinct fastq_1) as fastq_1
            , array_agg(distinct fastq_2) as fastq_2
        from prepared_samplesheet
        group by library, id
    )
select * from final;

-- Add new samples to the samples table, ignoring existing records
insert or ignore into samples
    by name
    select nextval('sample') as "sample", *
    from new_samples;

copy (
    select "sample"
        , library
        , id
        , family[1] as family
        , relationship[1] as relationship
        , donor[1] as donor
        , collection_interval_category[1] as collection_interval_category
        , collection_interval[1] as collection_interval
        , species_stdized[1] as species_stdized
        , unnest(fastq_1) as fastq_1
        , unnest(fastq_2) as fastq_2
    from samples

    -- Filter out records (usually controls) with duplicate select fields on primary key
    where length(family) < 2
      and length(relationship) < 2
      and length(donor) < 2
      and length(collection_interval_category) < 2
      and length(collection_interval) < 2
      and length(species_stdized) < 2
) to '{{ output }}' (format csv);