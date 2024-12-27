create temp table dp_for_plots as

with selected_positions as (

    select
        -- distinct on (chromosome, position)
        chromosome, position
    
    from
        -- read_parquet(getenv('FILTERED_CALLS'))
        -- read_parquet('results/pseudogenomes/species=Bifidobacterium_bifidum/family=B002/filtered_calls.parquet')
        read_parquet('/viper/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_bifidum/positions.parquet')
        -- read_csv('/viper/ptmp/thosi/wv-multilib/results/trees/B001-Bifidobacterium_bifidum/msa-positions-varying.fas', delim='\t')
    
    where position in []

),

subset_positions as (

    select
        species
        , family
        , id
        , chromosome
        , position
        , reference
        , alternate
        , info_ADF
        , info_ADR

    from
        main

    where
        (chromosome, position) in (
            select
                (chromosome, position)
            from
                selected_positions
        )

),

cleaned as (

    select
        species
        , family
        , id
        , chromosome
        , position
        , case
            when alternate[1] is null then [reference]
            when alternate[1] is not null then list_concat([reference], alternate)
            else null
        end as allele
        , info_ADF
        , info_ADR
    
    from
        subset_positions

),

unnested as (

    select
        species
        , family
        , id
        , chromosome
        , position
        , unnest(allele) as allele
        , unnest(info_ADF) as forward
        , unnest(info_ADR) as reverse

    from
        cleaned
        
),

final as (

    select
        species
        , family
        , id
        , chromosome || ':' || position as coords
        , allele
        , forward
        , reverse

    from
        unnested

    where
        forward > 0
        or reverse > 0

    order by
        chromosome
        , position
        , family
        , id

)

select * from final;

copy dp_for_plots to '/dev/stdout' (format csv);

copy dp_for_plots to '/viper/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_bifidum/allele_depth_plot_data.csv' (format csv);
