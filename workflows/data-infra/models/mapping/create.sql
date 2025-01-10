create type variant_type as enum ('snp', 'indel');

create type contigs as enum (
    select contig
    from read_csv('{{ contigs_seed }}')
    where species = '{{ species }}'
);

create table variants (
    -- ~*~ Required ~*~
    contig contigs,
    start_pos uinteger,
    sample uinteger,
    alleles varchar[],
    -- ~*~ Extra ~*~
    -- For precise variants, END is POS + length of REF allele - 1, and the for imprecise variants the corresponding best estimate.
    end_pos uinteger,
    --
    qual decimal(6, 1),
    -- ~*~ Extra: info ~*~
    --
    info_AD usmallint[],
    --
    info_ADF usmallint[],
    --
    info_ADR usmallint[],
    -- Depth of coverage
    info_DP usmallint,
    -- end-placement probability score
    info_EPP decimal(4, 1)[],
    -- Mapping quality
    info_MQ decimal(4, 1),
    -- mean mapping quality of observed reference, alternate alleles
    info_MQM decimal(4, 1)[],
    -- strand-bias probability score
    info_SP decimal(4, 1)[],
    -- variant type
    info_TYPE variant_type[],
    -- ~*~ Extra: format ~*~
    --
    format_GT varchar,
    -- phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field) (Integers)
    format_PL float[],
    --
    format_SP decimal(4, 1),
    primary key (contig, start_pos, "sample", end_pos)
);