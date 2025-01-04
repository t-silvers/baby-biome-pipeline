create type vc_types as enum ('snp', 'indel');

create type contigs as enum (
    select contig 
    from read_csv('{{ contigs_seed }}')
    where species = '{{ species }}'
);

create table variants (
    -- Required
    contig contigs,
    start_pos uinteger,
    sample uinteger,
    alleles varchar[],
    
    -- Extra
    end_pos uinteger,
    qual decimal(6, 1),
    
    -- Extra: info
    info_AD usmallint[],
    info_ADF usmallint[],
    info_ADR usmallint[],
    info_DP usmallint,
    info_MQ decimal(4, 1),
    info_TYPE vc_types[],

    -- Extra: format
    format_GT varchar,
    format_PL float[],
    format_SP usmallint,

    primary key (contig, start_pos, "sample", end_pos)
);