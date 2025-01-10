create type relationship as enum (
    select relationship from read_csv('{{ relationship_seed }}')
);

create type species as enum (
    select species from read_csv('{{ species_seed }}')
);

create type timepoint as enum (
    select timepoint from read_csv('{{ timepoint_seed }}')
);

create table samples (
    sample uinteger,
    library varchar,
    id varchar,
    family varchar[],
    relationship relationship[],
    donor varchar[],
    timepoint varchar[],
    collection_interval_category timepoint[],
    collection_interval interval[],
    species varchar[],
    species_stdized species[],
    plate usmallint[],
    well varchar[],
    barcode_1 varchar[],
    barcode_2 varchar[],
    notes varchar[],
    fastq_1 varchar[],
    fastq_2 varchar[],
    primary key (library, id)
);

create sequence "sample";