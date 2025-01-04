create type tools as enum (
    select tool from read_csv('{{ maptools_seed }}')
);

create table mapping_progress (
    sample uinteger,
    tool tools[],
    primary key ("sample")
);