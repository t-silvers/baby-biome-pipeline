create table reference_genomes (
    sample uinteger,
    reference_genome varchar,
    new_est_reads uinteger,
    fraction_total_reads float,
    primary key (sample)
);