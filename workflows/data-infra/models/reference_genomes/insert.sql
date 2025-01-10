insert or ignore into reference_genomes
    by name
    select * from read_csv('{{ input }}', hive_partitioning = false);

copy (
    select * from reference_genomes
    where new_est_reads > pow(10, cast('{{ read_pow }}' as int))
      and fraction_total_reads > cast('{{ read_frac }}' as float)
) to '{{ output }}' (format csv);