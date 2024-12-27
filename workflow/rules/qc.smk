rule total_reads:
    input:
        expand(
            workflow.source_path('results/{species}/multiqc/multiqc_data/multiqc_fastp.yaml'),
            species=all_species
        )
    output:
        (data_dir / 'data/qc/total_reads.duckdb').as_posix(),
    params:
        duckdbrc=duckdbrc_slurm,
        model=config['models']['qc']['fastp']
    resources:
        cpus_per_task=8,
        mem_mb=16_000
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1100))GB"

        yq -o=csv '.[] | [key, .summary.before_filtering.total_reads, .summary.after_filtering.total_reads]' {input} |\
          duckdb -init {params.duckdbrc} {output} -c ".read {params.model}"
        '''


rule samtools_stats:
    input:
        expand(
            workflow.source_path('results/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml'),
            species=all_species
        )
    output:
        (data_dir / 'data/qc/samtools_stats.duckdb').as_posix(),
    params:
        duckdbrc=duckdbrc_slurm,
        model=config['models']['qc']['samtools_stats']
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        yq -o=csv '.[] | [key, .reads_mapped_and_paired_percent, .reads_properly_paired_percent]' {input} |\
          duckdb -init {params.duckdbrc} {output} -c ".read {params.model}"
        '''


rule variant_quality_scores:
    input:
        expand(
            workflow.source_path('results/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml'),
            species=all_species
        )
    output:
        (data_dir / 'data/qc/variant_quality_scores.duckdb').as_posix(),
    params:
        duckdbrc=duckdbrc_slurm,
        model=config['models']['qc']['bcftools_stats_vqc']
    resources:
        cpus_per_task=8,
        mem_mb=16_000
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1100))GB"

        yq -o=csv 'to_entries[] | .key as $sample | .value | to_entries[] | [$sample, .key, .value]' {input} |\
          duckdb -init {params.duckdbrc} {output} -c ".read {params.model}"
        '''