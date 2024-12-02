rule total_reads:
    input:
        expand(
            'results/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
            species=isolated_species
        )
    output:
        'results/data/qc/total_reads.duckdb',
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
          duckdb -init config/duckdbrc-slurm {output} \
            -c ".read workflow/scripts/models/mqc_fastp.sql"
        '''


rule samtools_stats:
    input:
        expand(
            'results/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml',
            species=isolated_species
        )
    output:
        'results/data/qc/samtools_stats.duckdb',
    localrule: True
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        yq -o=csv '.[] | [key, .reads_mapped_and_paired_percent, .reads_properly_paired_percent]' {input} |\
          duckdb -init config/duckdbrc-local {output} \
            -c ".read workflow/scripts/models/mqc_samtools_stats.sql"
        '''


rule variant_quality_scores:
    input:
        expand(
            'results/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml',
            species=isolated_species
        )
    output:
        'results/data/qc/variant_quality_scores.duckdb',
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
          duckdb -init config/duckdbrc-slurm {output} \
            -c ".read workflow/scripts/models/mqc_bcftools_stats_vqc.sql"
        '''