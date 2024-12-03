checkpoint bactmap_samplesheet:
    input:
        'results/samplesheet.csv',
        'results/identification.csv',
    output:
        'results/samplesheets/bactmap_{species}.csv',
    localrule: True
    run:
        import pandas as pd

        samplesheet = pd.read_csv(input[0])
        identification = pd.read_csv(input[1])
        
        species_mask = identification['species'] == wildcards.species

        (
            identification
            [species_mask]
            .merge(samplesheet, on='sample')
            .filter(['sample', 'fastq_1', 'fastq_2'])
            .to_csv(output[0], index=False)
        )


rule bactmap:
    input:
        input='results/samplesheets/bactmap_{species}.csv',
    output:
        'results/bactmap/{species}/pipeline_info/pipeline_report.txt',
        'results/bactmap/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
    params:
        pipeline='bactmap',
        profile='singularity',
        nxf='-work-dir results/bactmap/{species}/work -config ' + bactmap_config,
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
        outdir='results/bactmap/{species}',
    handover: True
    localrule: True
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.04.4',
        'jdk/17.0.10'
    container:
        'docker://nfcore/bactmap'
    wrapper:
        'https://raw.githubusercontent.com/fm-key-lab/snakemake-wrappers/nf-core/bio/nf-core'


rule:
    """Collect mapping output.
    
    Collect mapping output such that the sample wildcard can be
    resolved by downstream rules.
    """
    input:
        ancient('results/bactmap/{species}/pipeline_info/pipeline_report.txt'),
        ancient('results/bactmap/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml'),
        ancient('results/bactmap/{species}/multiqc/multiqc_data/multiqc_fastp.yaml'),
        ancient('results/bactmap/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml'),
    output:
        touch('results/bactmap/{species}/fastp/{sample}_1.trim.fastq.gz'),
        touch('results/bactmap/{species}/fastp/{sample}_2.trim.fastq.gz'),
        touch('results/bactmap/{species}/samtools/{sample}.sorted.bam'),
        touch('results/bactmap/{species}/variants/{sample}.vcf.gz'),
        touch('results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'),
    localrule: True


rule:
    input:
        'results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'
    output:
        'results/variants/species={species}/family={family}/sample={sample}/annot_vcf.parquet'
    resources:
        cpus_per_task=4,
        runtime=5
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {input} convert -o {output}'



def species_family_vcfs(wildcards):
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[1]
    )

    mapping_samplesheet = pd.read_csv(
        checkpoints.bactmap_samplesheet
        .get(species=wildcards.species)
        .output[0]
    )

    samples = (
        samplesheet
        .filter(['sample', 'family'])
        .query(
            f"family == {wildcards['family']}"
        )
        .merge(
            mapping_samplesheet,
            how='right'
        )
        ['sample']
    )

    return expand(
        'results/variants/species={{species}}/family={{family}}/sample={sample}/annot_vcf.parquet',
        sample=samples
    )


rule:
    input:
        species_family_vcfs
    output:
        'results/variants/species={species}/family={family}/variants.duckdb'
    params:
        glob="'results/variants/species={species}/family={family}/sample=*/annot_vcf.parquet'"
    resources:
        cpus_per_task=24,
        mem_mb=120_000,
        runtime=10,
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export VCFS={params.glob};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output} -c ".read ' + 
            workflow.source_path('../scripts/models/annotated_vcf_parquet.sql') + 
            '"'
        )