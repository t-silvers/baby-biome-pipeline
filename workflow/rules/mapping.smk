checkpoint mapping_samplesheet:
    input:
        samplesheet='results/samplesheet.csv',
        identification='results/identification.csv',
    output:
        species_samplesheet='results/samplesheet-{species}.csv',
    localrule: True
    run:
        import pandas as pd

        samplesheet = pd.read_csv(input['samplesheet'])
        identification = pd.read_csv(input['identification'])
        
        species_mask = identification['species'] == wildcards.species

        (
            identification
            [species_mask]
            .merge(samplesheet, on='sample')
            .filter(['sample', 'fastq_1', 'fastq_2'])
            .to_csv(output['species_samplesheet'], index=False)
        )


rule bactmap:
    input:
        input='results/samplesheet-{species}.csv',
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
    output:
        'results/{species}/pipeline_info/pipeline_report.txt',
        'results/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
    params:
        pipeline='bactmap',
        profile='singularity',
        nxf='-work-dir results/{species}/work -config ' + bactmap_config,
        outdir='results/{species}',
    handover: True
    localrule: True
    envmodules:
        'apptainer/1.3.2',
        'nextflow/21.10',
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
        ancient('results/{species}/pipeline_info/pipeline_report.txt'),
        ancient('results/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml'),
        ancient('results/{species}/multiqc/multiqc_data/multiqc_fastp.yaml'),
        ancient('results/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml'),
    output:
        touch('results/{species}/fastp/{sample}_1.trim.fastq.gz'),
        touch('results/{species}/fastp/{sample}_2.trim.fastq.gz'),
        touch('results/{species}/samtools/{sample}.sorted.bam'),
        touch('results/{species}/variants/{sample}.vcf.gz'),
    localrule: True


def species_vcfs(wildcards):
    import pandas as pd

    species_samplesheet = (
        checkpoints.mapping_samplesheet
        .get(species=wildcards.species)
        .output['species_samplesheet']
    )

    samples = pd.read_csv(species_samplesheet)['sample']

    return expand(
        'results/{{species}}/variants/{sample}.vcf.gz',
        sample=samples
    )


rule:
    input:
        species_vcfs
    output:
        'results/data/variants/{species}.duckdb'
    params:
        glob="'results/{species}/variants/*.filtered.vcf.gz'"
    resources:
        cpus_per_task=24,
        mem_mb=120_000,
        runtime=10,
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export  MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
                FILTERED_VCFS={params.glob}

        duckdb -init config/duckdbrc-slurm {output} -c ".read workflow/scripts/models/annotated_vcfs.sql"
        '''