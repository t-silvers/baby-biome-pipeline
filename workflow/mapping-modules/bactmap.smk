rule bactmap_samplesheet:
    input:
        data_dir / 'data/samplesheet.csv',
        data_dir / 'data/identification/reference_genomes.csv',
    output:
        data_dir / 'data/samplesheets/bactmap_{species}.csv',
    log:
        log_dir / 'smk/mapping/bactmap_samplesheet_{species}.log'
    resources:
        njobs=1,
    run:
        import pandas as pd

        BACTMAP_COLS = ['sample', 'fastq_1', 'fastq_2']

        samplesheet = pd.read_csv(input[0])
        identification = pd.read_csv(input[1])
        
        species_mask = identification['reference_genome'] == wildcards.species

        (
            identification
            [species_mask]
            .merge(samplesheet, on='sample')
            .filter(BACTMAP_COLS)
            .to_csv(output[0], index=False)
        )


rule bactmap:
    input:
        data_dir / 'data/samplesheets/bactmap_{species}.csv',
    output:
        data_dir / 'results/bactmap/{species}/pipeline_info/pipeline_report.html',
        data_dir / 'results/bactmap/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
        data_dir / 'results/bactmap/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml',
        data_dir / 'results/bactmap/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml',
    params:
        # Dirs
        outdir=lambda wildcards: data_dir / f'results/bactmap/{wildcards.species}',
        workdir=lambda wildcards: log_dir / f'nxf/bactmap/{wildcards.species}/work',
        
        # Generic params
        config=config['mapping']['bactmap']['config'],
        profile=config['mapping']['bactmap']['profiles'],
        
        # Pipeline params
        extra=config['mapping']['bactmap']['extra'],
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
    log:
        log_dir / 'smk/mapping/bactmap_{species}.log'
    handover: True
    resources:
        njobs=295,
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run nf-core/bactmap \
          -config {params.config} \
          -profile {params.profile} \
          -resume \
          -work-dir {params.workdir} \
          --input {input} \
          --outdir {params.outdir} \
          --reference {params.reference} \
          {params.extra}
        '''


def aggregate_bactmap(wildcards):
    import pandas as pd

    species = (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .query('reference_genome in @wc_params["species"]')
        ['reference_genome']
        .unique()
    )

    return expand(
        data_dir / 'results/bactmap/{species}/pipeline_info/pipeline_report.html',
        species=species
    )


# PHONY
rule all_bactmap:
    input:
        aggregate_bactmap


rule bactmap_vcf:
    """Collect bcftools output from bactmap."""
    input:
        data_dir / 'results/bactmap/{species}/pipeline_info/pipeline_report.html',
    output:
        touch(data_dir / 'results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'),
    resources:
        njobs=1
    localrule: True