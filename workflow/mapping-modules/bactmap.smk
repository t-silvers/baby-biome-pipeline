import pathlib


rule bactmap_samplesheet:
    input:
        results / 'samplesheets/samplesheet.csv',
        results / 'samplesheets/reference_genomes.csv',
    output:
        results / 'samplesheets/bactmap_{species}.csv',
    log:
        logdir / 'smk/mapping/bactmap_samplesheet_{species}.log'
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

        # TODO: Remove samples that have already been analyzed; otherwise,
        #       relies on nxf dependency tracking. Use run db (sample,tool,...).


rule bactmap:
    input:
        results / 'samplesheets/bactmap_{species}.csv',
    output:
        results / 'bactmap/{species}/pipeline_info/pipeline_report.html',
        results / 'bactmap/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
        results / 'bactmap/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml',
        results / 'bactmap/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml',
    params:
        # Dirs
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent.parent,
        workdir=lambda wildcards: logdir / f'nxf/bactmap/{wildcards.species}/work',
        
        # Generic params
        config=config['mapping']['bactmap']['config'],
        profile=config['mapping']['bactmap']['profiles'],
        version=config['mapping']['bactmap']['version'],

        # Pipeline params
        extra=config['mapping']['bactmap']['extra'],
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
    log:
        logdir / 'smk/mapping/bactmap_{species}.log'
    handover: True
    resources:
        njobs=max_submit,
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
        results / 'bactmap/{species}/pipeline_info/pipeline_report.html',
        species=species
    )


# PHONY
rule all_bactmap:
    input:
        aggregate_bactmap


rule bactmap_vcf:
    """Collect bcftools output from bactmap."""
    input:
        results / 'bactmap/{species}/pipeline_info/pipeline_report.html',
    output:
        touch(results / 'bactmap/{species}/variants/{sample}.filtered.vcf.gz'),
    resources:
        njobs=1
    localrule: True