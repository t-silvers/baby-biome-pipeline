rule bactmap_samplesheet:
    input:
        data_dir / 'data/samplesheets/main.csv',
        data_dir / 'data/identification/reference_genomes.csv',
    output:
        data_dir / 'data/samplesheets/bactmap_{species}.csv',
    log:
        log_dir / 'smk/mapping/bactmap_samplesheet_{species}.log'
    resources:
        njobs=50
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
        njobs=295
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


rule bactmap_vcf_to_parquet:
    input:
        data_dir / 'results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'
    output:
        data_dir / 'data/variants/species={species}/{partitions}/{sample}.filtered.vcf.parquet',
    resources:
        cpus_per_task=4,
        runtime=5,
        njobs=1
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {input} convert -o {output}'


rule bactmap_vcf_clean:
    input:
        data_dir / 'data/variants/species={species}/{partitions}/{sample}.filtered.vcf.parquet',
    output:
        data_dir / 'data/variants/species={species}/{partitions}/{sample}.cleaned.vcf.parquet',
    params:
        alt_density_window=config['mapping']['alt_density_window_half_size'],
        model=workflow.source_path(models['vcfs']['clean']),
    resources:
        cpus_per_task=4,
        mem_mb=8_000,
        runtime=15,
        njobs=1
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)


def aggregate_bactmap_vcfs(wildcards):
    import pandas as pd

    def covars_to_partitions(df, covars=None):
        """Form hive-partitioned path substring from df columns."""
        EXCLUDE = ['species', 'sample']

        if 'partitions' in df.columns:
            raise ValueError

        if covars is None:
            covars = [col for col in df.columns if col not in EXCLUDE]
        
        def row_to_partitions(row):
            return '/'.join([f'{k}={v}' for k, v in row.items()])

        return df.filter(covars).apply(row_to_partitions, axis=1)
        
    def cleaned_vcf_pq(df):
        path = 'data/variants/species={species}/{partitions}/{sample}.cleaned.vcf.parquet'.format(**df.to_dict())
        return (data_dir / path).as_posix()

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .rename(columns={'reference_genome': 'species'})
        .assign(partitions=lambda df: covars_to_partitions(df, covars=['family', 'id', 'library']))
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: cleaned_vcf_pq(df))
        .values
        .flatten()
    )


# PHONY
rule all_bactmap_vcfs:
    input:
        aggregate_bactmap_vcfs