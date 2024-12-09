checkpoint bactmap_samplesheet:
    input:
        'resources/samplesheets/main.csv',
        'resources/reference_genomes.csv',
    output:
        'resources/samplesheets/bactmap_{species}.csv',
    log:
        'logs/smk/mapping/bactmap_samplesheet_{species}.log'
    resources:
        njobs=50
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
        input='resources/samplesheets/bactmap_{species}.csv',
    output:
        'results/bactmap/{species}/pipeline_info/pipeline_report.html',
    params:
        nxf=config['mapping']['bactmap']['nxf_args'] + ' -work-dir logs/nxf/taxprofiler/{species}/work',
        nxf_log='logs/nxf/bactmap/{species}.log',
        outdir='results/bactmap/{species}',
        pipeline='bactmap',
        profile=config['mapping']['bactmap']['profiles'],
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species]
    log:
        'logs/smk/mapping/bactmap_{species}.log'
    handover: True
    resources:
        njobs=200
    envmodules:
        'apptainer/1.3.2',
        'nextflow/24.04.4',
        'jdk/17.0.10'
    container:
        'docker://nfcore/bactmap'
    wrapper:
        'https://raw.githubusercontent.com/fm-key-lab/snakemake-wrappers/nf-core/bio/nf-core'


def agg_bactmap(wildcards):
    import pandas as pd

    references = pd.read_csv(
        checkpoints.reference_identification
        .get(**wildcards)
        .output[0]
    )
    
    species = references['species'].unique()

    return expand(
        'results/bactmap/{species}/pipeline_info/pipeline_report.html',
        species=species
    )


rule bactmap_reports:
    input:
        agg_bactmap


rule bactmap_vcf:
    """Collect mapping output.
    
    Collect mapping output such that the sample wildcard can be
    resolved by downstream rules.
    """
    input:
        ancient('results/bactmap/{species}/pipeline_info/pipeline_report.html'),
    output:
        touch('results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'),
    resources:
        njobs=1


rule vcf_to_parquet:
    input:
        'results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'
    output:
        'data/variants/species={species}/family={family}/id={id}/library={library}/{sample}.filtered.vcf.parquet',
    resources:
        cpus_per_task=4,
        runtime=5,
        njobs=1
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {input} convert -o {output}'


rule vcf_clean:
    input:
        'data/variants/species={species}/family={family}/id={id}/library={library}/{sample}.filtered.vcf.parquet'
    output:
        'data/variants/species={species}/family={family}/id={id}/library={library}/{sample}.cleaned.vcf.parquet'
    params:
        model=config['models']['variants']['vcf_clean']
    resources:
        cpus_per_task=4,
        mem_mb=8_000,
        runtime=15,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export VCF_PQ={input};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' -c ".read {params.model}" > {output}'
        )


def aggregate_vcfs(wildcards):
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[1]
    )

    mapping = pd.read_csv(
        checkpoints
        .bactmap_samplesheet
        .get(species=wildcards.species)
        .output[0]
    )

    return (
        samplesheet
        [samplesheet['family'] == wildcards['family']]
        .filter(['sample', 'id', 'library'])
        .merge(mapping)
        .transpose()
        .apply(lambda df: 'data/variants/species={{species}}/family={{family}}/id={id}/library={library}/{sample}.cleaned.vcf.parquet'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule snvs_db:
    input:
        aggregate_vcfs
    output:
        'data/variants/species={species}/family={family}/snvs.duckdb'
    params:
        glob="'data/variants/species={species}/family={family}/id=*/library=*/*.cleaned.vcf.parquet'",
        model=config['models']['variants']['snvs']
    log:
        'logs/smk/mapping/snvs_db_{species}_{family}.log'
    resources:
        cpus_per_task=32,
        mem_mb=240_000,
        runtime=10,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export VCFS={params.glob};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output} -c ".read {params.model}"'
        )


def aggregate_snvs_db(wildcards):
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[1]
    )

    references = pd.read_csv(
        checkpoints.reference_identification
        .get(**wildcards)
        .output[0]
    )

    mapping = pd.concat([
        pd.read_csv(checkpoints.bactmap_samplesheet.get(species=spp).output[0])
        for spp in references['species'].unique() 
        if spp in config['wildcards']['species'].split('|')
    ])

    return (
        samplesheet
        [samplesheet['family'].isin(config['wildcards']['families'].split('|'))]
        .filter(['sample', 'family'])
        .merge(mapping)
        .merge(references)
        .filter(['species', 'family'])
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: 'data/variants/species={species}/family={family}/snvs.duckdb'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule all_snvs:
    input:
        aggregate_snvs_db
