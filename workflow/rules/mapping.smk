localrules: bactmap_samplesheet, bactmap, bactmap_vcf


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
        nxf=config['mapping']['bactmap']['nxf_args'] + ' -work-dir logs/nf/taxprofiler/{species}/work',
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


rule bactmap_output:
    input:
        expand(
            'results/bactmap/{species}/pipeline_info/pipeline_report.html',
            species=config['wildcards']['species'].split('|')
        )


rule bactmap_vcf:
    """Collect mapping output.
    
    Collect mapping output such that the sample wildcard can be
    resolved by downstream rules.
    """
    input:
        'results/bactmap/{species}/pipeline_info/pipeline_report.html',
    output:
        touch('results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'),
    resources:
        njobs=1


rule vcf_to_parquet:
    input:
        'results/bactmap/{species}/variants/{sample}.filtered.vcf.gz'
    output:
        'results/tmp_data/species={species}/sample={sample}/filtered_vcf.parquet'
    resources:
        cpus_per_task=2,
        runtime=5,
        njobs=1
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {input} convert -o {output}'


def collect_vcfs(wildcards):
    import pandas as pd

    mapping_samplesheet = pd.read_csv(
        checkpoints.bactmap_samplesheet
        .get(species=wildcards.species)
        .output[0]
    )

    samples = mapping_samplesheet['sample']

    return expand(
        'results/tmp_data/species={{species}}/sample={sample}/filtered_vcf.parquet',
        sample=samples
    )


rule variants_db:
    input:
        collect_vcfs
    output:
        'results/variants/{species}.duckdb'
    params:
        glob_pq="'results/tmp_data/species={species}/sample=*/filtered_vcf.parquet'",
        glob_vcf="'results/bactmap/{species}/variants/*.filtered.vcf.gz'",
    resources:
        cpus_per_task=8,
        mem_mb=32_000,
        runtime=15,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        (
            'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
            'export VCFS_PQ={params.glob_pq};' +
            'export VCFS={params.glob_vcf};' +
            'duckdb -init ' + 
            workflow.source_path('../../config/duckdbrc-slurm') +
            ' {output} -c ".read ' + 
            workflow.source_path('../scripts/models/annotated_vcf_parquet.sql') + 
            '"'
        )