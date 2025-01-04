wildcard_constraints:
    family=r'B\d+',
    library=config['wildcards']['libraries'],
    mapping_tool='bactmap|legacy_mapping|snippy|sarek_bcftools|sarek_deepvariant|sarek_freebayes|sarek_haplotypecaller',
    sample=r'\d+',
    species=config['wildcards']['species']


include: '../mapping-modules/bactmap.smk'
include: '../mapping-modules/legacy_mapping.smk'
include: '../mapping-modules/sarek.smk'
include: '../mapping-modules/snippy.smk'


def vcf_template_from_tool(wildcards):

    # TODO: Replace with a "standardizing" bcftools annotate command within each module

    VCF_TEMPLATE = (results / '{pipeline}/{{species}}/variants/{{sample}}.{ext}').as_posix()

    if wildcards.mapping_tool == 'bactmap':
        return VCF_TEMPLATE.format(pipeline='bactmap', ext='filtered.vcf.gz')
    elif wildcards.mapping_tool == 'legacy_mapping':
        return VCF_TEMPLATE.format(pipeline='legacy_mapping', ext='calls.view.vcf.gz')
    elif wildcards.mapping_tool == 'snippy':
        return VCF_TEMPLATE.format(pipeline='snippy', ext='snps.raw.vcf')
    elif wildcards.mapping_tool == 'sarek_bcftools':
        return VCF_TEMPLATE.format(pipeline='sarek', ext='bcftools.vcf.gz')
    elif wildcards.mapping_tool == 'sarek_deepvariant':
        return VCF_TEMPLATE.format(pipeline='sarek', ext='deepvariant.vcf.gz')
    elif wildcards.mapping_tool == 'sarek_freebayes':
        return VCF_TEMPLATE.format(pipeline='sarek', ext='freebayes.vcf.gz')
    elif wildcards.mapping_tool == 'sarek_haplotypecaller':
        return VCF_TEMPLATE.format(pipeline='sarek', ext='haplotypecaller.vcf.gz')
    else:
        raise ValueError

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: Adapt this logic to non-multi-library data (and non-module?) use case

rule mapping_db:
    output:
        data / 'mapping.duckdb',
    params:
        model=workflow.source_path(models['workflow']['mapping']['create']),
        maptools_seed=workflow.source_path(seeds['types']['maptools']),
    log:
        logdir / 'smk/mapping/mapping_db.log'
    resources:
        njobs=1,
    run:
        transform(params['model'], params, db=output[0], log=log[0])


rule update_mapping_db:
    input:
        ancient(data / 'mapping.duckdb'),
    output:
        results / 'samplesheets/mapping_progress.csv',
    params:
        model=workflow.source_path(models['workflow']['mapping']['insert']),
    log:
        logdir / 'smk/mapping/update_identification_db.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        runtime=5,
        njobs=1
    run:        
        params.update({'output': output[0]})
        transform(params['model'], params, db=input[0], log=log[0])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


rule vcf_to_parquet:
    input:
        vcf_template_from_tool
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.parquet',
    resources:
        cpus_per_task=4,
        runtime=5,
        njobs=1,
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {input} convert -o {output}'


rule clean_vcf:
    input:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.parquet',
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.cleaned.parquet',
    params:
        alt_density_window=config['mapping']['alt_density_window_half_size'],
        model=lambda wildcards: workflow.source_path(models['vcfs'][wildcards.mapping_tool]),
    resources:
        cpus_per_task=4,
        mem_mb=8_000,
        runtime=15,
        njobs=1,
    run:
        params.update({'input': input[0], 'output': output[0]})
        transform(params['model'], params)


def aggregate_vcfs(wildcards):
    import pandas as pd

    vcfs = [
        f'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}'
        for mapping_tool in config['tools']['mapping'].split('|')
    ]

    def cleaned_vcf_pq(df):
        return '|'.join(list(map(lambda x: data_path_from_template(x + '.cleaned.parquet', df.to_dict()), vcfs)))

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .rename(columns={'reference_genome': 'species'})
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: cleaned_vcf_pq(df))
        .str.split('|')
        .explode()
        .values
        .flatten()
    )


# PHONY
rule all_mapping:
    input:
        aggregate_vcfs
