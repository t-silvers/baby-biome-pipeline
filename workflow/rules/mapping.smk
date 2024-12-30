wildcard_constraints:
    ext='filtered.vcf.gz|calls.view.vcf.gz|snps.vcf',
    mapping_tool=config['tools']['mapping']


include: '../mapping-modules/bactmap.smk'
include: '../mapping-modules/legacy_mapping.smk'
include: '../mapping-modules/sarek.smk'
include: '../mapping-modules/snippy.smk'


rule vcf_to_parquet:
    input:
        results / '{mapping_tool}/{species}/variants/{sample}.{ext}'
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.{ext}.parquet',
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
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.{ext}.parquet',
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.{ext}.cleaned.parquet',
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

    # TODO: Standardize vcf outputs?
    VCF_TEMPLATE = 'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.{ext}'
    VCF_EXTS = {
        'bactmap': 'filtered.vcf.gz',
        'legacy_mapping': 'calls.view.vcf.gz',
        'snippy': 'snps.vcf',
    }

    '''
    results / 'bactmap/             {species}/  variants/                           {sample}.filtered.vcf.gz'

    results / 'legacy_mapping/      {species}/  variants/                           {sample}.calls.view.vcf.gz'
    
    results / 'sarek/               {species}/  variantcalling/{vc_tool}/{sample}/  {sample}.{vc_tool}.vcf.gz'
    
    results / 'snippy/              {species}/  variants/                           {sample}.snps.vcf'





    VCF_TEMPLATE = 'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.{ext}'
    VCF_EXTS = {
        'bactmap': 'filtered.vcf.gz',
        'legacy_mapping': 'calls.view.vcf.gz',
        'snippy': 'snps.vcf',
    }


    results / '{mapping_tool}/{species}/variants/{sample}.{ext}'
    
    SAREK_VARIANT_CALLING_TOOLS = ['bcftools', 'deepvariant', 'freebayes', 'haplotypecaller']

    results / 'sarek/{{species}}/variantcalling/{vc_tool}/{{sample}}/{{sample}}.{vc_tool}.vcf.gz',
    vc_tool=[
        tool for tool in config['mapping']['sarek']['tools'].split(',') 
        if tool in SAREK_VARIANT_CALLING_TOOLS
    ]
    '''

    vcfs = [
        VCF_TEMPLATE.format(mapping_tool=tool, ext=VCF_EXTS[tool]) 
        for tool in config['tools']['mapping'].split('|')
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
