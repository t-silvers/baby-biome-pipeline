def mapping_sentinel(wildcards):
    """Get tool-specific paths to variant call results.

    This sentinel can be imperfect. Add check for individual files?
    """
    if wildcards.mapping_tool == 'bactmap':
        return ancient('results/bactmap/{species}/pipeline_info/pipeline_report.html')
    elif wildcards.mapping_tool.startswith('sarek_'):
        return ancient('results/sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml')
    elif wildcards.mapping_tool == 'snippy':
        return rules.all_snippy.input
    else:
        raise ValueError


def vcf_path_from_tool(wildcards):
    """Get tool-specific paths to variant call results."""
    if wildcards.mapping_tool == 'bactmap':
        return f'results/bactmap/{wildcards.species}/variants/{wildcards.sample}.filtered.vcf.gz'
    elif wildcards.mapping_tool == 'legacy_mapping':
        return f'results/legacy_mapping/{wildcards.species}/variants/{wildcards.sample}.calls.view.vcf.gz'
    elif wildcards.mapping_tool == 'snippy':
        return f'results/snippy/{wildcards.species}/variants/{wildcards.sample}.snps.raw.vcf'
    elif wildcards.mapping_tool == 'sarek_bcftools':
        return f'results/sarek/{wildcards.species}/variant_calling/bcftools/{wildcards.sample}/{wildcards.sample}.bcftools.vcf.gz'
    elif wildcards.mapping_tool == 'sarek_deepvariant':
        return f'results/sarek/{wildcards.species}/variant_calling/deepvariant/{wildcards.sample}/{wildcards.sample}.deepvariant.vcf.gz'
    elif wildcards.mapping_tool == 'sarek_freebayes':
        return f'results/sarek/{wildcards.species}/variant_calling/freebayes/{wildcards.sample}/{wildcards.sample}.freebayes.vcf.gz'
    elif wildcards.mapping_tool == 'sarek_haplotypecaller':
        return f'results/sarek/{wildcards.species}/variant_calling/haplotypecaller/{wildcards.sample}/{wildcards.sample}.haplotypecaller.vcf.gz'
    else:
        raise ValueError


rule vcf_to_parquet:
    input:
        mapping_sentinel,
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.parquet',
    params:
        vcf=lambda wildcards: vcf_path_from_tool(wildcards)
    envmodules:
        'vcf2parquet/0.4.1'
    shell:
        'vcf2parquet -i {params.vcf} convert -o {output}'


rule filter_variants:
    input:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.parquet',
    output:
        data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}.cleaned.parquet',
    params:
        # Files
        input=input[0],
        output=output[0],

        # Params
        ad=config['params']['variant_filter']['ad_ge'],
        ad_strand=config['params']['variant_filter']['ad_strand_ge'],
        alt_density_window=config['params']['variant_filter']['alt_density_window_half_size'],
        dp=config['params']['variant_filter']['dp_ge'],
        epp=config['params']['variant_filter']['epp_lt'],
        maf=config['params']['variant_filter']['maf_ge'],
        mq=config['params']['variant_filter']['mq_ge'],
        quality=config['params']['variant_filter']['quality_ge'],
        sp=config['params']['variant_filter']['sp_lt'],
    run:
        transform(models[f'clean_vcf_{wildcards.mapping_tool}'], params)