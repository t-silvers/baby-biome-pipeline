vcf_params = config['params']['variant_filter']


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
    elif wildcards.mapping_tool == 'legacy_mapping':
        raise NotImplementedError
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


checkpoint all_mapping:
    input:
        rules.all_vcfs.input
    output:
        'results/vcfs.csv',
    params:
        fields='tool|species|family|id|library|sample',
        glob=data / 'variants/tool=*/species=*/family=*/id=*/library=*/*.parquet',
        output=output[0],
        pat='variants/tool=(.*)/species=(.*)/family=(.*)/id=(.*)/library=(.*)/(.*).parquet',
    log:
        'logs/smk/all_mapping.log'
    localrule: True
    run:
        transform(models['vcfs']['manifest'], params, log=log[0])


def vcf_clean_model_from_tool(wildcards):
    if wildcards.mapping_tool in ['bactmap', 'sarek_bcftools']:
        return models['vcfs']['clean']['bcftools']
    elif wildcards.mapping_tool in ['sarek_freebayes', 'snippy']:
        return models['vcfs']['clean']['freebayes']
    elif wildcards.mapping_tool == 'sarek_haplotypecaller':
        return models['vcfs']['clean']['haplotypecaller']
    elif wildcards.mapping_tool in ['legacy_mapping', 'sarek_deepvariant']:
        raise NotImplementedError
    else:
        raise ValueError


def vcf_pq(wildcards):
    import pandas as pd
    vcfs = pd.read_csv(
        checkpoints.all_mapping.get(**wildcards).output[0]
    ).astype(str)
    vcf_row = vcfs[
        (vcfs['tool'] == str(wildcards.mapping_tool)) &
        (vcfs['species'] == str(wildcards.species)) &
        (vcfs['family'] == str(wildcards.family)) &
        (vcfs['id'] == str(wildcards.id)) &
        (vcfs['library'] == str(wildcards.library)) &
        (vcfs['sample'] == str(wildcards.sample))
    ]
    if vcf_row.shape[0] != 1:
        print(vcf_row)
        raise ValueError
    return vcf_row['file'].values[0]


rule filter_variants:
    input:
        vcf_pq
    output:
        multiext(
            (data / 'variants/tool={mapping_tool}/species={species}/family={family}/id={id}/library={library}/{sample}').as_posix(),
            '.clean.parquet', '.filtered.parquet'
        )
    params:
        # Files
        input=input[0],
        cleaned=output[0],
        filtered=output[1],

        # Params
        ad=vcf_params['ad_ge'],
        ad_strand=vcf_params['ad_strand_ge'],
        alt_density_window=vcf_params['alt_density_window_half_size'],
        dp=vcf_params['dp_ge'],
        epp=vcf_params['epp_lt'],
        maf=vcf_params['maf_ge'],
        mq=vcf_params['mq_ge'],
        quality=vcf_params['quality_ge'],
        sp=vcf_params['sp_lt'],
    run:
        for model in [
            vcf_clean_model_from_tool(wildcards),
            models['vcfs']['filter'],
        ]:
            transform(model, params)


def aggregate_filtered_variants(wildcards):
    import pandas as pd

    path_template = (
        data /
        'variants/tool={mapping_tool}/species={species}/family={family}' /
        'id={id}/library={library}/{sample}.filtered.parquet'
    ).as_posix()

    return (
        pd.read_csv(checkpoints.all_mapping.get(**wildcards).output[0])
        .filter(['tool', 'species', 'family', 'id', 'library', 'sample'])
        .rename(columns={'tool': 'mapping_tool'})
        .drop_duplicates()
        .transpose()
        .apply(lambda df: path_template.format(**df.to_dict()))
        .values
        .flatten()
    )


# def variants_db_path(wildcards):
#     return (
#         data / 'variants' /
#         f'species={wildcards.species}' /
#         f'family={wildcards.family}' /
#         f'tool={wildcards.mapping_tool}' /
#         'filtered.duckdb'
#     )


# rule create_variants_db:
#     output:
#         data / '.variants_{species}_{family}_{mapping_tool}_duckdb',
#     params:
#         db=lambda wildcards: variants_db_path(wildcards),
#         contigs=seeds['types']['contigs'],
#         species=lambda wildcards: wildcards.species,
#     log:
#         'logs/smk/create_variants_db_{species}_{family}_{mapping_tool}.log'
#     run:
#         transform(
#             models['variants']['create_db'],
#             params, db=params['db'], log=log[0]
#         )
#         shell('touch ' + output[0])


# rule update_variants_db:
#     input:
#         aggregate_filtered_variants,
#         data / '.variants_{species}_{family}_{mapping_tool}_duckdb',
#     output:
#         data / 'variants/species={species}/family={family}/tool={mapping_tool}/snvs-filtered.parquet',
#     params:
#         # Files
#         db=lambda wildcards: variants_db_path(wildcards),
#         max_temp_directory_size='300GB',
#         output=output[0],
#         vcf_pq_glob=lambda wildcards: data / 'variants' / f'tool={wildcards.mapping_tool}' / f'species={wildcards.species}' / f'family={wildcards.family}' / '**/*.filtered.parquet',

#         # Params
#         frac_cohort_core=config['params']['variant_filter']['frac_cohort_core'],
#     log:
#         'logs/smk/update_variants_db_{species}_{family}_{mapping_tool}.log'
#     run:
#         transform(
#             models['variants']['insert_db'],
#             params, db=params['db'], log=log[0]
#         )


# def aggregate_filtered_snvs(wildcards):
#     import pandas as pd

#     path_template = (
#         # TODO: temp
#         data / 'variants/species={species}/family={family}/tool=bactmap/snvs-filtered.parquet'
#         # data / 'variants/species={species}/family={family}/tool={mapping_tool}/snvs-filtered.parquet'
#     ).as_posix()

#     return (
#         pd.read_csv(checkpoints.all_mapping.get(**wildcards).output[0])
#         .filter(['species', 'family', 'tool'])
#         .rename(columns={'tool': 'mapping_tool'})
#         .drop_duplicates()
#         .transpose()
#         .apply(lambda df: path_template.format(**df.to_dict()))
#         .values
#         .flatten()
#     )