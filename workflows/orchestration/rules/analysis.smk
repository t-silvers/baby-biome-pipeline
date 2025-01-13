module widevariant:
    snakefile: '../../analysis/Snakefile'
    config: config


use rule * from widevariant


localrules:
    bactmap,
    bactmap_samplesheet,
    sarek,
    sarek_samplesheet,
    snippy_samplesheet,
    srst2_samplesheet,
    taxprofiler,
    taxprofiler_samplesheet,


checkpoint taxprofiler_samplesheet:
    input:
        'results/.identification_cache',
        sample_info='resources/samplesheet.csv',
    output:
        samplesheet='resources/samplesheets/taxprofiler.csv',
    params:
        # Files
        db=data / 'identification_cache.duckdb',
        input=input['sample_info'],
        output=output['samplesheet'],
        
        # Params
        # NOTE: Must conform to EBI ENA controlled vocabulary
        instrument_platform='ILLUMINA',
    log:
        'logs/smk/taxprofiler_samplesheet.log'
    run:
        transform(
            models['nfcore_inputs']['taxprofiler'],
            params, db=params['db'], log=log[0], readonly=True
        )


rule reference_genomes_from_bracken:
    input:
        data / '.reference_genomes_duckdb',
        rules.all_taxprofiler.input,
        sample_info='resources/samplesheet.csv',
    output:
        samplesheet='resources/samplesheet-with-ref.csv',
    params:
        # Files
        db=data / 'reference_genomes.duckdb',
        bracken_glob=data / 'identification/tool=taxprofiler/db=*/family=*/id=*/library=*/*.bracken.parquet',
        sample_info=input['sample_info'],
        samplesheet=output['samplesheet'],

        # Params
        read_frac=config['params']['reference_genome']['minimum_fraction_reference'],
        read_pow=config['params']['reference_genome']['minimum_readspow_reference'],
        available_genomes=config['wildcards']['genomes'],
    log:
        'logs/smk/reference_genomes_from_bracken.log'
    run:
        transform(models['reference_genomes']['update_db'], params, db=params['db'], log=log[0])


checkpoint bactmap_samplesheet:
    input:
        'results/.mapping_cache',
        sample_info='resources/samplesheet-with-ref.csv',
    output:
        samplesheet='resources/samplesheets/bactmap_{species}.csv',
    params:
        # Files
        db=data / 'mapping_cache.duckdb',
        input=input['sample_info'],
        output=output['samplesheet'],
        
        # Params
        reference_genome=lambda wildcards: wildcards.species,
    log:
        'logs/smk/bactmap_samplesheet_{species}.log'
    run:
        transform(
            models['nfcore_inputs']['bactmap'],
            params, db=params['db'], log=log[0], readonly=True
        )


def trimmed_reads(wildcards):
    return f'results/bactmap/{wildcards.species}/fastp/*.trim.fastq.gz'


checkpoint srst2_samplesheet:
    input:
        rules.all_bactmap.input,
        sample_info='resources/samplesheet-with-ref.csv',
    output:
        samplesheet='resources/samplesheets/srst2_{species}.csv',
    params:
        glob=lambda wildcards: trimmed_reads(wildcards),
        pat='fastp/(\\d+)_(1|2).trim.fastq.gz$',
        output=output['samplesheet']
    log:
        'logs/smk/srst2_samplesheet_{species}.log'
    run:
        transform(models['srst2']['samplesheet'], params, log=log[0])


use rule srst2 from widevariant with:
    input:
        multiext(
            'results/bactmap/{species}/fastp/{sample}',
            '_1.trim.fastq.gz', '_2.trim.fastq.gz',
        ),


checkpoint sarek_samplesheet:
    input:
        'results/.mapping_cache',
        sample_info='resources/samplesheet-with-ref.csv',
    output:
        samplesheet='resources/samplesheets/sarek_{species}.csv',
    params:
        # Files
        db=data / 'mapping_cache.duckdb',
        input=input['sample_info'],
        output=output['samplesheet'],
        
        # Params
        reference_genome=lambda wildcards: wildcards.species,
    log:
        'logs/smk/sarek_samplesheet_{species}.log'
    run:
        transform(
            models['nfcore_inputs']['sarek'],
            params, db=params['db'], log=log[0], readonly=True
        )


checkpoint snippy_samplesheet:
    input:
        sample_info='resources/samplesheet-with-ref.csv',
    output:
        samplesheet='resources/samplesheets/snippy.csv',
    log:
        'logs/smk/snippy_samplesheet.log'
    localrule: True
    shell:
        'cp {input} {output}'


def variants_db_path(wildcards):
    return (
        data / 'variants' / 
        f'species={wildcards.species}' / 
        f'family={wildcards.family}' / 
        f'tool={wildcards.mapping_tool}' / 
        'all.duckdb'
    )


rule create_variants_db:
    output:
        data / '.variants/species={species}/family={family}/tool={mapping_tool}/.done',
    params:
        db=lambda wildcards: variants_db_path(wildcards),
        contigs=seeds['types']['contigs'],
        species=lambda wildcards: wildcards.species,
    log:
        'logs/smk/create_variants_db_{species}_{family}_{mapping_tool}.log'
    run:
        transform(
            models['variants']['create_db'],
            params, db=params['db'], log=log[0]
        )
        shell('touch ' + output[0])


rule update_variants_db:
    input:
        rules.all_mapping.input,
        data / '.variants_{species}_{family}_{mapping_tool}_duckdb',
    output:
        'results/pseudogenomes/species={species}/family={family}/tool={mapping_tool}/snvs.parquet',
    params:
        # Files
        db=lambda wildcards: variants_db_path(wildcards),
        max_temp_directory_size='300GB',
        output=output[0],
        vcf_pq_glob=lambda wildcards: data / 'variants' / f'tool={wildcards.mapping_tool}' / f'species={wildcards.species}' / f'family={wildcards.family}' / '**/*.cleaned.parquet',

        # Params
        frac_cohort_core=config['params']['variant_filter']['frac_cohort_core'],
    log:
        'logs/smk/update_variants_db_{species}_{family}_{mapping_tool}.log'
    run:
        transform(
            models['variants']['insert_db'],
            params, db=params['db'], log=log[0]
        )