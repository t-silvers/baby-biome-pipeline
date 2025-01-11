module widevariant:
    snakefile: '../../analysis/Snakefile'
    config: config


use rule * from widevariant


localrules:
    bactmap,
    bactmap_samplesheet,
    reference_identification,
    sarek,
    sarek_samplesheet,
    taxprofiler,
    taxprofiler_samplesheet,


rule taxprofiler_samplesheet:
    input:
        'resources/samplesheet.csv',
        data / '.identification_duckdb',
    output:
        'resources/samplesheets/taxprofiler.csv',
    params:
        # Files
        db=data / 'identification.duckdb',
        input=input[0],
        output=output[0],
        
        # Params
        # NOTE: Must conform to EBI ENA controlled vocabulary
        instrument_platform='ILLUMINA',
    log:
        'logs/smk/taxprofiler_samplesheet.log'
    run:
        transform(models['taxprofiler_samplesheet'], params, db=params['db'], log=log[0], readonly=True)


# TODO: Specify db name
rule prepare_identification:
    input:
        rules.all_taxprofiler.input
    output:
        'resources/identification-raw.csv',
    params:
        bracken_glob=data / 'identification/tool=taxprofiler/db=*/family=*/id=*/library=*/*.bracken.parquet',
        output=output[0],
    log:
        'logs/smk/prepare_identification.log'
    run:
        transform(models['prepare_identification'], params, log=log[0])


rule update_identification_db:
    input:
        'resources/identification-raw.csv',
        data / '.identification_duckdb',
    output:
        'resources/reference-pp.csv',
    params:
        # Files
        db=data / 'identification.duckdb',
        input=input[0],
        output=output[0],
        
        # Params
        read_frac=config['params']['reference_genome']['minimum_fraction_reference'],
        read_pow=config['params']['reference_genome']['minimum_readspow_reference'],
    log:
        'logs/smk/update_identification_db.log'
    run:
        transform(models['update_identification_db'], params, db=params['db'], log=log[0])


checkpoint reference_identification:
    input:
        'resources/samplesheet.csv',
        'resources/reference-pp.csv',
    output:
        samplesheet_with_reference='resources/samplesheet-with-ref.csv',
    log:
        'logs/smk/reference_identification.log'
    run:
        import pandas as pd

        samplesheet = pd.read_csv(input[0])
        ref_genomes = pd.read_csv(input[1])
        
        merged = samplesheet.merge(ref_genomes, on='sample')

        def plate_and_seq_id_agree():
            ids = {
                spp: merged[spp].str.split('_', expand=True)[0]
                for spp in ['species_stdized', 'reference_genome']
            }
            return ids['species_stdized'] == ids['reference_genome']

        def reference_is_available():
            return merged['reference_genome'].isin(wc_params['species_ref'])

        mask = plate_and_seq_id_agree() & reference_is_available()
        merged[mask].to_csv(output[0], index=False)


rule bactmap_samplesheet:
    input:
        'resources/samplesheet-with-ref.csv',
    output:
        'resources/samplesheets/bactmap_{species}.csv',
    log:
        'logs/smk/bactmap_samplesheet_{species}.log'
    run:
        import pandas as pd

        BACTMAP_COLS = ['sample', 'fastq_1', 'fastq_2']

        (
            pd.read_csv(input[0])
            .query('reference_genome == @wildcards.species')
            .filter(BACTMAP_COLS)
            .to_csv(output[0], index=False)
        )


use rule srst2 from widevariant with:
    input:
        multiext(
            'results/bactmap/{species}/fastp/{sample}',
            '_1.trim.fastq.gz', '_2.trim.fastq.gz',
        ),


rule sarek_samplesheet:
    input:
        'resources/samplesheet-with-ref.csv',
    output:
        'resources/samplesheets/sarek_{species}.csv',
    log:
        'logs/smk/sarek_samplesheet_{species}.log'
    run:
        import pandas as pd

        SAREK_COLS = ['patient', 'sample', 'lane', 'fastq_1', 'fastq_2']

        (
            pd.read_csv(input[0])
            .query('reference_genome == @wildcards.species')

            # TODO: Consider other "patient" groupings; however, using
            #       e.g. `.rename(columns={'family': 'patient'})` will not
            #       work here (must instead be unique).
            .assign(patient=lambda df: df['sample'])

            # TODO: Check availability of lane information
            .assign(lane='lane1')

            .filter(SAREK_COLS)
            .to_csv(output[0], index=False)
        )


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
        frac_cohort_share_pos=1,
    log:
        'logs/smk/update_variants_db_{species}_{family}_{mapping_tool}.log'
    run:
        transform(models['update_variants_db'], params, db=params['db'], log=log[0])