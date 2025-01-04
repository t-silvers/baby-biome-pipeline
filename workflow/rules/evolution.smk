# include: '../evolution-modules/raxml_ng.smk'
# include: '../evolution-modules/report.smk'


# TODO: Replace with checkpoint and/or parsing of provided sample sheet
COHORT_KEY = {
    'ecoli_bactmap_b012': {
        'glob': 'tool=bactmap/species=Escherichia_coli/family=B012',
        'species': 'Escherichia_coli'
    }
}


# checkpoint define_cohorts:
#     input:
#         results / 'samplesheets/samplesheet.csv',
#         results / 'samplesheets/reference_genomes.csv',
#         aggregate_vcfs,
#     output:
#         results / 'samplesheets/cohorts.csv',
#     params:
#         vars=['species', 'family']
#     log:
#         logdir / 'smk/evolution/define_cohorts.log'
#     resources:
#         njobs=1,
#     run:
#         import pandas as pd

#         # TODO: Replace columns with hive-partitioning-formatted string,
#         #       e.g., "species={species}/family={family}". Implement
#         #       sorting/ordering.

#         (
#             pd.read_csv(input[1])
#             .filter(['sample', 'reference_genome'])
#             .rename(columns={'reference_genome': 'species'})
#             .merge(
#                 pd.read_csv(input[0]),
#                 on='sample'
#             )
#             .filter(params['vars'])
#             .drop_duplicates()
#             .to_csv(output[0], index=False)
#         )


rule variants_db:
    output:
        data / 'pseudogenomes/species={species}/cohort={cohort}/variants.duckdb',
    params:
        model=workflow.source_path(models['vcfs']['create']),
        contigs_seed=workflow.source_path(seeds['types']['contigs']),
    log:
        logdir / 'smk/evolution/variants_db_tool_{species}_{cohort}.log'
    resources:
        njobs=1,
    run:
        params.update({'species': wildcards.species})
        transform(params['model'], params, db=output[0], log=log[0])


rule update_variants_db:
    input:
        data / 'pseudogenomes/species={species}/cohort={cohort}/variants.duckdb',
        aggregate_vcfs,
    output:
        results / 'variants/species={species}/cohort={cohort}/shared_snvs.parquet',
    params:
        model=workflow.source_path(models['vcfs']['insert']),
        max_temp_directory_size='300GB',
        
        # Params
        vcf_pq_glob=lambda wildcards: data / 'variants' / COHORT_KEY[wildcards.cohort]['glob'] / '**/*.cleaned.parquet',
        frac_cohort_share_pos=1,
    log:
        logdir / 'smk/evolution/update_variants_db_{species}_{cohort}.log'
    resources:
        cpus_per_task=lambda wildcards, attempt: attempt * 4,
        mem_mb=lambda wildcards, attempt: attempt * 4_000,
        runtime=lambda wildcards, attempt: attempt * 10,
        njobs=1,
    run:
        params.update({'output': output[0]})
        transform(params['model'], params, db=input[0], log=log[0])


def aggregate_snvs(wildcards):
    return [
        results / 'variants/species={species}/cohort={cohort}/shared_snvs.parquet'.format(
            species=COHORT_KEY[cohort]['species'], cohort=cohort
        ) for cohort in list(COHORT_KEY.keys())
    ]


# PHONY
rule all_snvs:
    input:
        aggregate_snvs
