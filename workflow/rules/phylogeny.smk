# rule filter_variant_calls:
#     input:
#         'results/variants/{species}.duckdb',
#     output:
#         'results/pseudogenomes/species={species}/family={family}/filtered_calls.parquet',
#     params:
#         dp=config['pseudogenome']['total_read_depth'],
#         maf=config['pseudogenome']['maf'],
#         qual=config['pseudogenome']['quality'],
#         sdp=config['pseudogenome']['strand_read_depth'],
#     resources:
#         cpus_per_task=16,
#         mem_mb=48_000,
#         runtime=20,
#         njobs=1
#     envmodules:
#         'duckdb/nightly'
#     shell:
#         (
#             'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
#             'export FAMILY={wildcards.family};' +
#             'export DP={params.dp};' +
#             'export SNPGAP=3;' +
#             'export MAF={params.maf};' +
#             'export QUAL={params.qual};' +
#             'export STRAND_DP={params.sdp};' +
#             'duckdb -readonly -init ' + 
#             workflow.source_path('../../config/duckdbrc-slurm') +
#             ' {input} -c ".read ' + 
#             workflow.source_path('../scripts/filter_genotype_calls.sql') + 
#             '" > {output}'
#         )


# rule calls_to_msas:
#     input:
#         'results/pseudogenomes/species={species}/family={family}/filtered_calls.parquet',
#     output:
#         'results/pseudogenomes/species={species}/family={family}/annot_msa.duckdb',
#     resources:
#         cpus_per_task=12,
#         mem_mb=16_000,
#         runtime=5,
#         njobs=1
#     envmodules:
#         'duckdb/nightly'
#     shell:
#         (
#             'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
#             'export FILTERED_CALLS={input};' +
#             'duckdb -init ' + 
#             workflow.source_path('../../config/duckdbrc-slurm') +
#             ' {output} -c ".read ' + 
#             workflow.source_path('../scripts/models/msa.sql') + 
#             '"'
#         )


# rule calls_to_fasta:
#     input:
#         'results/pseudogenomes/species={species}/family={family}/annot_msa.duckdb',
#     output:
#         'results/pseudogenomes/species={species}_family={family}_msa.fas',
#     resources:
#         cpus_per_task=2,
#         mem_mb=4_000,
#         runtime=5,
#         njobs=1
#     envmodules:
#         'duckdb/nightly'
#     shell:
#         (
#             'export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB";' +
#             'duckdb -readonly -init ' + 
#             workflow.source_path('../../config/duckdbrc-slurm') +
#             ' {input} -c ".read ' + 
#             workflow.source_path('../scripts/calls_to_msa.sql') + 
#             '" > {output}'
#         )


# def collect_msas(wildcards, ):
#     import pandas as pd

#     samplesheet = pd.read_csv(
#         checkpoints.samplesheet
#         .get(**wildcards)
#         .output[1]
#     )

#     identification = pd.read_csv(
#         checkpoints.reference_identification
#         .get(**wildcards)
#         .output[0]
#     )

#     identification = identification[
#         identification['species'].isin(
#             config['wildcards']['species'].split('|')
#         )
#     ]

#     return list(
#         samplesheet
#         .filter(['family', 'sample'])
#         .merge(identification, on='sample')
#         .filter(['species', 'family'])
#         .dropna()
#         .drop_duplicates()
#         .transpose()
#         .apply(lambda df: 'results/pseudogenomes/species={species}_family={family}_msa.fas'.format(**df.to_dict()))
#         .values
#         .flatten()
#     )