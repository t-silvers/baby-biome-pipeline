include: '../phylogeny-modules/raxml_ng.smk'
include: '../phylogeny-modules/report.smk'


cohort_to_glob = {
    'B001-Bifidobacterium_adolescentis': "'../../../data/variants/species=Bifidobacterium_adolescentis/family=B001/id=*/library=*/*.cleaned.vcf.parquet'",
    'B001-Bifidobacterium_bifidum': "'../../../data/variants/species=Bifidobacterium_bifidum/family=B001/id=*/library=*/*.cleaned.vcf.parquet'",
    'B002-Bifidobacterium_bifidum': "'../../../data/variants/species=Bifidobacterium_bifidum/family=B002/id=*/library=240704_B002_B001_Lib_AVITI_reseq/*.cleaned.vcf.parquet'",
    'B001-Bifidobacterium_longum': "'../../../data/variants/species=Bifidobacterium_longum/family=B001/id=*/library=*/*.cleaned.vcf.parquet'",
    'B002-Bifidobacterium_longum': "'../../../data/variants/species=Bifidobacterium_longum/family=B002/id=*/library=240704_B002_B001_Lib_AVITI_reseq/*.cleaned.vcf.parquet'",
    # 'B001-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B001/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B002-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B002/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B006-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B006/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B009-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B009/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B012-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B012/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B013-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B013/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B016-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B016/id=*/library=*/*.cleaned.vcf.parquet'",
    # 'B017-Escherichia_coli': "'../../../data/variants/species=Escherichia_coli/family=B017/id=*/library=*/*.cleaned.vcf.parquet'",
}



rule msa:
    input:
        (data_dir / 'data/samplesheet.duckdb').as_posix(),
        aggregate_snvs_db
    output:
        multiext(
            workflow.source_path('results/trees/{cohort}/msa'),
            '-positions.parquet',
            '-metadata.txt',
            '-positions-varying.fas',
            '-varying.fas',
            '-positions.fas',
            '-len.fas',
            '.fas',
        )
    params:
        # TODO: Make filtering params here explicit
        duckdbrc=duckdbrc_slurm,
        model=config['models']['msa'],
        glob=lambda wc: cohort_to_glob[wc.cohort],
        trees_dir=workflow.source_path('results/trees/{cohort}'),
    resources:
        cpus_per_task=4,
        mem_mb=31_000,
        runtime=119,
        njobs=1
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
               VCF_PQS={params.glob};

        mkdir -p {params.trees_dir} && \
        cd {params.trees_dir} && \
        duckdb -readonly -init {params.duckdbrc} ../../../{input[0]} -c ".read ../../../{params.model}"
        '''


rule snp_sites:
    input:
        workflow.source_path('results/trees/{cohort}/msa.fas'),
    output:
        workflow.source_path('results/trees/{cohort}/msa-filtered.fas'),
    resources:
        njobs=1
    envmodules:
        'snp-sites/2.5.1'
    shell:
        'snp-sites -o {output} {input}'


def aggregate_trees(wildcards):
    return expand(
        workflow.source_path('results/trees/{cohort}/msa-filtered.raxml.bestTree'),
        cohort=list(cohort_to_glob.keys())
    )


rule all_trees:
    input:
        aggregate_trees