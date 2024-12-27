def ref_and_pe_fastqs(wildcards):
    import pandas as pd

    # TODO: May want to refactor dependency on samplesheet
    #       by adding fastq paths to reference_identification
    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[0]
    )

    samplesheet = samplesheet[samplesheet['sample'] == int(wildcards.sample)]

    species, r1, r2 = (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .merge(samplesheet, on='sample', how='inner')
        .filter(['reference_genome', 'fastq_1', 'fastq_2'])
        .values
        .flatten()
    )

    ref = config['public_data']['reference'][species]

    return {'ref': ref, 'R1': r1, 'R2': r2}


rule snippy:
    input:
        unpack(ref_and_pe_fastqs)
    output:
        multiext(
            (data_dir / 'results/snippy/{species}/variants/{sample}/{sample}.snps').as_posix(),
            '.aligned.fa', '.bed', '.csv', '.filt.vcf', '.gff', '.html', '.log', '.subs.vcf', '.tab', '.txt', '.vcf',
        )
    params:
        # Dirs
        outdir=lambda wildcards: data_dir / 'results/snippy' / wildcards.species / 'variants' / wildcards.sample,

        # Pipeline params
        mapqual=30,
        basequal=20,
        mincov=10,
        minfrac=".95",
        minqual=25,
        extra="--quiet --cleanup",
    log:
        log_dir / 'smk/mapping/snippy_{species}/{sample}.log'
    resources:
        cpus_per_task=4,
        mem_mb=2_000,
        njobs=1,
        runtime=15,
    envmodules:
        'snippy/4.6.0'
    shell:
        '''
        snippy \
          --force \
          --mapqual {params.mapqual} \
          --basequal {params.basequal} \
          --mincov {params.mincov} \
          --minfrac {params.minfrac} \
          --minqual {params.minqual} \
          --cpus {resources.cpus_per_task} \
          --ram "$(({resources.mem_mb} / 1200))" \
          --R1 {input.R1} \
          --R2 {input.R2} \
          --ref {input.ref} \
          --outdir {params.outdir} \
          --prefix "{wildcards.sample}.snps" \
          {params.extra}
        '''


def aggregate_snippy_vcfs(wildcards):
    import pandas as pd

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
        .apply(lambda df: (data_dir / 'results/snippy/{species}/variants/{sample}/{sample}.snps.vcf'.format(**df.to_dict())).as_posix())
        .values
        .flatten()
    )


# PHONY
rule all_snippy_vcfs:
    input:
        aggregate_snippy_vcfs
