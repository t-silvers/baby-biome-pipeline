import pathlib

snippy_report = '--report' in config['mapping']['snippy']['extra']


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

    # TODO: Remove samples that have already been analyzed; otherwise,
    #       relies on nxf dependency tracking. Use run db (sample,tool,...).

    return {'ref': ref, 'R1': r1, 'R2': r2}


def get_snippy_cpus(wildcards, attempt):
    if snippy_report:
        return 32
    else:
        return attempt * 4


def get_snippy_mem(wildcards, attempt):
    if snippy_report:
        return 16_000
    else:
        return attempt * 1_000


def get_snippy_time(wildcards, attempt):
    if snippy_report:
        return attempt * 30
    else:
        return 5


# TODO: Use `snippy-multi`-generated submission script?
rule snippy:
    """Call variants using Snippy.
    
    Args
      --mincov
        hard threshold for minimum coverage
      --minfrac
        hard threshold for minimum fraction
      --force
        overwrite existing output; required since smk will create the `outdir` on init.
        Will also overwrite the "temp" directory if interrupted and resumed.
    """
    input:
        unpack(ref_and_pe_fastqs)
    output:
        multiext(
            (results / 'snippy/{species}/variants/{sample}.snps').as_posix(),
            '.aligned.fa', '.bed', '.csv', '.filt.vcf', '.gff', '.html', '.log', '.subs.vcf', '.tab', '.txt', '.vcf',
        )
    params:
        # Dirs
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent / pathlib.Path(output[0]).stem.split('.')[0],

        # Pipeline params
        basequal=config['mapping']['snippy']['basequal'],
        extra=config['mapping']['snippy']['extra'],
        mapqual=config['mapping']['snippy']['mapqual'],
        mincov=config['mapping']['snippy']['mincov'], # Hard threshold
        minfrac=config['mapping']['snippy']['minfrac'], # Hard threshold
        minqual=config['mapping']['snippy']['minqual'],
    log:
        logdir / 'smk/mapping/snippy_{species}/{sample}.log'
    resources:
        cpus_per_task=get_snippy_cpus,
        mem_mb=16_000,
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
          {params.extra} && \

        cd {params.outdir} && \
        mv {wildcards.sample}.snps* .. && \
        cd .. && \
        rm -rf {params.outdir}
        '''