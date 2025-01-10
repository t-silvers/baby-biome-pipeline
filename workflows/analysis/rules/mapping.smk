import pathlib

include: './legacy-mapping.smk'


rule bactmap:
    input:
        'resources/samplesheets/bactmap_{species}.csv',
    output:
        'results/bactmap/{species}/pipeline_info/pipeline_report.html',
        'results/bactmap/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
        'results/bactmap/{species}/multiqc/multiqc_data/multiqc_samtools_stats_samtools.yaml',
        'results/bactmap/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml',
    params:
        pipeline=config['params']['bactmap']['pipeline'],

        # Dirs
        outdir=subpath(output[0], ancestor=2),
        workdir=lambda wildcards: f'logs/nxf/bactmap_{wildcards.species}_work',
        
        # Generic params
        config=config['params']['bactmap']['config'],
        profile=config['params']['bactmap']['profiles'],
        version=config['params']['bactmap']['version'], # NOTE: Unused when running downloaded pipeline

        # Pipeline params
        extra=config['params']['bactmap']['extra'],
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
    log:
        'logs/smk/bactmap_{species}.log'
    handover: True
    envmodules:
        'apptainer',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run {params.pipeline} \
          -config {params.config} \
          -profile {params.profile} \
          -resume \
          -work-dir {params.workdir} \
          --input {input} \
          --outdir {params.outdir} \
          --reference {params.reference} \
          {params.extra}
        '''


rule sarek:
    input:
        'resources/samplesheets/sarek_{species}.csv',
    output:
        'results/sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml',
    params:
        pipeline=config['params']['sarek']['pipeline'],

        # Dirs
        outdir=subpath(output[0], ancestor=2),
        workdir=lambda wildcards: f'logs/nxf/sarek_{wildcards.species}_work',

        # Generic params
        config=config['params']['sarek']['config'],
        profile=config['params']['sarek']['profiles'],
        version=config['params']['sarek']['version'],

        # Pipeline params
        extra=config['params']['sarek']['extra'],
        fasta=lambda wildcards: config['public_data']['reference'][wildcards.species],
        tools=config['params']['sarek']['tools'],
    log:
        'logs/smk/mapping/sarek_{species}.log'
    handover: True
    envmodules:
        'apptainer',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run {params.pipeline} \
          -config {params.config} \
          -profile {params.profile} \
          -resume \
          -work-dir {params.workdir} \
          --input {input} \
          --outdir {params.outdir} \
          --igenomes_ignore \
          --genome null \
          --fasta {params.fasta} \
          --fasta_fai {params.fasta}.fai \
          --known_indels false \
          --known_snps false \
          --tools {params.tools} \
          {params.extra}
        '''


# TODO: Use `snippy-multi`-generated submission script instead?
def ref_and_pe_fastqs(wildcards):
    import pandas as pd

    sample = int(wildcards.sample)
    
    species, r1, r2 = (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .query('sample == @sample')
        .filter(['reference_genome', 'fastq_1', 'fastq_2'])
        .values
        .flatten()
    )

    if species != wildcards.species:
        raise KeyError

    ref = config['public_data']['reference'][wildcards.species]

    return {'ref': ref, 'R1': r1, 'R2': r2}


rule snippy:
    """Call variants using Snippy.
    
    Args
      --mapqual     soft threshold for minimum mapping quality
      --basequal    soft threshold for minimum base quality
      --mincov      hard threshold for minimum coverage
      --minfrac     hard threshold for minimum fraction
      --minqual     soft threshold for minimum QUAL field
      --force       overwrite existing output; required since smk will create 
                    the `outdir` on init. Will also overwrite the "temp" 
                    directory if interrupted and resumed.
    """
    input:
        unpack(ref_and_pe_fastqs)
    output:
        multiext(
            'results/snippy/{species}/variants/{sample}.snps',
            '.aligned.fa', '.bed', '.csv', '.filt.vcf', '.gff', '.html', '.log', '.raw.vcf', '.subs.vcf', '.tab', '.txt', '.vcf',
        )
    params:
        # Dirs
        outdir=subpath(output[0], strip_suffix='.snps.aligned.fa'),

        # Pipeline params
        report='--report' in config['params']['snippy']['extra'],
        basequal=config['params']['snippy']['basequal'],
        extra=config['params']['snippy']['extra'],
        mapqual=config['params']['snippy']['mapqual'],
        mincov=config['params']['snippy']['mincov'],
        minfrac=config['params']['snippy']['minfrac'],
        minqual=config['params']['snippy']['minqual'],
    log:
        'logs/smk/snippy_{species}_{sample}.log'
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