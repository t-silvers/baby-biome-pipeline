# NOTE: `threads:` directive is parsed by several wrappers.

legacy_mapping_params = config['params']['legacy_mapping']


def fastqs(wildcards):
    import pandas as pd

    sample = int(wildcards.sample)

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .query('sample == @sample')
        .filter(['fastq_1', 'fastq_2'])
        .values
        .flatten()
    )


# TODO: Not included in dependency graph
# TODO: Move to qc.smk
rule fastqc:
    input:
        fastqs
    output:
        multiext(
            'results/legacy_mapping/{species}/fastqc/{sample}',
            '_R1_fastqc.html', '_R1_fastqc.zip', '_R2_fastqc.html', '_R2_fastqc.zip'
        ),
    params:
        outdir=subpath(output[0], parent=True),
        extra='--quiet',
    envmodules:
        'widevariant-legacy/1.0'
    shell:
        # TODO: A bit hacky...
        '''
        TMPDIR=$(mktemp -d -p {params.outdir} -t tmp.XXXXXX)

        fastqc {params.extra} --outdir "$TMPDIR" {input}

        for ext in html zip; do
          mv "$TMPDIR"/*_R1*_fastqc.$ext "{params.outdir}/{wildcards.sample}_R1_fastqc.$ext"
          mv "$TMPDIR"/*_R2*_fastqc.$ext "{params.outdir}/{wildcards.sample}_R2_fastqc.$ext"
        done
        rm -rf "$TMPDIR"
        '''


rule cutadapt_pe:
    input:
        fastqs
    output:
        fastq1='results/legacy_mapping/{species}/prepare/trimmed.{sample}_R1.fastq',
        fastq2='results/legacy_mapping/{species}/prepare/trimmed.{sample}_R2.fastq',
        qc='results/legacy_mapping/{species}/prepare/trimmed.{sample}.qc.txt',
    params:
        adapters=legacy_mapping_params['cutadapt_pe']['adapters'],
        extra=legacy_mapping_params['cutadapt_pe']['extra'],
    # NOTE: Wrapper looks for a `threads` attribute
    threads: 4
    envmodules:
        'cutadapt/4.9'
    wrapper:
        'v3.13.8/bio/cutadapt/pe'


rule sickle_pe:
    """ using sickle.

    Args:
        -g:
            If passed, output is gzipped.
        --quality,-q: (int)
            Threshold for trimming based on average quality in a window
        --readlen,-l: (int)
            Threshold to keep a read based on length after trimming
    """
    input:
        r1='results/legacy_mapping/{species}/prepare/trimmed.{sample}_R1.fastq',
        r2='results/legacy_mapping/{species}/prepare/trimmed.{sample}_R2.fastq',
    output:
        r1='results/legacy_mapping/{species}/prepare/processed.{sample}_R1.fastq.gz',
        r2='results/legacy_mapping/{species}/prepare/processed.{sample}_R2.fastq.gz',
        rs='results/legacy_mapping/{species}/prepare/processed.{sample}.single.fastq',
    params:
        qual_type=legacy_mapping_params['sickle_pe']['qual_type'],
        extra=legacy_mapping_params['sickle_pe']['extra'],
    envmodules:
        'sickle/1.33'
    wrapper:
        'v3.13.8/bio/sickle/pe'


def bowtie2_align_input(wildcards):
    return {
        'sample': [
            f'results/legacy_mapping/{wildcards.species}/prepare/processed.{wildcards.sample}_R1.fastq.gz',
            f'results/legacy_mapping/{wildcards.species}/prepare/processed.{wildcards.sample}_R2.fastq.gz',
        ],
        'idx': multiext(
            config['public_data']['reference-bowtie2'][wildcards.species],
            '.1.bt2',
            '.2.bt2',
            '.3.bt2',
            '.4.bt2',
            '.rev.1.bt2',
            '.rev.2.bt2',
        )
    }


# TODO: By default, bowtie2 prints a SAM header with @HD, @SQ and @PG lines. 
# When one or more --rg arguments are specified, bowtie2 will also print an 
# @RG line that includes all user-specified --rg tokens separated by tabs.
rule bowtie2_align:
    """
    Notes:
        -I/--minins <int> Default: 0 (essentially imposing no minimum)
            The minimum fragment length for valid paired-end alignments. E.g. 
            if -I 60 is specified and a paired-end alignment consists of two 
            20-bp alignments in the appropriate orientation with a 20-bp gap 
            between them, that alignment is considered valid (as long as -X 
            is also satisfied). A 19-bp gap would not be valid in that case. 
            If trimming options -3 or -5 are also used, the -I constraint is 
            applied with respect to the untrimmed mates.

            The larger the difference between -I and -X, the slower Bowtie 2 
            will run. This is because larger differences between -I and -X 
            require that Bowtie 2 scan a larger window to determine if a concordant 
            alignment exists. For typical fragment length ranges (200 to 400 
            nucleotides), Bowtie 2 is very efficient.
        
        -X/--maxins <int>
            The maximum fragment length for valid paired-end alignments. E.g. 
            if -X 100 is specified and a paired-end alignment consists of two 
            20-bp alignments in the proper orientation with a 60-bp gap between 
            them, that alignment is considered valid (as long as -I is also 
            satisfied). A 61-bp gap would not be valid in that case. If trimming 
            options -3 or -5 are also used, the -X constraint is applied with 
            respect to the untrimmed mates, not the trimmed mates.

            The larger the difference between -I and -X, the slower Bowtie 2 
            will run. This is because larger differences between -I and -X 
            require that Bowtie 2 scan a larger window to determine if a 
            concordant alignment exists. For typical fragment length ranges 
            (200 to 400 nucleotides), Bowtie 2 is very efficient.
    """
    input:
        unpack(bowtie2_align_input),
    output:
        'results/legacy_mapping/{species}/samtools/mapped.{sample}.bam',
        metrics='results/legacy_mapping/{species}/samtools/mapped.{sample}.metrics.txt',
    params:
        extra=legacy_mapping_params['bowtie2_align']['extra'],
    # NOTE: Wrapper looks for a `threads` attribute
    # NOTE: Use at least two threads
    threads: 4
    envmodules:
        'widevariant-legacy/1.0'
    # TODO: Correct parsing of `params.extra` in wrapper by func `snakemake_wrapper_utils.samtools.get_samtools_opts`?
    wrapper:
        'v3.13.8/bio/bowtie2/align'


rule samtools_sort:
    input:
        'results/legacy_mapping/{species}/samtools/mapped.{sample}.bam',
    output:
        'results/legacy_mapping/{species}/samtools/mapped.{sample}.sorted.bam',
    envmodules:
        'widevariant-legacy/1.0'
    shell:
        'samtools view -bSF4 {input} | samtools sort - -o {output}'


rule picard_add_read_groups:
    """Add read groups.

    Notes:
        - Read groups are missing from headers. picard>=3 (?) require @RG in headers.
        - The GATK CLI is transitioning
        - Wrapper uses deprecated opts that are breaking (namely, `get_java_opts`)

    Args:
        --INPUT,-I <String>
            Input file (BAM or SAM or a GA4GH url).  Required. 
        --OUTPUT,-O <File>
            Output file (SAM, BAM or CRAM).  Required. 
        --RGLB,-LB <String>
            Read-Group library  Required. 
        --RGPL,-PL <String>
            Read-Group platform (e.g. ILLUMINA, SOLID)  Required. 
        --RGPU,-PU <String>
            Read-Group platform unit (eg. run barcode)  Required. 
            Platform unit (e.g., flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.
        --RGSM,-SM <String>
            Read-Group sample name  Required.
            Use pool name where a pool is being sequenced. (see Refs)
    
    Reference:
        - https://samtools.github.io/hts-specs/SAMv1.pdf
        - https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    """
    input:
        'results/legacy_mapping/{species}/samtools/mapped.{sample}.sorted.bam',
    output:
        'results/legacy_mapping/{species}/samtools/fixed-rg.{sample}.sorted.bam',
    params:
        extra=lambda wildcards: '-RGLB lib1 -RGPL ILLUMINA -RGPU ' + wildcards.sample + ' -RGSM mpg_L31206',
    envmodules:
        'widevariant-legacy/1.0'
    # wrapper:
    #     'v3.13.8/bio/picard/addorreplacereadgroups'
    shell:
        'picard AddOrReplaceReadGroups -I {input} -O {output} {params.extra}'


# TODO: update snakemake version in container to use wrapper.
rule picard_markduplicates:
    input:
        bams='results/legacy_mapping/{species}/samtools/fixed-rg.{sample}.sorted.bam',
    output:
        bam='results/legacy_mapping/{species}/samtools/dedup.{sample}.sorted.bam',
        metrics='results/legacy_mapping/{species}/samtools/dedup.{sample}.metrics.txt',
    params:
        extra=legacy_mapping_params['picard_markduplicates']['extra'],
    envmodules:
        'widevariant-legacy/1.0'
    # wrapper:
    #     'v3.13.8/bio/picard/markduplicates'
    shell:
        'picard MarkDuplicates -I {input} -O {output.bam} {params.extra} --METRICS_FILE {output.metrics}'


rule samtools_index:
    input:
        'results/legacy_mapping/{species}/samtools/dedup.{sample}.sorted.bam',
    output:
        'results/legacy_mapping/{species}/samtools/dedup.{sample}.sorted.bam.bai',
    threads: 4  # This value - 1 passed to `-@` arg
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.13.8/bio/samtools/index'


def samtools_mpileup_input(wildcards):
    return {
        'bam': f'results/legacy_mapping/{wildcards.species}/samtools/dedup.{wildcards.sample}.sorted.bam',
        'reference_genome': config['public_data']['reference'][wildcards.species],
    }


rule samtools_mpileup:
    input:
        unpack(samtools_mpileup_input),
    output:
        'results/legacy_mapping/{species}/samtools/{sample}.mpileup.gz',
    params:
        extra=legacy_mapping_params['samtools_mpileup']['extra'],
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.13.8/bio/samtools/mpileup'


def bcftools_mpileup_input(wildcards):
    return {
        'alignments': f'results/legacy_mapping/{wildcards.species}/samtools/dedup.{wildcards.sample}.sorted.bam',
        'ref': config['public_data']['reference'][wildcards.species],
        'index': config['public_data']['reference'][wildcards.species] + '.fai',
    }


rule bcftools_mpileup:
    """
    Note that using "samtools mpileup" to generate BCF or VCF files has been
    removed.  To output these formats, please use "bcftools mpileup" instead.
    """
    input:
        unpack(bcftools_mpileup_input),
    output:
        pileup='results/legacy_mapping/{species}/variants/{sample}.mpileup.bcf',
    params:
        uncompressed_bcf=False,
        extra=legacy_mapping_params['bcftools_mpileup']['extra'],
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.14.0/bio/bcftools/mpileup'


rule bcftools_call:
    input:
        'results/legacy_mapping/{species}/variants/{sample}.mpileup.bcf',
    output:
        'results/legacy_mapping/{species}/variants/{sample}.calls.bcf',
    params:
        uncompressed_bcf=False,
        
        # NOTE: valid options include -c/--consensus-caller or -m/--multiallelic-caller
        caller=legacy_mapping_params['bcftools_call']['caller'],
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.13.8/bio/bcftools/call'


rule bcftools_query:
    input:
        'results/legacy_mapping/{species}/variants/{sample}.calls.bcf',
    output:
        'results/legacy_mapping/{species}/variants/{sample}.variant_quals.csv',
    envmodules:
        'widevariant-legacy/1.0'
    shell:
        'bcftools query -f "%POS,%REF,%ALT,%FQ\n" {input} > {output}'


rule bcftools_view:
    input:
        'results/legacy_mapping/{species}/variants/{sample}.calls.bcf',
    output:
        'results/legacy_mapping/{species}/variants/{sample}.calls.view.vcf.gz',
    params:
        extra=legacy_mapping_params['bcftools_view']['extra'],
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.13.8/bio/bcftools/view'


# TODO: Not included in dependency graph
rule tabix_index:
    input:
        'results/legacy_mapping/{species}/variants/{sample}.calls.view.vcf.gz',
    output:
        'results/legacy_mapping/{species}/variants/{sample}.calls.view.vcf.gz.tbi',
    params:
        legacy_mapping_params['tabix_index'],
    envmodules:
        'widevariant-legacy/1.0'
    wrapper:
        'v3.13.8/bio/tabix/index'