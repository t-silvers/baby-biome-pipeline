def all_fastqs(wildcards):
    import pandas as pd

    return (
        pd.read_csv(
            checkpoints.samplesheet.get(**wildcards).output[1]
        )
        .filter(like='fastq', axis=1)
        .melt(value_name='path')
        .dropna()
        .drop_duplicates(subset=['path'])
        ['path']
        .values
        .flatten()
    )


rule sequencing_qc_falco:
    input:
        all_fastqs
    output:
        touch('results/qc/falco/sequencing_qc.done')
    params:
        outdir='results/qc/falco'
    resources:
        cpus_per_task=8,
        mem_mb=4_000,
        runtime=45
    envmodules:
        'falco/1.2.5'
    shell:
        'falco --outdir {params.outdir} --threads {resources.cpus_per_task} {input}'


def paired_fastqs(wildcards):
    import pandas as pd

    return (
        pd.read_csv(
            checkpoints.fastqs.get(**wildcards).output[0]
        )
        .dropna()
        .drop_duplicates(subset=['fastq_1', 'fastq_2'])
        .query(
            f"sample == {wildcards['sample']}"
        )
        .filter(like='fastq', axis=1)
        .values
        .flatten()
    )


rule read_adapter_trimming:
    input:
        paired_fastqs
    output:
        



        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            $out_fq1 \\
            $out_fq2 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """

        def fail_fastq = params.save_trimmed_fail ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.trim.fastq.gz \\
            --out2 ${prefix}_2.trim.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $options.args \\
            2> ${prefix}.fastp.log

        echo \$(fastp --version 2>&1) | sed -e "s/fastp //g" > ${software}.version.txt
        """



--cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50