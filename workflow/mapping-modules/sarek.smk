import pathlib

SAREK_SENTINEL = 'sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml'

SAREK_VARIANT_CALLING_TOOLS = {
    'mpileup': 'bcftools',
    'deepvariant': 'deepvariant',
    'freebayes': 'freebayes',
    'haplotypecaller': 'haplotypecaller',
}


rule sarek_samplesheet:
    input:
        results / 'samplesheets/samplesheet.csv',
        results / 'samplesheets/reference_genomes.csv',
        results / 'samplesheets/mapping_progress.csv',
    output:
        results / 'samplesheets/sarek_{species}.csv',
    log:
        logdir / 'smk/mapping/sarek_samplesheet_{species}.log'
    resources:
        njobs=1,
    run:
        import pandas as pd

        SAREK_COLS = ['patient', 'sample', 'lane', 'fastq_1', 'fastq_2']

        samplesheet = pd.read_csv(input[0])
        identification = pd.read_csv(input[1])

        species_mask = identification['reference_genome'] == wildcards.species

        sarek_input = (
            identification
            [species_mask]
            .filter(['sample'])
            .merge(samplesheet, on='sample')

            # TODO: Consider other "patient" groupings; however, using
            #       e.g. `.rename(columns={'family': 'patient'})` will not
            #       work here (must instead be unique).
            .assign(patient=lambda df: df['sample'])

            # TODO: Check availability of lane information
            .assign(lane='lane1')

            .filter(SAREK_COLS)
        )

        # NOTE: Remove samples that have already been analyzed; otherwise,
        #       relies on nxf dependency tracking, which will fail here.

        progress = pd.read_csv(input[2])

        # TODO: Check individual tools
        excluded = progress[progress['tool'].str.startswith('sarek_')]['sample'].to_list()

        sarek_input = sarek_input[
            ~sarek_input['sample'].isin(excluded)
        ]

        sarek_input.to_csv(output[0], index=False)


rule sarek:
    input:
        results / 'samplesheets/sarek_{species}.csv',
    output:
        results / SAREK_SENTINEL,
    params:
        pipeline=config['mapping']['sarek']['pipeline'],

        # Dirs
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent.parent,
        workdir=lambda wildcards: logdir / f'nxf/sarek_{wildcards.species}_work',

        # Generic params
        config=config['mapping']['sarek']['config'],
        profile=config['mapping']['sarek']['profiles'],
        version=config['mapping']['sarek']['version'],

        # Pipeline params
        extra=config['mapping']['sarek']['extra'],
        fasta=lambda wildcards: config['public_data']['reference'][wildcards.species],
        tools=config['mapping']['sarek']['tools'],
    log:
        logdir / 'smk/mapping/sarek_{species}.log'
    handover: True
    resources:
        njobs=max_submit,
    envmodules:
        # TODO: Dynamically resolve envmodules
        # 'apptainer/1.3.2',
        # 'apptainer/1.1.7',
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
          {params.extra} && \

        cd {params.outdir} && \
        mv variant_calling variants && \
        mv variants/*/*/*.vcf.gz variants
        '''


def aggregate_sarek(wildcards):
    import pandas as pd

    species = (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .query('reference_genome in @wc_params["species"]')
        ['reference_genome']
        .unique()
    )

    return expand(results / SAREK_SENTINEL, species=species)


# PHONY
rule all_sarek:
    input:
        aggregate_sarek


rule sarek_vcf:
    """Collect vcf output from sarek."""
    input:
        results / SAREK_SENTINEL,
    output:
        expand(
            results / 'sarek/{{species}}/variants/{{sample}}.{vc_tool}.vcf.gz',
            vc_tool=[
                SAREK_VARIANT_CALLING_TOOLS[tool]
                for tool in config['mapping']['sarek']['tools'].split(',')
                if tool in list(SAREK_VARIANT_CALLING_TOOLS.keys())
            ]
        )
    resources:
        njobs=1
    localrule: True