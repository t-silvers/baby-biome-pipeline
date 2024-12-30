import pathlib

SAREK_VARIANT_CALLING_TOOLS = ['bcftools', 'deepvariant', 'freebayes', 'haplotypecaller']


rule sarek_samplesheet:
    input:
        results / 'samplesheets/samplesheet.csv',
        results / 'samplesheets/reference_genomes.csv',
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

        (
            identification
            [species_mask]
            .merge(samplesheet, on='sample')

            # TODO: Consider other "patient" groupings
            .rename(columns={'family': 'patient'})
            
            # TODO: Check availability of lane information
            .assign(lane='lane1')

            .filter(SAREK_COLS)
            .to_csv(output[0], index=False)
        )

        # TODO: Remove samples that have already been analyzed; otherwise,
        #       relies on nxf dependency tracking. Use run db (sample,tool,...).


rule sarek:
    input:
        results / 'samplesheets/sarek_{species}.csv',
    output:
        results / 'sarek/{species}/pipeline_info/pipeline_report.html',
    params:
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
        'apptainer/1.3.2',
        'nextflow/24.10',
        'jdk/17.0.10'
    shell:
        '''
        nextflow run nf-core/sarek \
          -config {params.config} \
          -profile {params.profile} \
          -r {params.version} \
          -resume \
          -work-dir {params.workdir} \
          --input {input} \
          --outdir {params.outdir} \
          --igenomes_ignore \
          --genome null \
          --fasta {params.fasta} \
          --tools {params.tools} \
          {params.extra}
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

    return expand(
        results / 'sarek/{species}/pipeline_info/pipeline_report.html',
        species=species
    )


# PHONY
rule all_sarek:
    input:
        aggregate_sarek


rule sarek_vcf:
    """Collect vcf output from sarek."""
    input:
        results / 'sarek/{species}/pipeline_info/pipeline_report.html',
    output:
        touch(expand(
            results / 'sarek/{{species}}/variantcalling/{vc_tool}/{{sample}}/{{sample}}.{vc_tool}.vcf.gz',
            vc_tool=[
                tool for tool in config['mapping']['sarek']['tools'].split(',') 
                if tool in SAREK_VARIANT_CALLING_TOOLS
            ]
        ))
    resources:
        njobs=1
    localrule: True