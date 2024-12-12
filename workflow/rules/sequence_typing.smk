import json


with open(workflow.source_path('../../config/genome_size.json'), 'r') as f:
    genome_sizes = json.load(f)


checkpoint bactopia_samplesheet:
    input:
        'resources/samplesheets/main.csv',
        'resources/reference_genomes.csv',
    output:
        'resources/samplesheets/bactopia_{species}.txt'
    params:
        runtype='paired-end',
        extra=None
    log:
        'logs/smk/sequence_typing/bactopia_samplesheet_{species}.log'
    resources:
        njobs=50
    run:
        import pandas as pd

        pd.set_option('future.no_silent_downcasting', True)


        samplesheet = pd.read_csv(input[0])
        identification = pd.read_csv(input[1])
        
        species_mask = identification['species'] == wildcards.species

        (
            identification
            [species_mask]
            .merge(samplesheet.drop('species', axis=1), on='sample')
            .assign(
                runtype=params.runtype,
                extra=params.extra,
                genome_size=lambda df: df['species'].replace(genome_sizes)
            )
            .rename(columns={'fastq_1': 'r1', 'fastq_2': 'r2'})
            .filter(['sample', 'runtype', 'genome_size', 'species', 'r1', 'r2', 'extra'])
            .to_csv(output[0], sep='\t', index=False)
        )


rule bactopia:
    input:
        'resources/samplesheets/bactopia_{species}.txt'
    output:
        'results/bactopia/{species}/merged-results/mlst.tsv'
    params:
        cfg=config['sequence_typing']['bactopia']['nfconfig'],
        conda_pkgs=config['sequence_typing']['bactopia']['conda_pkgs_dirs'],
        covg=config['sequence_typing']['bactopia']['coverage'],
        data=config['sequence_typing']['bactopia']['datasets_cache'],
        outdir='results/bactopia/{species}',
        profile=config['sequence_typing']['bactopia']['profiles']
    log:
        'logs/smk/sequence_typing/bactopia_{species}.log'
    handover: True
    resources:
        njobs=200
    envmodules:
        'bactopia/3.1.0',
        'apptainer/1.3.2',
        'nextflow/24.04.4',
        'jdk/17.0.10',
        'anaconda/3/2023.03'
    container:
        'docker://bactopia/bactopia:latest'
    shell:
        '''
        export CONDA_PKGS_DIRS={params.conda_pkgs}
        mkdir -p $CONDA_PKGS_DIRS

        CONDADIR="/ptmp/$USER/.bactopia/conda"
        mkdir -p $CONDADIR

        NXF_WORKDIR="logs/nxf/bactopia/{wildcards.species}/work"
        mkdir -p $NXF_WORKDIR

        bactopia \
          -profile {params.profile} \
          -work-dir $NXF_WORKDIR \
          -qs 295 \
          --max_cpus 8 \
          --condadir $CONDADIR \
          --datasets_cache {params.data} \
          --samples {input} \
          --coverage {params.covg} \
          --outdir {params.outdir} \
          --nfconfig {params.cfg}

        # TODO: Must ensure is most recently created directory
        mv results/bactopia/{wildcards.species}/bactopia-runs/bactopia-*/merged-results results/bactopia/{wildcards.species}
        '''


def agg_bactopia(wildcards):
    import pandas as pd

    references = pd.read_csv(
        checkpoints.reference_identification
        .get(**wildcards)
        .output[0]
    )
    
    species = [
        spp for spp in references['species'].unique()
        if spp in config['wildcards']['species'].split('|')
    ]

    return expand(
        'results/bactopia/{species}/merged-results/mlst.tsv',
        species=species
    )


rule bactopia_mlst:
    input:
        agg_bactopia

