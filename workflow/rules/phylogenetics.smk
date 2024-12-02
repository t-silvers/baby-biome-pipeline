rule:
    input:
        db='data/sample_info.duckdb',
        calls='results/{species}/pseudogenomes/annot_filtered_calls.csv',
    output:
        'results/{species}/pseudogenomes/msa.fas',
    resources:
        cpus_per_task=12,
        mem_mb=16_000,
        runtime=5,
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export  MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
                FILTERED_CALLS={input.calls}

        duckdb -readonly -init config/duckdbrc-slurm {input.db} \
          -c ".read workflow/scripts/calls_to_msa.sql" > {output}
        '''


rule filter_invariant_sites:
    input:
        'results/{species}/pseudogenomes/msa.fas',
    output:
        'results/{species}/pseudogenomes/msa_filtered.fas',
    localrule: True
    envmodules:
        'snp-sites/2.5.1'
    shell:
        '''
        snp-sites -o {output} {input}
        '''


rule raxml_ng:
    input:
        'results/{species}/pseudogenomes/msa_filtered.fas',
    params:
        extra='--all --model GTR+G --bs-trees 200',
        # outgroup=lambda wildcards: '--outgroup ' + config['outgroup'][wildcards.donor][wildcards.species]['ID'],
        outgroup='',
        prefix='results/{species}/raxml/msa_filtered',
    output:
        multiext(
            'results/{species}/raxml/msa_filtered.raxml',
            '.reduced.phy',
            '.rba',
            '.bestTreeCollapsed',
            '.bestTree',
            '.mlTrees',
            '.support',
            '.bestModel',
            '.bootstraps',
            '.log'
        )
    resources:
        cpus_per_task=48,
        mem_mb=64_000,
        runtime=120,
    envmodules:
        'raxml-ng/1.2.2_MPI'
    shell:
        '''
        export OMP_PLACES=threads

        raxml-ng \
          {params.extra} \
          {params.outgroup} \
          --msa {input} \
          --threads {resources.cpus_per_task} \
          --prefix {params.prefix} \
          --redo
        
        touch {output}
        '''
