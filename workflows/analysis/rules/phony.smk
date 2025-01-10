def aggregate_bracken(wildcards) -> list[str]:
    import pandas as pd

    def bracken_pq(df) -> str:
        return data_path_from_template(
            'identification/tool=taxprofiler/db={{db_name}}/family={family}/id={id}/library={library}/{sample}.bracken.parquet',
            df.to_dict()
        )

    def bracken_db() -> dict:
        bracken_dbs = (
            pd.read_csv(config['params']['taxprofiler']['databases'])
            .query('tool == "bracken"')
        )

        if bracken_dbs.shape[0] != 1:
            raise ValueError('Expected exactly one bracken database.')

        return bracken_dbs.to_dict('records')[0]

    samplesheet = pd.read_csv(
        checkpoints.samplesheet
        .get(**wildcards)
        .output[0]
    )

    return (
        samplesheet
        .transpose()
        .apply(lambda df: bracken_pq(df).format(**bracken_db()))
        .values
        .flatten()
    )


def aggregate_srst2(wildcards) -> list[str]:
    import pandas as pd

    # TODO: Could replace with `mlst --list` check
    available_schemas = list(config['public_data']['mlst'].keys())

    def srst2_results(df) -> str:
        return data_path_from_template(
            'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}.txt', 
            df.to_dict()
        )

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .filter(['reference_genome', 'family', 'id', 'library', 'sample'])
        .rename(columns={'reference_genome': 'species'})
        .query('species in @available_schemas')
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: srst2_results(df))
        .values
        .flatten()
    )


def aggregate_bactmap(wildcards) -> list[str]:
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
        'results/bactmap/{species}/pipeline_info/pipeline_report.html',
        species=species
    )


def aggregate_sarek(wildcards) -> list[str]:
    import pandas as pd

    species = (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .query('reference_genome in @wc_params["species_ref"]')
        ['reference_genome']
        .unique()
    )

    return expand(
        'results/sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml',
        species=species
    )


def aggregate_snippy(wildcards) -> list[str]:
    import pandas as pd

    return (
        pd.read_csv(
            checkpoints.reference_identification
            .get(**wildcards)
            .output[0]
        )
        .filter(['reference_genome', 'sample'])
        .rename(columns={'reference_genome': 'species'})
        .dropna()
        .drop_duplicates()
        .transpose()
        .apply(lambda df: 'results/snippy/{species}/variants/{sample}.snps.aligned.fa'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule all_taxprofiler:
    input:
        aggregate_bracken


rule all_srst2:
    input:
        aggregate_srst2


rule all_bactmap:
    input:
        aggregate_bactmap


rule all_sarek:
    input:
        aggregate_sarek


rule all_snippy:
    input:
        aggregate_snippy