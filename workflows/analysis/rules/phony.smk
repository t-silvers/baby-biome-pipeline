import pathlib
from abc import ABC, abstractmethod
from typing import Union

import pandas as pd
from snakemake.checkpoints import Checkpoint, CheckpointJob


BRACKEN_TEMPLATE = 'identification/tool=taxprofiler/db={db_name}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.bracken.parquet'

SNIPPY_TEMPLATE = 'results/snippy/{species}/variants/{sample}.snps.aligned.fa'

SRST2_TEMPLATE = 'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}.txt'

BACTMAP_TEMPLATE = 'results/bactmap/{species}/pipeline_info/pipeline_report.html'

SAREK_TEMPLATE = 'results/sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml'

CLEANED_VCF_TEMPLATE = 'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.cleaned.parquet'


def aggregate_bactmap(wildcards) -> list[str]:
    return _aggregate_nfcore_mapping(BACTMAP_TEMPLATE, checkpoints.reference_identification, **wildcards)


def aggregate_bracken(wildcards) -> list[str]:
    return BrackenAggregator.from_data(read_samplesheet(checkpoints.samplesheet, **wildcards))


def aggregate_sarek(wildcards) -> list[str]:
    return _aggregate_nfcore_mapping(SAREK_TEMPLATE, checkpoints.reference_identification, **wildcards)


def aggregate_snippy(wildcards) -> list[str]:
    return SnippyAggregator.from_data(read_reference(checkpoints.reference_identification, **wildcards))


def aggregate_srst2(wildcards) -> list[str]:
    return SRST2Aggregator.from_data(read_reference(checkpoints.reference_identification, **wildcards))


def aggregate_vcfs(wildcards) -> list[str]:
    return VCFAggregator.from_data(read_reference(checkpoints.reference_identification, **wildcards))


rule all_bactmap:
    input:
        aggregate_bactmap


rule all_taxprofiler:
    input:
        aggregate_bracken


rule all_sarek:
    input:
        aggregate_sarek


rule all_snippy:
    input:
        aggregate_snippy


rule all_srst2:
    input:
        aggregate_srst2


rule all_mapping:
    input:
        aggregate_vcfs


def read_reference(rule: Checkpoint, **wildcards) -> pd.DataFrame:
    return read_checkpoint('samplesheet_with_reference', rule, **wildcards)


def read_samplesheet(rule: Checkpoint, **wildcards) -> pd.DataFrame:
    return read_checkpoint('samplesheet', rule, **wildcards)


def read_checkpoint(key: Union[str, int], rule: Checkpoint, **wildcards) -> pd.DataFrame:
    """Load output data file from checkpoint as dataframe."""
    job: CheckpointJob = rule.get(**wildcards)
    return pd.read_csv(job.output[key])


def _aggregate_nfcore_mapping(path_template: str, rule: Checkpoint, **wildcards):
    genomes = wc_params['species_ref']
    species = (
        read_reference(rule, **wildcards)
        .query('reference_genome in @genomes')
        ['reference_genome']
        .unique()
    )
    return expand(path_template, species=species)


class PathAggregator(ABC):
    """Base class for output path aggregators."""

    def __init__(self, directory: Union[str, pathlib.Path], sep: str = '|'):
        self.directory = pathlib.Path(directory)
        self.sep = sep
        self.path_templates = self.prepare_templates()

    @classmethod
    def from_data(cls, data: pd.DataFrame) -> list[str]:
        aggregator = cls()
        return aggregator(data)

    def __call__(self, data: pd.DataFrame) -> list[str]:
        """Format string template with values from a pandas DataFrame"""
        return (
            self.prepare_data(data)
            .drop_duplicates()
            .transpose()
            .apply(lambda df: self._fmt_paths(df))
            .str.split(self.sep)
            .explode()
            .values
            .flatten()
        )

    @abstractmethod
    def prepare_templates(self) -> list[str]:
        """Prepare path templates."""

    @abstractmethod
    def prepare_data(self, data) -> pd.DataFrame:
        """Prepare input data."""

    def _fmt_paths(self, df) -> str:
        return self.sep.join([(self.directory / t.format(**df.to_dict())).as_posix() for t in self.path_templates])


class BrackenAggregator(PathAggregator):
    """Collect bracken results from taxprofiler output."""

    _bracken_template = BRACKEN_TEMPLATE

    def __init__(self):
        super().__init__(data)
    
    @property
    def bracken_dbs(self) -> list[str]:
        """Get db names used for bracken."""
        return (
            pd.read_csv(
                config['params']['taxprofiler']['databases']
            )
            .query('tool == "bracken"')
            ['db_name']
            .unique()
            .tolist()
        )

    def prepare_templates(self) -> list[str]:
        return [self._bracken_template.format(db_name=x) for x in self.bracken_dbs]

    def prepare_data(self, data) -> pd.DataFrame:
        return data


class SRST2Aggregator(PathAggregator):
    """Collect srst2 results."""

    _srst2_template = SRST2_TEMPLATE

    def __init__(self):
        super().__init__(data)

    @property
    def sequence_type_schemas(self) -> list[str]:
        # TODO: Could replace with `mlst --list` check
        return list(config['public_data']['mlst'].keys())

    def prepare_templates(self) -> list[str]:
        return [self._srst2_template]

    def prepare_data(self, data) -> pd.DataFrame:
        return (
            data
            .filter(['reference_genome', 'family', 'id', 'library', 'sample'])
            .rename(columns={'reference_genome': 'species'})
            .query('species in @self.sequence_type_schemas')
            .dropna()
        )


class SnippyAggregator(PathAggregator):
    """Collect snippy results."""

    _snippy_template = SNIPPY_TEMPLATE

    def __init__(self):
        super().__init__('')

    def prepare_templates(self) -> list[str]:
        return [self._snippy_template]

    def prepare_data(self, data) -> pd.DataFrame:
        return (
            data
            .filter(['reference_genome', 'sample'])
            .rename(columns={'reference_genome': 'species'})
            .dropna()
        )


class VCFAggregator(PathAggregator):
    """Collect cleaned VCF parquet files from all mapping tools output."""

    _vcfs_template = CLEANED_VCF_TEMPLATE

    def __init__(self):
        super().__init__(data)
    
    @property
    def mapping_tools(self):
        return config['tools']['mapping'].split('|')
    
    def prepare_templates(self) -> list[str]:
        return [self._vcfs_template.format(mapping_tool=x) for x in self.mapping_tools]

    def prepare_data(self, data) -> pd.DataFrame:
        return data.rename(columns={'reference_genome': 'species'}).dropna()