from abc import ABC, abstractmethod
import pathlib
from typing import Union

import pandas as pd
from snakemake.checkpoints import Checkpoint, CheckpointJob


BRACKEN_FIELDS = ['db_name', 'family', 'id', 'library', 'sample']
BRACKEN_TEMPLATE = 'identification/tool=taxprofiler/db={db_name}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.bracken.parquet'

SNIPPY_FIELDS = ['species', 'sample']
SNIPPY_TEMPLATE = 'results/snippy/{species}/variants/{sample}.snps.aligned.fa'

SRST2_FIELDS = ['species', 'family', 'id', 'library', 'sample']
SRST2_TEMPLATE = 'identification/tool=srst2/species={species}/family={family}/id={id}/library={library}/{sample}.txt'

# CLEANED_VCF_FIELDS = ['mapping_tool', 'species', 'family', 'id', 'library', 'sample']
# CLEANED_VCF_TEMPLATE = 'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.cleaned.parquet'

BACTMAP_FIELDS = ['species']
BACTMAP_TEMPLATE = 'results/bactmap/{species}/pipeline_info/pipeline_report.html'

SAREK_FIELDS = ['species']
SAREK_TEMPLATE = 'results/sarek/{species}/pipeline_info/nf_core_sarek_software_mqc_versions.yml'

VCF_FIELDS = ['mapping_tool', 'species', 'family', 'id', 'library', 'sample']
VCF_TEMPLATE = 'variants/tool={mapping_tool}/species={{species}}/family={{family}}/id={{id}}/library={{library}}/{{sample}}.parquet'


rule dummy_bracken:
    input:
        'resources/samplesheets/taxprofiler.csv',
    output:
        'results/.bracken.done',
    localrule: True
    shell:
        'touch {output}'


def aggregate_bracken(wildcards) -> list[str]:
    sample_info = read_sample_info(checkpoints.taxprofiler_samplesheet, **wildcards)
    samplesheet = read_samplesheet(checkpoints.taxprofiler_samplesheet, **wildcards)

    if not samplesheet.empty():
        data = samplesheet.filter(['sample']).merge(sample_info)
        return BrackenAggregator.from_data(data)
    
    else:
        return 'results/.bracken.done'


def aggregate_snippy(wildcards) -> list[str]:
    samplesheet = read_samplesheet(checkpoints.snippy_samplesheet, **wildcards)
    return SnippyAggregator.from_data(samplesheet)


def aggregate_srst2(wildcards) -> list[str]:
    srst2_schema = list(config['public_data']['mlst'].keys())
    data = pd.concat([
        _aggregate_over_wcs(checkpoints.srst2_samplesheet, species=__)
        for __ in srst2_schema
    ])
    return SRST2Aggregator.from_data(data)


def aggregate_vcfs(wildcards) -> list[str]:
    mapping_tools = config['tools']['mapping'].split('|')

    data_list = []
    for tool in mapping_tools:
        if tool == 'bactmap':
            cj = checkpoints.bactmap_samplesheet
        elif tool.startswith('sarek_'):
            cj = checkpoints.sarek_samplesheet
        elif tool == 'snippy':
            # TODO: fix for cases where no species wildcard present
            # cj = checkpoints.snippy_samplesheet
            raise NotImplementedError
        else:
            raise NotImplementedError

        # TODO: temp
        for species in ['Escherichia_coli']:
            sample_info = read_sample_info(cj, species=species)
            samplesheet = read_samplesheet(cj, species=species)
            data_list.append(
                samplesheet.filter(['sample']).merge(sample_info)
            )

    data = pd.concat(data_list)
    return VCFAggregator.from_data(data)


def _aggregate_over_wcs(checkpoint_job, **wildcards):
    sample_info = read_sample_info(checkpoint_job, **wildcards)
    samplesheet = read_samplesheet(checkpoint_job, **wildcards)
    return samplesheet.filter(['sample']).merge(sample_info)


def aggregate_bactmap(wildcards) -> list[str]:
    return _aggregate_nfcore_mapping(BACTMAP_TEMPLATE, checkpoints.bactmap_samplesheet, **wildcards)


def aggregate_sarek(wildcards) -> list[str]:
    return _aggregate_nfcore_mapping(SAREK_TEMPLATE, checkpoints.sarek_samplesheet, **wildcards)


rule all_bracken:
    input:
        aggregate_bracken


rule all_snippy:
    input:
        aggregate_snippy


rule all_srst2:
    input:
        aggregate_srst2


rule all_vcfs:
    input:
        aggregate_vcfs


# NF-core pipelines

rule all_bactmap:
    input:
        aggregate_bactmap


rule all_taxprofiler:
    input:
        'results/taxprofiler/multiqc/multiqc_report.html'


rule all_sarek:
    input:
        aggregate_sarek


def read_sample_info(rule: Checkpoint, **wildcards) -> pd.DataFrame:
    job: CheckpointJob = rule.get(**wildcards)
    return pd.read_csv(job.rule.input['sample_info'])


def read_samplesheet(rule: Checkpoint, **wildcards) -> pd.DataFrame:
    job: CheckpointJob = rule.get(**wildcards)
    return pd.read_csv(job.output['samplesheet'])


def _aggregate_nfcore_mapping(path_template: str, rule: Checkpoint, **wildcards):
    sample_info = read_sample_info(rule, **wildcards)
    species = sample_info['reference_genome'].unique()
    species = [sp for sp in species if not read_samplesheet(rule, species=sp).empty]
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
        return (
            data
            # TODO: Temp filtering as patch
            [data['id'].str.startswith('B')]
            .dropna(subset=['family'])
        )


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
    """Collect VCF parquet files from all mapping tools output."""

    _vcfs_template = VCF_TEMPLATE

    def __init__(self):
        super().__init__(data)

    @property
    def mapping_tools(self):
        return config['tools']['mapping'].split('|')

    def prepare_templates(self) -> list[str]:
        return [self._vcfs_template.format(mapping_tool=x) for x in self.mapping_tools]

    def prepare_data(self, data) -> pd.DataFrame:
        return (
            data
            # NOTE: Robust to missing data not affecting required fields.
            .dropna(subset=['reference_genome', 'family', 'id', 'library', 'sample'])
            .rename(columns={'reference_genome': 'species'})
        )