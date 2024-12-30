import os
import json
import pathlib
import subprocess
import tempfile
from typing import Optional

from jinja2 import Template


# directories
run_id = str(config['run_id'])

outdir = pathlib.Path(config['directories']['out'])
data = outdir / 'data'
resources = outdir / 'resources'
results = outdir / 'results' / run_id

logdir = pathlib.Path(config['directories']['log'])
logdir = logdir / run_id

# wildcards
wc_params = {k: v.split('|') for k, v in config['wildcards'].items()}
# TODO: Check for required wildcards

# executor
max_submit = int(os.environ['SLURM_MAX_SUBMIT_JOB_LIMIT'])

# etl config
data_infra_path = workflow.source_path('../../config/data_infra.json')
with open(data_infra_path, 'r') as f:
    data_infra_cfg = json.load(f)

models = data_infra_cfg['models']
seeds = data_infra_cfg['seeds']


# TODO: Add SQL flavor flexibility via sqlglot
class DuckDB:
    """Run SQL-based transforms from a template using DuckDB."""

    def __init__(self):
        self._check()

    def __call__(self, path: str, db: Optional[str] = 'memory') -> None:
        if db == 'memory':
            cmd = [self.exe, '-init', self.rc, '-c ".read', path, '"']
        else:
            db_parent = pathlib.PurePosixPath(db).parent
            if not db_parent.exists():
                raise FileNotFoundError(
                    f'Parent directory of db, {db_parent}, does not exist.'
                )
            cmd = [self.exe, '-init', self.rc, db, '-c ".read', path, '"']
        subprocess.run(' '.join(cmd), shell=True, check=True)

    @classmethod
    def from_template(cls, template: str, params: dict, db: Optional[str] = 'memory', log: Optional[str] = None) -> None:
        duckdb = cls()
        sql = duckdb.write_sql(template, params)
        if log:
            with open(log, 'w') as f:
                f.write(sql)
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.sql', delete=True) as f:
            f.write(sql)
            f.flush()
            duckdb(f.name)

    @property
    def rc(self):
        return workflow.source_path('../../config/duckdbrc-slurm')

    @property
    def exe(self):
        return os.environ['DUCKDB']

    @staticmethod
    def write_sql(template, params):
        with open(template) as f:
            t = Template(f.read())
        return t.render(params)

    def _check(self):
        if 'DUCKDB' not in os.environ:
            raise EnvironmentError('DUCKDB environment variable not set.')

        if self.exe is None:
            raise FileNotFoundError('DuckDB executable not found.')

        if self.rc is None:
            raise FileNotFoundError('DuckDB rc file not found.')


transform = DuckDB.from_template
"""Convenience function for running dynamically generated SQL scripts for ETL using DuckDB."""


# TODO: Deprecate
def covars_to_partitions(df, covars=None):
    """Form hive-partitioned path substring from df columns."""
    EXCLUDE = ['species', 'sample']

    if 'partitions' in df.columns:
        raise ValueError

    if covars is None:
        covars = [col for col in df.columns if col not in EXCLUDE]
    
    def row_to_partitions(row):
        return '/'.join([f'{k}={v}' for k, v in row.items()])

    return df.filter(covars).apply(row_to_partitions, axis=1)


def data_path_from_template(template: str, params: dict) -> str:
    return path_from_template(data, template, params)


def path_from_template(directory: pathlib.Path, template: str, params: dict) -> str:
    return (directory / template.format(**params)).as_posix()