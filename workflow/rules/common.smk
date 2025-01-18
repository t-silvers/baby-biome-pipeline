import copy
import functools
import os
import pathlib
import subprocess
import tempfile
from typing import Optional

from jinja2 import Template
import yaml


ETL_DIR = '../baby-biome-data-infra'


class DuckDB:
    """Run SQL-based transforms from a template using DuckDB."""

    def __init__(self):
        self._check()

    def __call__(self, path: str, db: Optional[str] = 'memory', readonly: Optional[bool] = False) -> None:
        if isinstance(db, pathlib.Path):
            db = db.as_posix()
        cmd = copy.deepcopy(self._cli)
        if db != 'memory':
            db_parent = pathlib.Path(db).parent
            if not db_parent.exists():
                raise FileNotFoundError(f'Parent directory of db, {db_parent}, does not exist.')
            if readonly:
                cmd += ['-readonly']
            cmd += [db]
        cmd += ['-c ".read', path, '"']
        subprocess.run(' '.join(cmd), shell=True, check=True)

    @classmethod
    def from_template(cls, template: str, params: dict, db: Optional[str] = 'memory', readonly: Optional[bool] = False, log: Optional[str] = None) -> None:
        duckdb = cls()
        sql = duckdb.write_sql(template, params)
        if log:
            with open(log, 'w') as f:
                f.write(sql)
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.sql', delete=True) as f:
            f.write(sql)
            f.flush()
            duckdb(f.name, db=db)

    @property
    def rc(self):
        return config['duckdbrc_slurm']

    @property
    def exe(self):
        return os.environ['DUCKDB']

    @property
    def _cli(self) -> list:
        return [self.exe, '-init', self.rc]

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


def source_etl(d, seed_or_model):
    def func(x):
        return workflow.source_path(ETL_DIR + '/' + seed_or_model + '/' + x)
    if isinstance(d, dict):
        return {k: source_etl(v, seed_or_model) for k, v in d.items()}
    elif isinstance(d, list):
        return [source_etl(v, seed_or_model) for v in d]
    else:
        return func(d)


# Directories
data = pathlib.Path(config['directories']['data'])
resources = pathlib.Path(config['directories']['resources'])


# ETL
etl_cfg = workflow.source_path(ETL_DIR + '/' + 'config.yml')

with open(etl_cfg, 'r') as f:
    data_infra = yaml.safe_load(f)

models, seeds = [source_etl(data_infra[d], d) for d in ['models', 'seeds']]