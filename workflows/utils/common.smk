import copy
import os
import json
import pathlib
import subprocess
import tempfile
from typing import Optional

from jinja2 import Template


def check_dir(p: str) -> pathlib.Path:
    d = pathlib.Path(p)
    d.mkdir(parents=True, exist_ok=True)
    return d


data, resources = [check_dir(config['directories'][x]) for x in ['data', 'resources']]


# TODO: Use schema to check for required wildcards
wc_params = {k: v.split('|') for k, v in config['wildcards'].items()}


def parse_data_tools():
    with open(workflow.source_path('../data-infra/config.json'), 'r') as f:
        dcfg = json.load(f)

    return [
        {k: workflow.source_path(f'../data-infra/{x}/{v}') for k, v in dcfg[x].items()}
        for x in ['models', 'seeds']
    ]


models, seeds = parse_data_tools()


# TODO: Add SQL flavor flexibility via sqlglot
class DuckDB:
    """Run SQL-based transforms from a template using DuckDB."""

    def __init__(self):
        self._check()

    def __call__(self, path: str, db: Optional[str] = 'memory', readonly: Optional[bool] = False) -> None:
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