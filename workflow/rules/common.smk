import os

group_home = os.environ['GROUP_HOME']

bactmap_config = group_home + '/config/nextflow/bactmap.config'

taxprofiler_cfg = group_home + '/config/nextflow/taxprofiler.config'

archive = '/r/' + os.environ['USER'][0] + '/' + os.environ['USER']

archive_dir = archive + '/20241202_wv'

isolated_species = config['wildcards']['species'].split('|')