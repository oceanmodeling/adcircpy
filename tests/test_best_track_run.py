from pathlib import Path
import shutil
from unittest.mock import patch

from adcircpy.cmd import best_track_run
from tests import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'

download_mesh(
        url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
        directory=INPUT_DIRECTORY,
)


def test_best_track_run():
    cmd = [
        'best_track_run',
        f'{INPUT_DIRECTORY / "fort.14"}',
        'Sandy2012',
        '--spinup-days=0.5',
        '--crs=EPSG:4326',
        '--output-directory=/tmp/test',
        '--constituents=all',
        '--overwrite',
        '--timestep=10.',
        '--tau0-gen',
        f'--stations-file={Path(__file__).parent}/stations.txt',
        '--elev-stat=6.',
        '--run-days=0.5',
    ]
    if shutil.which('padcirc') is None:
        if shutil.which('adcirc') is None:
            cmd.append('--nproc=1')
        else:
            cmd.append('--skip-run')

    with patch('sys.argv', cmd):
        best_track_run.main()
