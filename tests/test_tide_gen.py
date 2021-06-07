from pathlib import Path
from unittest.mock import patch

from adcircpy.cmd import tide_gen
from tests import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'

download_mesh(
    url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
    directory=INPUT_DIRECTORY,
)


def test_tide_gen():
    cmd = [
        'tide_gen',
        f'{INPUT_DIRECTORY / "fort.14"}',
        '2021-02-26T00:00:00',
        '15',
        '--mesh-crs=epsg:4326',
        '--output-file=/dev/null',
    ]

    with patch('sys.argv', cmd):
        tide_gen.main()
