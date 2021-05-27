from os import PathLike
from pathlib import Path

from adcircpy.forcing.waves import WaveForcing


class WaveWatch3DataForcing(WaveForcing):
    def __init__(self, filename: PathLike, nrs: int = 5, interval_seconds: int = 3600):
        if not isinstance(filename, Path):
            filename = Path(filename)
        self.filename = filename
        super().__init__(nrs=nrs, interval_seconds=interval_seconds)

    def write(self, directory: PathLike, overwrite: bool = False):
        # ww3data is just a netCDF file so needs no fort.22
        pass
