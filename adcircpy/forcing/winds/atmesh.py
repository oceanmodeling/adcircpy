from os import PathLike
from pathlib import Path

from adcircpy.forcing.winds.base import WindForcing


class AtmosphericMeshForcing(WindForcing):
    def __init__(self, filename: PathLike, nws: int = 5, interval_seconds: int = 3600):
        if not isinstance(filename, Path):
            filename = Path(filename)
        self.filename = filename
        super().__init__(nws=nws, interval_seconds=interval_seconds)

    def write(self, directory: PathLike, overwrite: bool = False):
        # atmesh is only a netCDF file so no fort.22 is needed
        pass
