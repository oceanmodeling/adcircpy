from os import PathLike

from adcircpy.forcing.winds import WindForcing


class AtmosphericMeshForcing(WindForcing):
    def write(self, directory: PathLike, overwrite: bool = False):
        # atmesh is only a netCDF file so no fort.22 is needed
        pass
