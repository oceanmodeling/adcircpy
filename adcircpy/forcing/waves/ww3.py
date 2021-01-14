from os import PathLike

from adcircpy.forcing.waves import WaveForcing


class WaveWatch3DataForcing(WaveForcing):
    def write(self, directory: PathLike, overwrite: bool = False):
        # ww3data is just a netCDF file so needs no fort.22
        pass

