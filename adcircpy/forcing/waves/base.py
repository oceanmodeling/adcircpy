from os import PathLike

from adcircpy.forcing.base import Forcing


class WaveForcing(Forcing):
    def __init__(self, nrs: int, interval_seconds: int):
        self.NRS = nrs
        self.RSTIMINC = interval_seconds

    def write(self, directory: PathLike, overwrite: bool = False):
        # TODO implement this
        raise NotImplementedError(
            'writing wave forcing to file is not yet implemented')
