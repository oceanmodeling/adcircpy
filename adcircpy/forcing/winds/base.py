from os import PathLike

from adcircpy.forcing.base import Forcing


class WindForcing(Forcing):
    def __init__(self, nws: int, interval_seconds: int):
        self.NWS = nws
        self.WTIMINC = interval_seconds

    def write(self, directory: PathLike, overwrite: bool):
        # TODO implement this
        raise NotImplementedError
