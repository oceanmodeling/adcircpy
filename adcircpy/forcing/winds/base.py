from abc import ABC, abstractmethod
from os import PathLike

from adcircpy.forcing.base import Forcing


class WindForcing(Forcing, ABC):
    def __init__(self, nws: int, interval_seconds: int):
        self.NWS = nws
        self.WTIMINC = interval_seconds

    @abstractmethod
    def write(self, directory: PathLike, overwrite: bool = False):
        raise NotImplementedError
