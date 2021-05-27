from abc import ABC, abstractmethod
from os import PathLike

from adcircpy.forcing.base import Forcing


class WaveForcing(Forcing, ABC):
    def __init__(self, nrs: int, interval_seconds: int):
        super().__init__(interval_seconds)
        self.NRS = nrs

    @abstractmethod
    def write(self, directory: PathLike, overwrite: bool = False):
        raise NotImplementedError
