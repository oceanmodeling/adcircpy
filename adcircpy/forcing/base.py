from abc import ABC, abstractmethod
from os import PathLike


class Forcing(ABC):
    @abstractmethod
    def write(self, directory: PathLike, overwrite: bool = False):
        raise NotImplementedError
