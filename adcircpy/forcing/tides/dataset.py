from abc import ABC, abstractmethod
from os import PathLike
from typing import List, Tuple

import numpy as np


class TidalDataset(ABC):
    def __init__(self, path: PathLike = None):
        """
        create a new tidal dataset object
        :param path: file path or URL pointing to dataset location
        """

        self.path = str(path) if path is not None else None

    def __call__(
        self, constituent: str, vertices: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        get tidal ampltidue and phase
        :param constituent: tidal constituent
        :param vertices: XY locations at which to sample (Mx2)
        :return: amplitude and phase arrays at given locations
        """
        return self.get_amplitude(constituent, vertices), self.get_phase(constituent, vertices)

    @abstractmethod
    def get_amplitude(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
        """
        generate tidal ampltidue
        :param constituent: tidal constituent
        :param vertices: XY locations at which to sample (Mx2)
        :return: amplitude at given locations
        """
        raise NotImplementedError

    @abstractmethod
    def get_phase(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
        """
        generate tidal phase
        :param constituent: tidal constituent
        :param vertices: XY locations at which to sample (Mx2)
        :return: phase at given locations
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def x(self) -> np.ndarray:
        """
        :return: 1D array of X values of vertices
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def y(self) -> np.ndarray:
        """
        :return: 1D array of Y values of vertices
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def constituents(self) -> List[str]:
        """
        :return: list of constituents available on the data source.
        """
        raise NotImplementedError

    @staticmethod
    def _assert_vertices(vertices: np.ndarray):
        """
        :param vertices: list of XY locations
        :return: whether vertices are in XY format (Mx2)
        """
        assert (
            len(vertices.shape) == 2 and vertices.shape[1] == 2
        ), 'vertices must be of shape Mx2'

    def __eq__(self, other: 'TidalDataset') -> bool:
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__
