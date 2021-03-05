from abc import ABC, abstractmethod
from os import PathLike

import numpy as np


class TidalDataset(ABC):
    CONSTITUENTS: [str] = NotImplementedError

    def __init__(self, path: PathLike = None):
        self.path = str(path) if path is not None else None

    @abstractmethod
    def __call__(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> (np.ndarray, np.ndarray):
        raise NotImplementedError

    @abstractmethod
    def get_amplitude(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def get_phase(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError

    @property
    @abstractmethod
    def x(self) -> np.ndarray:
        """ 1xN array of X values """
        raise NotImplementedError

    @property
    @abstractmethod
    def y(self) -> np.ndarray:
        """ 1xN array of Y values """
        raise NotImplementedError

    @staticmethod
    def _assert_vertices(vertices: np.ndarray):
        assert len(vertices.shape) == 2 and vertices.shape[1] == 2, \
            'vertices must be of shape Mx2'
