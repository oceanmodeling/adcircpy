from functools import lru_cache

# from collections.abc import Mapping
import pathlib

from adcircpy.outputs.maxele import Maxele


class OutputCollection(
    # Mapping
):
    def __init__(
        self,
        fort61=None,
        fort62=None,
        fort63=None,
        fort64=None,
        maxele=None,
        maxvel=None,
        crs=None,
    ):
        self._crs = crs
        self._maxele = maxele

    def __iter__(self):
        for name in self.get_output_types():
            yield self._container[name]

    def __len__(self):
        return len(self.get_output_types())

    def get_output(self, name):
        return self._container[name]

    def get_output_types(self):
        return [_ for _ in self._container.keys() if _ is not None]

    def _certify_output_type(self, inst, obj):
        # TODO: should use _filetype attribute instead or create a Enum class
        if isinstance(inst, obj):
            return inst
        elif isinstance(inst, (str, pathlib.Path)):
            return obj(inst, crs=self.crs)

    @property
    def maxele(self):
        return self._maxele

    @property
    def crs(self):
        return self._crs

    @property
    @lru_cache(maxsize=None)
    def _container(self):
        return {}

    @property
    def _maxele(self):
        return self._container['maxele']

    @property
    def _crs(self):
        return self.__crs

    @_maxele.setter
    def _maxele(self, maxele):
        self._container['maxele'] = self._certify_output_type(maxele, Maxele)

    @_crs.setter
    def _crs(self, crs):
        self.__crs = crs
