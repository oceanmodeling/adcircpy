from abc import abstractmethod
from datetime import datetime

from pyproj import CRS


class WindForcing:
    def __init__(self, start_date: datetime, end_date: datetime, crs: CRS, nws: int = 20):
        self.start_date = start_date
        self.end_date = end_date
        self.NWS = nws
        self.crs = crs

    @property
    def start_date(self):
        return self._start_date

    @start_date.setter
    def start_date(self, start_date: datetime):
        self._start_date = start_date

    @property
    @abstractmethod
    def _start_date(self):
        raise NotImplementedError

    @_start_date.setter
    @abstractmethod
    def _start_date(self, start_date):
        raise NotImplementedError

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self, end_date: datetime):
        self._end_date = end_date

    @property
    @abstractmethod
    def _end_date(self):
        raise NotImplementedError

    @_end_date.setter
    @abstractmethod
    def _end_date(self, end_date):
        raise NotImplementedError

    @property
    def NWS(self):
        try:
            return self.__nws
        except AttributeError:
            return 20

    @NWS.setter
    def NWS(self, nws: int):
        assert nws in [19, 20]
        self.__nws = int(nws)

    def transform_to(self, crs: CRS):
        if self.crs != crs:
            # TODO implement reprojection
            pass
