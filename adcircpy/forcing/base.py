from abc import ABC, abstractmethod
from datetime import timedelta
from os import PathLike


class Forcing(ABC):
    def __init__(self, interval: timedelta):
        self.interval = interval

    @abstractmethod
    def write(self, directory: PathLike, overwrite: bool = False):
        raise NotImplementedError

    @property
    def interval(self) -> timedelta:
        return self.__interval

    @interval.setter
    def interval(self, interval: timedelta):
        if interval is not None and not isinstance(interval, timedelta):
            interval = timedelta(seconds=interval)
        self.__interval = interval

    def __eq__(self, other: 'Forcing') -> bool:
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__
