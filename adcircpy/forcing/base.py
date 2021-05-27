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
        if not isinstance(interval, timedelta):
            interval = timedelta(seconds=interval)
        self.__interval = interval
