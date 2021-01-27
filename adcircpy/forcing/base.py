from abc import ABC, abstractmethod
from datetime import timedelta
from os import PathLike


class Interval:
    def __init__(self):
        self.value = None

    def __set__(self, instance, value: int):
        if not isinstance(value, timedelta):
            value = timedelta(seconds=value)
        self.value = value

    def __get__(self, instance, owner):
        return self.value


class Forcing(ABC):
    def __init__(self, interval: timedelta):
        self.interval = Interval()
        self.interval = interval

    @abstractmethod
    def write(self, directory: PathLike, overwrite: bool = False):
        raise NotImplementedError
