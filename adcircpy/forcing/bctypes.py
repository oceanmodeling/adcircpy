from abc import ABCMeta, abstractmethod


class BoundaryCondition(metaclass=ABCMeta):
    @property
    @abstractmethod
    def btype(self):
        raise NotImplementedError


class EtaBc(BoundaryCondition):
    @property
    def btype(self):
        return 'iettype'


class VelBc(BoundaryCondition):
    @property
    def btype(self):
        return 'ifltype'


class TempBc(BoundaryCondition):
    @property
    def btype(self):
        return 'itetype'


class SalBc(BoundaryCondition):
    @property
    def btype(self):
        return 'isatype'


class TraceBc(BoundaryCondition):
    @property
    def btype(self):
        return 'itrtype'
