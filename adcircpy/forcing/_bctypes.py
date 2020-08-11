from abc import ABCMeta, abstractmethod


class _BoundaryCondition(metaclass=ABCMeta):
    @property
    @abstractmethod
    def btype(self):
        raise NotImplementedError


class _EtaBc(_BoundaryCondition):
    @property
    def btype(self):
        return 'iettype'


class _VelBc(_BoundaryCondition):
    @property
    def btype(self):
        return 'ifltype'


class _TempBc(_BoundaryCondition):
    @property
    def btype(self):
        return 'itetype'


class _SalBc(_BoundaryCondition):
    @property
    def btype(self):
        return 'isatype'


class _TraceBc(_BoundaryCondition):
    @property
    def btype(self):
        return 'itrtype'
