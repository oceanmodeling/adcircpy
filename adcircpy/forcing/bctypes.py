import abc


class BoundaryCondition(metaclass=abc.ABCMeta):

    @property
    @abc.abstractproperty
    def btype(self):
        """ """


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
