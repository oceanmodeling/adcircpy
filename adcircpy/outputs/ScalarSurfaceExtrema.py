# global imports
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# local imports
from AdcircPy.Mesh import UnstructuredMesh


class ScalarSurfaceExtrema(UnstructuredMesh):
    """
    Subclass representing a general Output Surface Extrema, which is
    instantiated by the _OutputFactory class and normally used for ASCII files
    in which the actual output tye is unclear.
    """
    def __init__(self, xy, values, elements, SpatialReference, vertical_datum,
                 model_times=None, node_id=None, element_id=None,
                 **Boundaries):
        UnstructuredMesh.__init__(self, xy, values,  elements, SpatialReference,
                                  node_id=node_id, element_id=element_id,
                                  vertical_datum=vertical_datum, **Boundaries)
        self._model_times = model_times

    def make_plot(self, title='Surface Extrema', axes=None, timesteps=False,
                  vmin=None, vmax=None, **kwargs):
        if axes is None:
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        else:
            ax1 = axes
        super(ScalarSurfaceExtrema, self).make_plot(
            self, axes=ax1, title=title, vmin=vmin, vmax=vmax, **kwargs)
        if len(self._model_times) > 0 and timesteps is True:
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            _ax = ax2.tripcolor(self._Tri, self._times)
            ax2.set_title('Time in seconds after coldstart after which the '
                          + 'extrema happened.')
            ax2.axis('scaled')
            plt.colorbar(_ax)

    @staticmethod
    def is_netcdf(path):
        try:
            Dataset(path)
            return True
        except OSError:
            return False

    @property
    def values(self):
        return self._z

    @property
    def model_times(self):
        return self._model_times

    @property
    def _model_times(self):
        return self.__model_times

    @_model_times.setter
    def _model_times(self, model_times):
        model_times = np.asarray(model_times)
        assert model_times.shape == self.values.shape
        self.__model_times = model_times
