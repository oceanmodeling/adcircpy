from adcircpy.outputs.base import _ScalarSurfaceOutputTimeseries


class Fort63(_ScalarSurfaceOutputTimeseries):
    _filetype = "fort.63"
    _cmap = 'jet'
    _levels = 256

    def subset(self, dst, overwrite=False):
        zeta = self._ptr['zeta'][:].data
        zeta_subset = zeta[:, idx_subset]
        zout = dst.createVariable('zeta', 'f8', ('time', 'node'))
        zout[:] = zeta_subset[:]
