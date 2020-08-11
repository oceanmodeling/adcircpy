from adcircpy.outputs.base import _ScalarSurfaceOutputTimeseries


class Fort63(_ScalarSurfaceOutputTimeseries):
    _filetype = "fort.63"
    _cmap = 'jet'
    _levels = 256
