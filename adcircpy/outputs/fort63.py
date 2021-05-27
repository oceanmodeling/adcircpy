from adcircpy.outputs.base import ScalarSurfaceOutputTimeseries


class Fort63(ScalarSurfaceOutputTimeseries):
    _filetype = 'fort.63'
    _cmap = 'jet'
    _levels = 256
