from datetime import timedelta
from AdcircPy.Model._GlobalOutputs import _GlobalOutputs


class ElevationGlobalOutput(_GlobalOutputs):
    def __init__(self, sampling_frequency=timedelta(0), spinup=False,
                 harmonic_analysis=False, netcdf=True):
        super(ElevationGlobalOutput, self).__init__(sampling_frequency, spinup,
                                                    harmonic_analysis, netcdf)
