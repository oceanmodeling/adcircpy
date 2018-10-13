from datetime import timedelta
from AdcircPy.Model._ModelOutputs import _ModelOutputs

class ElevationGlobalOutput(_ModelOutputs):
  def __init__(self, sampling_frequency=timedelta(minutes=6), netcdf=True, spinup=False, harmonic_analysis=False):
    super(ElevationGlobalOutput, self).__init__(sampling_frequency, netcdf, spinup, harmonic_analysis)
