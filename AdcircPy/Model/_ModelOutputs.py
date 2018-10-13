from datetime import timedelta

class _ModelOutputs(object):
  def __init__(self, sampling_frequency, netcdf, spinup, harmonic_analysis):
    self.spinup = spinup
    self.netcdf = netcdf
    self.sampling_frequency = sampling_frequency
    self.harmonic_analysis=harmonic_analysis
    self._init_sampling_frequency()

  def _init_sampling_frequency(self):
    if self.sampling_frequency is None:
      self.sampling_frequency=timedelta(minutes=6)
    elif isinstance(self.sampling_frequency, timedelta)==False:
      raise TypeError('sampling_frequency must be a datetime.timedelta object.')