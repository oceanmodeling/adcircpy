from collections import OrderedDict
from datetime import timedelta

class _StationOutputs(OrderedDict):
  def __init__(self, sampling_frequency, netcdf, spinup, **kwargs):
    super(_StationOutputs, self).__init__(**kwargs)
    self.init_sampling_frequency(sampling_frequency)
    self.spinup = spinup
    self.netcdf = netcdf

  def init_sampling_frequency(self, sampling_frequency):
    if sampling_frequency is None:
      self.sampling_frequency=timedelta(minutes=6)
    elif isinstance(sampling_frequency, timedelta):
      self.sampling_frequency=sampling_frequency
    else:
      raise TypeError('sampling_frequency must be a datetime.timedelta object.')