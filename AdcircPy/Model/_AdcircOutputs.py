from collections import OrderedDict
from datetime import timedelta

class _AdcircOutputs(OrderedDict):
  def __init__(self, sampling_frequency, netcdf, **kwargs):
    super(_AdcircOutputs, self).__init__(**kwargs)
    self.init_sampling_frequency(sampling_frequency)
    self.netcdf = netcdf

  def init_sampling_frequency(self, sampling_frequency):
    if sampling_frequency is None:
      self.sampling_frequency=timedelta(minutes=6)
    elif isinstance(sampling_frequency, timedelta):
      self.sampling_frequency=sampling_frequency
    else:
      raise TypeError('sampling_frequency must be a datetime.timedelta object.')