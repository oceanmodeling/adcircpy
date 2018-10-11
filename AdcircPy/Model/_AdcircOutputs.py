from datetime import timedelta

def init_sampling_frequency(self, sampling_frequency):
  if sampling_frequency is None:
    self.sampling_frequency=timedelta(minutes=6)
  elif isinstance(sampling_frequency, timedelta):
    self.sampling_frequency=sampling_frequency
  else:
    raise TypeError('sampling_frequency must be a datetime.timedelta object.')