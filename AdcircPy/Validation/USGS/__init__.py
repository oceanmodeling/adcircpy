from collections import OrderedDict
from AdcircPy.Validation.USGS import _HighWaterMarks

class HighWaterMarks(OrderedDict):
  url = 'https://stn.wim.usgs.gov/STNServices/HWMs/FilteredHWMs.json'
  params = dict()
  params['EventType']   = 2 # 2 for hurricane
  params['EventStatus'] = 0 # 0 for completed

  def __call__(self, **kwargs):
    dict.__init__(self, **kwargs)
  
  @classmethod
  def from_event_name(cls, eventName, target_datum, filter_dict=None):
    return _HighWaterMarks.from_event_name(cls, eventName, target_datum, filter_dict)

  @classmethod
  def from_csv(cls, path):
    return _HighWaterMarks.from_csv(cls, path)

  @classmethod
  def get_event_list(cls):
    return _HighWaterMarks.get_event_list(cls)

  @classmethod
  def get_event_id_from_name(cls, eventName):
    return _HighWaterMarks.get_event_id_from_name(cls, eventName)

  def _set_vertical_datum(self, target_datum: str):
    """ Allowed datums can be found on VDatum's website. """
    _HighWaterMarks._set_vertical_datum(self, target_datum)

  def _filter(self, excellent=False, good=False, fair=False, poor=False, riverine=False, non_still_water=False, keep_undefined=False):
    return _HighWaterMarks._filter(self, excellent, good, fair, poor, riverine, non_still_water, keep_undefined)
 
  def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None,**kwargs):
    return _HighWaterMarks.make_plot(self, axes, vmin, vmax, extent, epsg, **kwargs)