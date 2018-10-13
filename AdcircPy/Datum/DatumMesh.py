import numpy as np
from scipy.interpolate import griddata
from AdcircPy.Datum.VDatum import VDatum
from AdcircPy.Model import AdcircMesh as _AdcircMesh

class DatumMesh(object):
  """
  This class should be used to instantiate UnstructuredGrids that 
  represent Datum 
  """

  def __init__(self, x, y, values, elements, source_hdatum, source_vdatum, target_hdatum, target_vdatum, description):
    self.x = x
    self.y = y
    self.values = values
    self.elements = elements
    self.source_hdatum = source_hdatum
    self.source_vdatum = source_vdatum
    self.target_hdatum = target_hdatum
    self.target_vdatum = target_vdatum
    self.description   = description

  @classmethod
  def build_datum_mesh(cls, AdcircMesh, source_hdatum, target_hdatum, source_vdatum, target_vdatum, vdatum_jar_path):
    if isinstance(AdcircMesh, _AdcircMesh)==False:
      AdcircMesh = _AdcircMesh.from_fort14(AdcircMesh)
    xyz = VDatum.jar_wrapper(AdcircMesh.get_xyz(), source_hdatum, source_vdatum, target_hdatum, target_vdatum, vdatum_jar_path, return_nodata=True)
    x = xyz[:,0]
    y = xyz[:,1]
    values = np.ma.masked_equal(xyz[:,2], -999999.)
    values = AdcircMesh.values - values
    masked, = np.where(values.mask)
    not_masked,= np.where(~values.mask)
    pad_values = griddata((x[not_masked], y[not_masked]), values[not_masked], (x[masked], y[masked]), method='nearest')
    for i, _idx in enumerate(masked):
      values[_idx] = pad_values[i]
    return cls(x, y, values, AdcircMesh.elements, source_hdatum, source_vdatum, target_hdatum, target_vdatum, AdcircMesh.description)

  def dump(self, path):
    with open(path,'w') as f:
      f.write('{:<24}'.format(self.description))
      f.write('{}:{} to {}:{} (add to convert)\n'.format(self.source_hdatum, self.source_vdatum, self.target_hdatum, self.target_vdatum))
      f.write('{:<9d}{:<9d}\n'.format(self.elements.shape[0], self.x.size))
      for i, value in enumerate(self.values):
        f.write('{:>10d}'.format(i+1))
        f.write('{:>16.10f}'.format(self.x[i]))
        f.write('{:>16.10f}'.format(self.y[i]))
        f.write('{:>16.10f}'.format(value))
        f.write('\n')
