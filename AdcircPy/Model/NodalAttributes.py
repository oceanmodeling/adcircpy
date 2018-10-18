import numpy as np

class NodalAttributes(dict):
  def __init__(self, spinup_attributes=None, runtime_attributes=None, **attributes):
    super(NodalAttributes, self).__init__(**attributes)
    self.spinup_attributes = spinup_attributes
    self.runtime_attributes = runtime_attributes

  
  @classmethod
  def from_fort13(cls, path, spinup_attributes=None, runtime_attributes=None):
    fort13={}
    with open(path, 'r') as f:
      f.readline().strip()
      NP = int(f.readline().split()[0])
      NAttr = int(f.readline().split()[0]) 
      i=0
      while i < NAttr:
        attribute_name = f.readline().strip()
        units = f.readline().strip()
        if units == '1':
            units='unitless'
        f.readline()
        defaults = [float(x) for x in f.readline().split()]
        fort13[attribute_name] = {'units'    : units,
                                  'defaults' : defaults}
        i+=1
      for i in range(NAttr):
        attribute_name = f.readline().strip()
        numOfNodes = int(f.readline())
        values = np.zeros((NP,len(fort13[attribute_name]['defaults'])))
        values[:]=np.nan
        j=0
        while j < numOfNodes:
            str = f.readline().split()
            node_number = int(str[0])-1
            node_values = [float(x) for x in str[1:]]
            values[node_number,:] = node_values
            j+=1
        values[np.where(np.isnan(values[:,0])),:] = fort13[attribute_name]['defaults']
        fort13[attribute_name]['values'] = values
    return cls(spinup_attributes, runtime_attributes, **fort13)

  @property
  def spinup_attributes(self):
    """ """
    return self._spinup_attributes

  @spinup_attributes.setter
  def spinup_attributes(self, attributes):
    self._check_attributes(attributes)
    if attributes is None:
      if len(self.keys())>0:
        self._init_default_spinup_attributes(attributes)
      else:
        self._spinup_attributes = list()
    else:
      if isinstance(attributes, str):
        self._spinup_attributes = list(attributes)
      elif isinstance(attributes, list):
        self._spinup_attributes = attributes

  @property
  def runtime_attributes(self):
    """ """
    return self._runtime_attributes

  @runtime_attributes.setter
  def runtime_attributes(self, attributes):
    self._check_attributes(attributes)
    if attributes is None:
      if len(self.keys())>0:
        self._init_default_runtime_attributes(attributes)
      else:
        self._runtime_attributes = list()
    else:
      if isinstance(attributes, str):
        self._runtime_attributes = list(attributes)
      elif isinstance(attributes, list):
        self._runtime_attributes = attributes
  
  def _check_attributes(self, attributes):
    if attributes is not None:
      if isinstance(attributes, str):
        attributes=list(attributes)
      if isinstance(attributes, list):
        for attribute in attributes:
          if attribute not in self.keys():
            raise IOError('Attribute \'{}\' not found in fort.13.'.format(attribute))

  def _init_default_spinup_attributes(self, attributes):
    self._spinup_attributes = list()
    if 'mannings_n_at_sea_floor' in self.keys():
      self._spinup_attributes.append('mannings_n_at_sea_floor')
    if 'primitive_weighting_in_continuity_equation' in self.keys():
      self._spinup_attributes.append('primitive_weighting_in_continuity_equation')
    if 'surface_submergence_state' in self.keys():
      self._spinup_attributes.append('surface_submergence_state')

  def _init_default_runtime_attributes(self, attributes):
    self._runtime_attributes = list()
    for attribute in self.keys():
      self._runtime_attributes.append(attribute)
