import numpy as np

class NodalAttributes(dict):
  def __init__(self, spinup_attributes=None, runtime_attributes=None, **attributes):
    super(NodalAttributes, self).__init__(**attributes)
    self.spinup_attributes = spinup_attributes
    self.runtime_attributes = runtime_attributes
    self._init_spinup_attributes()
    self._init_runtime_attributes()

  @classmethod
  def parse_fort13(cls, path, spinup_attributes=None, runtime_attributes=None):
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

  def __check_attributes(self, attributes):
    if isinstance(attributes, list):
      for attribute in attributes:
        if attribute not in self.keys():
          raise IOError('Attribute \'{}\' not found in fort.13.'.format(attribute))

  def _init_spinup_attributes(self):
    self.__check_attributes(self.spinup_attributes)
    if self.spinup_attributes is None:
      self.spinup_attributes = list()
      if 'mannings_n_at_sea_floor' in self.keys():
        self.spinup_attributes.append('mannings_n_at_sea_floor')
      if 'primitive_weighting_in_continuity_equation' in self.keys():
        self.spinup_attributes.append('primitive_weighting_in_continuity_equation')
      if 'surface_submergence_state' in self.keys():
        self.spinup_attributes.append('surface_submergence_state')

  def _init_runtime_attributes(self):
    self.__check_attributes(self.runtime_attributes)
    if self.runtime_attributes is None:
      self.runtime_attributes = list()
      for attribute in self.keys():
        self.runtime_attributes.append(attribute)
