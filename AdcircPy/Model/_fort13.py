import numpy as np

def _init_fort13(self, path):
  if path == None:
    self.fort13=None
  else:
    self.fort13={}
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
        self.fort13[attribute_name] = {'units'    : units,
                                  'defaults' : defaults}
        i+=1
      for i in range(NAttr):
        attribute_name = f.readline().strip()
        numOfNodes = int(f.readline())
        values = np.zeros((NP,len(self.fort13[attribute_name]['defaults'])))
        values[:]=np.nan
        j=0
        while j < numOfNodes:
            str = f.readline().split()
            node_number = int(str[0])-1
            node_values = [float(x) for x in str[1:]]
            values[node_number,:] = node_values
            j+=1
        values[np.where(np.isnan(values[:,0])),:] = self.fort13[attribute_name]['defaults']
        self.fort13[attribute_name]['values'] = values

    