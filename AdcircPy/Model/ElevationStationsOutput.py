

def from_fort15(cls, path, netcdf=True):
  stations=dict()
  with open(path,'r') as f:
    for line in f:
      if 'NOUTE' in line:
        print(line)

  return cls(stations, netcdf)

def from_csv(self, path):
  pass

def init_params(self):
  pass
