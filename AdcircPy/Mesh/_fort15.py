from AdcircPy.Mesh   import Mesh
from AdcircPy import TidalDB

def from_fort15(path):
    fort15 = dict()
    f = open(path)
    fort15['RUNDES']=f.readline().split('!')[0].strip()
    fort15['RUNID']=f.readline().split('!')[0].strip()
    fort15['NFOVER'] = int(f.readline().split('!')[0])
    fort15['NABOUT'] = int(f.readline().split('!')[0])
    fort15['NSCREEN'] = int(f.readline().split('!')[0])
    fort15['IHOT'] = int(f.readline().split('!')[0]) # required input
    fort15['ICS'] = int(f.readline().split('!')[0]) # required input
    fort15['IM'] = int(f.readline().split('!')[0]) # required input
    if fort15['IM'] == 21:
        fort15['IDEN'] = int(f.readline().split('!')[0])
    fort15['NOLIBF'] = int(f.readline().split('!')[0])
    fort15['NOLIFA'] = int(f.readline().split('!')[0])
    fort15['NOLICA'] = int(f.readline().split('!')[0])
    fort15['NOLICAT'] = int(f.readline().split('!')[0])
    number_of_nodal_attributes = int(f.readline().split('!')[0]) # skip total number of nodal attributes
    if number_of_nodal_attributes > 0:
        fort15['AttrName'] = list()
        for i in range(number_of_nodal_attributes):
            fort15['AttrName'].append(f.readline().split('!')[0].split()[0])
    else:
        fort15['AttrName'] = None
    fort15['NCOR'] = int(f.readline().split('!')[0])
    fort15['tidal_potential_and_self_attraction_flag'] = int(f.readline().split('!')[0])
    fort15['wind_forcing_type_flag'] = int(f.readline().split('!')[0])
    fort15['ramping_flag'] = int(f.readline().split('!')[0])
    fort15['gravitational_constant'] = float(f.readline().split('!')[0])
    fort15['GWCE_weighting_factor_flag'] = int(f.readline().split('!')[0])
    # Untested!
    if fort15['GWCE_weighting_factor_flag'] == -5:
        line = f.readline().split('!')
        fort15['Tau0FullDomainMin'] = float(line[0])
        fort15['Tau0FullDomainMax'] = float(line[1])
    fort15['timestep'] = float(f.readline().split('!')[0])
    fort15['offset_days'] = float(f.readline().split('!')[0])
    fort15['reference_time'] = float(f.readline().split('!')[0])

    if fort15['wind_forcing_type_flag'] > 1 and \
        (fort15['wind_forcing_type_flag']!=11 or fort15['wind_forcing_type_flag']!=9):
        
        fort15['meteorological_parameters'] = dict()
        
        if np.abs(fort15['wind_forcing_type_flag']) in [101, 301, 401, 111, 311, 411]:
            fort15['meteorological_parameters']['time_interval_rad_stress'] = float(f.readline().split('!')[0])
        
        elif np.abs(fort15['wind_forcing_type_flag']) in [2, 4, 5, 7, 10, 12, 15, 16]:
            fort15['meteorological_parameters']['time_interval_met_forcing'] = float(f.readline().split('!')[0])
        
        elif np.abs(fort15['wind_forcing_type_flag']) in [102, 302, 402, 104, 304, 404, 105, 305, 405, 107, 307, 407, 110, 310, 410, 112, 312, 412, 115, 315, 415, 116, 316, 416]:
            line = f.readline().split('!')[0].split()
            fort15['meteorological_parameters']['time_interval_met_forcing'] = float(line[0])
            fort15['meteorological_parameters']['time_interval_rad_stress'] = float(line[1])
        
        elif np.abs(fort15['wind_forcing_type_flag']) in [3]:
            line = f.readline().split('!')[0].split()
            fort15['meteorological_parameters']['time_interval_met_forcing'] = float(line[0])
            fort15['meteorological_parameters']['time_interval_rad_stress'] = float(line[1])
    return fort15




def generate_forcing_from_TPXO(self, Mesh):
    
    if isinstance(Mesh, ("".__class__, u"".__class__)):
        fname = Mesh
        Mesh  = Mesh.init_from_file(fname)
    boundary_TPXO = list()
    if Mesh.ocean_boundaries is not None:
        Tpxo = TidalDB.TPXO()
        for boundary in Mesh.ocean_boundaries:
            boundary_TPXO.append(Tpxo.get_constituents_at_lonlat(Mesh.x[boundary], Mesh.y[boundary], self.constituent_list))
        self.boundary_TPXO = boundary_TPXO
    return self

def generate_equilibrium_arguments(self, start_date, end_date):
    pass
