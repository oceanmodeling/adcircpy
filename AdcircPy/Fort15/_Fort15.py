from AdcircPy.Mesh   import Mesh
from AdcircPy.Fort15 import TidalDB

def __init__(self, **kwargs):
    self.header_1         = kwargs.pop("header_1", None)
    self.header_2         = kwargs.pop("header_2", None)
    self.reference_date   = kwargs.pop("reference_date", None)
    self.reference_time   = kwargs.pop("reference_time", None)
    self.RNDAY            = kwargs.pop("RNDAY", None)
    self.DRAMP            = kwargs.pop("DRAMP", None)
    self.nodal_attributes = kwargs.pop("nodal_attributes", None)
    self.constituent_list = kwargs.pop("constituent_list", None)
    self.boundary_TPXO    = kwargs.pop("boundary_TPXO", None)
    self._mode            = kwargs.pop("mode", None)
    self.eta_stations     = kwargs.pop("eta_stations", None)
    self.met_stations     = kwargs.pop("met_stations", None)
    self.vel_stations     = kwargs.pop("vel_stations", None)
    self.eta_start        = kwargs.pop("eta_start", None)
    self.eta_stop         = kwargs.pop("eta_stop", None)
    self.eta_steps        = kwargs.pop("eta_steps", None)
    self.vel_start        = kwargs.pop("vel_start", None)
    self.vel_stop         = kwargs.pop("vel_stop", None)
    self.vel_steps        = kwargs.pop("vel_steps", None)
    self.hotstart_steps   = kwargs.pop("hotstart_steps", None)
    self.output_filepath  = kwargs.pop("output_filepath", None)

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

def parse_fort15(path):
    if path is None:
        return
    else:
        fort15 = dict()
        f = open(path)
        fort15['run_description']=f.readline().split('!')[0].strip()
        fort15['run_id']=f.readline().split('!')[0].strip()
        fort15['error_override_flag'] = int(f.readline().split('!')[0])
        fort15['log_level_flag'] = int(f.readline().split('!')[0])
        fort15['verbose_level_flag'] = int(f.readline().split('!')[0])
        fort15['hotstart_flag'] = int(f.readline().split('!')[0]) # required input
        fort15['coordinate_system_flag'] = int(f.readline().split('!')[0]) # required input
        fort15['model_type_flag'] = int(f.readline().split('!')[0]) # required input
        if fort15['model_type_flag'] == 21:
            fort15['form_of_density_forcing_flag'] = int(f.readline().split('!')[0])
        fort15['bottom_stress_parameterization_flag'] = int(f.readline().split('!')[0])
        fort15['finite_amplitude_flag'] = int(f.readline().split('!')[0])
        fort15['advective_terms_flag'] = int(f.readline().split('!')[0])
        fort15['time_derivative_of_advective_terms_flag'] = int(f.readline().split('!')[0])
        number_of_nodal_attributes = int(f.readline().split('!')[0]) # skip total number of nodal attributes
        if number_of_nodal_attributes > 0:
            fort15['nodal_attributes'] = list()
            for i in range(number_of_nodal_attributes):
                fort15['nodal_attributes'].append(f.readline().split('!')[0].split()[0])
        else:
            fort15['nodal_attributes'] = None
        fort15['coriolis_flag'] = int(f.readline().split('!')[0])
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