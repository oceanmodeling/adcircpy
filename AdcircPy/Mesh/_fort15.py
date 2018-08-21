from collections import OrderedDict
from datetime import datetime
import numpy as np
from AdcircPy import Mesh
# from AdcircPy.core import TidalDB

def _init_fort15(self, path):
    if path is None:
        return
    _f = open(path)
    def __line():
        return _f.readline().split('!')[0].strip()
    f = dict()
    f['RUNDES'] = __line()
    f['RUNID'] = __line()
    f['NFOVER'] = int(__line())
    f['NABOUT'] = int(__line())
    f['NSCREEN'] = int(__line())
    f['IHOT'] = int(__line()) # required input
    f['ICS'] = int(__line()) # required input
    f['IM'] = int(__line()) # required input
    if f['IM'] == 21:
        f['IDEN'] = int(__line())
    f['NOLIBF'] = int(__line())
    f['NOLIFA'] = int(__line())
    f['NOLICA'] = int(__line())
    f['NOLICAT'] = int(__line())
    number_of_nodal_attributes = int(__line()) # skip total number of nodal attributes
    if number_of_nodal_attributes > 0:
        f['AttrName'] = list()
        for i in range(number_of_nodal_attributes):
            f['AttrName'].append(__line().split()[0])
    f['NCOR'] = int(__line())
    f['NTIP'] = int(__line())
    f['NWS'] = int(__line())
    f['NRAMP'] = int(__line())
    f['G'] = float(__line())
    f['TAU0'] = int(__line())
    # Untested!
    if f['TAU0'] == -5:
        line = _f.readline().split('!')
        f['Tau0FullDomainMin'] = float(line[0])
        f['Tau0FullDomainMax'] = float(line[1])
    f['DTDP'] = float(__line())
    f['STATIM'] = float(__line())
    f['REFTIM'] = float(__line())
    if f['NWS'] > 0:
        if f['NWS']==11 or f['NWS']==9:
            raise Exception('Invalid NWS number.')
        f['WTIMNC'] = _f.readline()
    f['RNDAY'] = float(__line())
    f['DRAMP'] = _f.readline()
    f['A00, B00, C00'] = _f.readline()
    if f['NOLIFA'] in [0,1]:
        f['H0'] = _f.readline()
    elif f['NOLIFA'] in [2,3]:
        f['H0, INTEGER, INTEGER, VELMIN'] = _f.readline()
    f['SLAM0, SFEA0'] = _f.readline()
    if f['NOLIBF']==0:
        f['TAU'] = _f.readline()
    elif f['NOLIBF']==1:
        f['CF'] = _f.readline()
    elif f['NOLIBF']==2:
        f['CF, HBREAK, FTHETA, FGAMMA'] = _f.readline()
    if f['IM'] in [0,1,2]:
        f['ESLM'] = _f.readline()
    elif f['IM'] == 10:
        f['ESLM, ESLC'] = _f.readline()
    f['CORI'] = float(__line())
    f['NTIF'] = int(__line())
    _constituents = list()
    for i in range(f['NTIF']):
        _constituents.append(__line()); __line()
    f['NTIF'] = {'NTIF'         : f['NTIF'],
                 'constituents' : _constituents}
    _NBFR = int(__line())
    _parameters=list()
    for i in range(_NBFR):
        _parameters.append((__line(), _f.readline()))
    f['BOUNTAG'] = OrderedDict()
    for _constituent, _constants in _parameters:
        _c = _constants.split()
        f['BOUNTAG'][_constituent] = {
            'forcing_frequency'     : float(_c[0]),
            'nodal_factor'          : float(_c[1]),
            'equilibrium_arguments' : float(_c[2]),
            'amplitude'             : list(),
            'phase'                 : list()}

    for _constituent in list(f['BOUNTAG'].keys()):
        __line()
        for _boundary in self.Boundaries.ocean_boundaries:
            for j in range(len(_boundary)):
                line = _f.readline().split()
                f['BOUNTAG'][_constituent]['amplitude'].append(float(line[0]))
                f['BOUNTAG'][_constituent]['phase'].append(float(line[1]))
        f['BOUNTAG'][_constituent]['amplitude'] = np.asarray(f['BOUNTAG'][_constituent]['amplitude'])
        f['BOUNTAG'][_constituent]['phase'] = np.asarray(f['BOUNTAG'][_constituent]['phase'])
    f['ANGINN'] = float(__line())
    if self.Boundaries.inflow_boundaries is not None:
        raise NotImplementedError('This fort.15 appears to include inflow boundaries which are not yet supported')
    f['StationOutputs'] = dict()
    def __stations():
        line = __line().split()
        _dict = { 'id'         : list(),
                  'coords'     : list(),
                  'format'     : int(line[0]),
                  'start_days' : float(line[1]),
                  'stop_days'  : float(line[2]),
                  'interval'   : float(line[3])}
        for i in range(int(__line())):
            line = _f.readline().split('!')
            if len(line)>1:
                _id = line[-1].strip()
            else:
                _id = i
            _dict['id'].append(_id)
            line = line[0].split()
            _dict['coords'].append((float(line[0]), float(line[1])))
        return _dict
    f['StationOutputs']['elevation'] = __stations()
    f['StationOutputs']['velocity'] = __stations()
    if f['IM'] == 10:
        f['StationOutputs']['NSTAC'] = __stations()
    if f['NWS'] != 0:
        f['StationOutputs']['meteo'] = __stations()
    f['NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE'] = _f.readline()
    f['NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV'] = _f.readline()
    if f['IM'] == 10:
        f['NOUTGC, TOUTSGC, TOUTFGC, NSPOOLGC'] = _f.readline()
    if f['NWS'] != 0:
        f['NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW'] = _f.readline()
    _NFREQ = int(__line())
    if _NFREQ > 0:
        for i in range(_NFREQ):
            f['StationOutputs']['HarmonicConstituents'] = _f.readline()
            _f.readline()
    f['THAS, THAF, NHAINC, FMV'] = _f.readline()
    f['NHASE, NHASV, NHAGE, NHAGV'] = _f.readline()
    f['NHSTAR, NHSINC'] = _f.readline()
    f['ITITER, ISLDIA, CONVCR, ITMAX'] = _f.readline()
    if f['IM']>1:
        raise Exception('3D not yet supported.')
    f['NCPROJ'] = _f.readline()
    f['NCINST'] = _f.readline()
    f['NCSOUR'] = _f.readline()
    f['NCHIST'] = _f.readline()
    f['NCREF']  = _f.readline()
    f['NCCOM']  = _f.readline()
    f['NCHOST'] = _f.readline()
    f['NCCONV'] = _f.readline()
    f['NCCONT'] = _f.readline()
    f['NCDATE'] = datetime.strptime(_f.readline().strip(), '%Y-%m-%d %H:%M:%S UTC')
    self.fort15 = Mesh.fort15(**f)

def _generate_forcing_from_TPXO(self):
    boundary_TPXO = list()
    if self.ocean_boundaries is not None:
        Tpxo = TidalDB.TPXO()
        for boundary in self.ocean_boundaries:
            boundary_TPXO.append(Tpxo.get_constituents_at_lonlat(self.x[boundary], self.y[boundary], self.constituent_list))
        self.fort15.boundary_TPXO = boundary_TPXO

def _generate_equilibrium_arguments(self, start_date, end_date):
    pass
