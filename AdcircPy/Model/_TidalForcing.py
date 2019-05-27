# global imports
from datetime import datetime, timedelta
import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset
from requests.structures import CaseInsensitiveDict

# local imports
from AdcircPy.utils import __init_TPXO_cache

# unittest imports
import unittest
import os
from AdcircPy import Model

if os.getenv('DOCKERENV') is None:
    tpxo_path = __init_TPXO_cache()
else:
    tpxo_path = None


class _TidalForcing(object):
    """
    tide_fac13.f reimplementation with TPXO included.
    """
    __orbital_frequency = {'M4':      0.0002810378050173,
                           'M6':      0.0004215567080107,
                           'MK3':     0.0002134400613513,
                           'S4':      0.0002908882086657,
                           'MN4':     0.0002783986019952,
                           'S6':      0.0004363323129986,
                           'M3':      0.0002107783537630,
                           '2MK3':    0.0002081166466594,
                           'M8':      0.0005620756090649,
                           'MS4':     0.0002859630068415,
                           'M2':      0.0001405189025086,
                           'S2':      0.0001454441043329,
                           'N2':      0.0001378796994865,
                           'Nu2':     0.0001382329037065,
                           'MU2':     0.0001355937006844,
                           '2N2':     0.0001352404964644,
                           'lambda2': 0.0001428049013108,
                           'T2':      0.0001452450073529,
                           'R2':      0.0001456432013128,
                           '2SM2':    0.0001503693061571,
                           'L2':      0.0001431581055307,
                           'K2':      0.0001458423172006,
                           'K1':      0.0000729211583579,
                           'O1':      0.0000675977441508,
                           'OO1':     0.0000782445730498,
                           'S1':      0.0000727220521664,
                           'M1':      0.0000702594512543,
                           'J1':      0.0000755603613800,
                           'RHO':     0.0000653117453487,
                           'Q1':      0.0000649585411287,
                           '2Q1':     0.0000623193381066,
                           'P1':      0.0000725229459750,
                           'Mm':      0.0000026392030221,
                           'Ssa':     0.0000003982128677,
                           'Sa':      0.0000001991061914,
                           'Msf':     0.0000049252018242,
                           'Mf':      0.0000053234146919}
    __orbital_frequency_units = 'rad/sec'
    __doodson_coefficient = {'M4':      4,
                             'M6':      6,
                             'MK3':     3,
                             'S4':      4,
                             'MN4':     4,
                             'S6':      6,
                             'M3':      3,
                             '2MK3':    3,
                             'M8':      8,
                             'MS4':     4,
                             'M2':      2,
                             'S2':      2,
                             'N2':      2,
                             'Nu2':     2,
                             'MU2':     2,
                             '2N2':     2,
                             'lambda2': 2,
                             'T2':      2,
                             'R2':      2,
                             '2SM2':    2,
                             'L2':      2,
                             'K2':      2,
                             'K1':      1,
                             'O1':      1,
                             'OO1':     1,
                             'S1':      1,
                             'M1':      1,
                             'J1':      1,
                             'RHO':     1,
                             'Q1':      1,
                             '2Q1':     1,
                             'P1':      1,
                             'Mm':      0,
                             'Ssa':     0,
                             'Sa':      0,
                             'Msf':     0,
                             'Mf':      0}
    __tidal_potential_amplitude = {'M2': 0.242334,
                                   'S2': 0.112841,
                                   'N2': 0.046398,
                                   'K2': 0.030704,
                                   'K1': 0.141565,
                                   'O1': 0.100514,
                                   'P1': 0.046843,
                                   'Q1': 0.019256}
    __earth_tidal_potential = {'M2': 0.693,
                               'S2': 0.693,
                               'N2': 0.693,
                               'K2': 0.693,
                               'K1': 0.736,
                               'O1': 0.695,
                               'P1': 0.706,
                               'Q1': 0.695}

    if tpxo_path is not None:
        __nc = Dataset(tpxo_path)
    else:
        __nc = None

    __constituents = set()
    __major8 = ['K1', 'O1', 'P1', 'Q1', 'M2', 'S2', 'N2', 'K2']

    def __init__(self, start_date, end_date, spinup_days=7.,
                 constituents='all'):
        self.__set_start_date(start_date)
        self.__set_spinup_days(spinup_days)
        self.__set_forcing_start_date()
        self.__set_end_date(end_date)
        self.__set_forcing_end_date()
        self.__set_tpxo_constituents()
        self.__set_constituents(constituents)
        self.__init_orbital_parameters()
        self.__set_nodal_factor()
        self.__set_greenwich_term()
        self.__set_tpxo_constituents()

    def __call__(self, constituent):
        return self.get_forcing_factors(constituent)

    def get_forcing_factors(self, constituent):
        try:
            tidal_potential_amplitude \
                = self.tidal_potential_amplitude[constituent]
        except KeyError:
            tidal_potential_amplitude = None
        try:
            earth_tidal_potential \
                = self.earth_tidal_potential[constituent]
        except KeyError:
            earth_tidal_potential = None

        return [tidal_potential_amplitude,
                self.orbital_frequency[constituent],
                earth_tidal_potential,
                self.nodal_factor[constituent],
                self.greenwich_term[constituent]]

    def _set_TPXO_interp(self, AdcircMesh, method='nearest'):
        """
        Method to generate the TPXO interpolation on the open boundaries.
        """
        assert isinstance(AdcircMesh, Model.AdcircMesh)
        if len(AdcircMesh.OceanBoundaries) == 0:
            raise RuntimeError('Cannot generate TPXO parameters for mesh '
                               'without defined ocean boundaries.')
        x = self.nc['lon_z'][:].reshape(self.nc['lon_z'].size)
        y = self.nc['lat_z'][:].reshape(self.nc['lat_z'].size)
        data = list()
        for ocean_boundary in AdcircMesh.OceanBoundaries:
            SpatialReference = ocean_boundary['SpatialReference']
            if not SpatialReference.IsGeographic():
                raise NotImplementedError(
                    'Tidal Run only supported for meshes in '
                    + ' geographic coordinates.')
            _x = ocean_boundary['vertices'][:, 0]
            _x = np.asarray([_ + 360. for _ in _x if _ < 0])
            _y = ocean_boundary['vertices'][:, 1]
            _idx = np.where(np.logical_and(
                            np.logical_and(np.min(_x) <= x, np.max(_x) >= x),
                            np.logical_and(np.min(_y) <= y, np.max(_y) >= y)))
            constituents = CaseInsensitiveDict()
            for constituent in self.constituents:
                if constituent in self.tpxo_constituents.keys():
                    constituents[constituent] = dict()
                    idx = self.tpxo_constituents[constituent]
                    ha = self.nc['ha'][idx, :, :].reshape(
                        self.nc['ha'][idx, :, :].size)
                    hp = self.nc['hp'][idx, :, :].reshape(
                        self.nc['hp'][idx, :, :].size)
                    ha = griddata((x[_idx], y[_idx]), ha[_idx], (_x, _y),
                                  method=method)
                    hp = griddata((x[_idx], y[_idx]), hp[_idx], (_x, _y),
                                  method=method)
                    constituents[constituent] = np.vstack([ha, hp]).T
            data.append(constituents)
        self.__TPXO_interp = data

    def __set_start_date(self, start_date):
        assert isinstance(start_date, datetime)
        self.__start_date = start_date

    def __set_spinup_days(self, spinup_days):
        spinup_days = float(spinup_days)
        spinup_days = np.abs(spinup_days)
        assert spinup_days > 0.
        self.__spinup_days = spinup_days

    def __set_forcing_start_date(self):
        forcing_start_date = self.start_date - timedelta(
            minutes=self.spinup_days*24.*60.)
        self.__forcing_start_date = forcing_start_date

    def __set_end_date(self, end_date):
        assert isinstance(end_date, datetime)
        if end_date <= self.start_date:
            raise Exception('end_date must be larger than start_date.')
        self.__end_date = end_date

    def __set_forcing_end_date(self):
        self.__forcing_end_date = self.end_date

    def __set_tpxo_constituents(self):
        tpxo_constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1',
                             'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1']
        self.__tpxo_constituents = dict()
        for index, constituent in enumerate(tpxo_constituents):
            self.tpxo_constituents[constituent] = index

    def __set_constituents(self, constituents):
        if constituents is None:
            self.__constituents = set()
        else:
            if isinstance(constituents, str):
                if constituents.lower() == 'major':
                    constituents = self.major8
                elif constituents.lower() == 'all':
                    constituents = self.tpxo_constituents
            constituents = list(constituents)
            self.__constituents = set()
            for constituent in constituents:
                if constituent not in self.tpxo_constituents:
                    raise Exception('Constituent requested ('
                                    + '{}'.format(constituent)
                                    + ') not found in TPXO database.')
                else:
                    self.__constituents.add(constituent)

    def __init_orbital_parameters(self):
        self.DYR = self.forcing_start_date.year - 1900.
        self.DDAY = (self.forcing_start_date.timetuple().tm_yday
                     + int((self.forcing_start_date.year-1901.)/4.) - 1)
        self.hour_middle = (self.forcing_start_date.hour
                            + (self.forcing_end_date
                               - self.forcing_start_date).total_seconds()
                            / 3600 / 2)
        self.DN = self.__get_lunar_node(self.hour_middle)  # HR
        self.N = np.deg2rad(self.DN)
        self.DP = self.__get_lunar_perigee(self.hour_middle)  # HR
        self.P = np.deg2rad(self.DP)
        self.DH = self.__get_solar_mean_longitude(
            self.forcing_start_date.hour)  # HR
        self.H = np.deg2rad(self.DH)
        self.DS = self.__get_lunar_mean_longitude(
            self.forcing_start_date.hour) % 360.  # HR
        self.S = np.deg2rad(self.DS)
        self.DP1 = self.__get_solar_perigee(self.forcing_start_date.hour)  # HR
        self.P1 = np.deg2rad(self.DP1)
        self.I = np.arccos(.9136949-.0356926*np.cos(self.N))
        self.DI = np.rad2deg(self.I)
        self.NU = np.arcsin(.0897056*np.sin(self.N)/np.sin(self.I))
        self.DNU = np.rad2deg(self.NU)
        self.XI = self.N-2.*np.arctan(.64412*np.tan(self.N/2)) - self.NU
        self.DXI = np.rad2deg(self.XI)

        self.DT = (180.+self.forcing_start_date.hour*(360./24))
        self.T = np.deg2rad(self.DT)
        self.NUP = np.arctan(np.sin(self.NU)
                             / (np.cos(self.NU)+.334766/np.sin(2.*self.I)))
        self.DNUP = np.rad2deg(self.NUP)
        self.DPC = (self.DP - self.DXI)
        self.PC = np.deg2rad(self.DPC)
        self.R = np.arctan(np.sin(2.*self.PC)
                           / ((1./6.)*(1./np.tan(.5*self.I))**2
                              - np.cos(2.*self.PC)))
        self.DR = np.rad2deg(self.R)
        self.NUP2 = np.arctan(np.sin(2.*self.NU)
                              / (np.cos(2.*self.NU)+.0726184
                                 / np.sin(self.I)**2))/2.
        self.DNUP2 = np.rad2deg(self.NUP2)
        self.Q = np.arctan2((5.*np.cos(self.I)-1.)*np.sin(self.PC),
                            (7.*np.cos(self.I)+1.)*np.cos(self.PC))
        self.DQ = np.rad2deg(self.Q)

    def __set_nodal_factor(self):
        self.__nodal_factor = CaseInsensitiveDict()
        for constituent in self.constituents:
            self.nodal_factor[constituent] \
                = self.__get_nodal_factor(constituent)

    def __set_greenwich_term(self):
        # greenwich terms are referenced to the spinup_date
        self.__greenwich_term = CaseInsensitiveDict()
        for constituent in self.constituents:
            self.greenwich_term[constituent] \
                = self.__get_greenwich_term(constituent) % 360.

    def __get_nodal_factor(self, constituent):
        if constituent == "M2":
            return self.__EQ78()
        elif constituent == "S2":
            return 1.0
        elif constituent == "N2":
            return self.__EQ78()
        elif constituent == "K1":
            return self.__EQ227()
        elif constituent == "M4":
            return (self.__EQ78())**2.
        elif constituent == "O1":
            return self.__EQ75()
        elif constituent == "M6":
            return (self.__EQ78())**3.
        elif constituent == "MK3":
            return self.__EQ78()*self.__EQ227()
        elif constituent == "S4":
            return 1.0
        elif constituent == "MN4":
            return (self.__EQ78())**2.
        elif constituent == "Nu2":
            return self.__EQ78()
        elif constituent == "S6":
            return 1.0
        elif constituent == "MU2":
            return self.__EQ78()
        elif constituent == "2N2":
            return self.__EQ78()
        elif constituent == "OO1":
            return self.__EQ77()
        elif constituent == "lambda2":
            return self.__EQ78()
        elif constituent == "S1":
            return 1.0
        elif constituent == "M1":
            return self.__EQ207()
        elif constituent == "J1":
            return self.__EQ76()
        elif constituent == "Mm":
            return self.__EQ73()
        elif constituent == "Ssa":
            return 1.0
        elif constituent == "Sa":
            return 1.0
        elif constituent == "Msf":
            return self.__EQ78()
        elif constituent == "Mf":
            return self.__EQ74()
        elif constituent == "RHO":
            return self.__EQ75()
        elif constituent == "Q1":
            return self.__EQ75()
        elif constituent == "T2":
            return 1.0
        elif constituent == "R2":
            return 1.0
        elif constituent == "2Q1":
            return self.__EQ75()
        elif constituent == "P1":
            return 1.0
        elif constituent == "2SM2":
            return self.__EQ78()
        elif constituent == "M3":
            return self.__EQ149()
        elif constituent == "L2":
            return self.__EQ215()
        elif constituent == "2MK3":
            return self.__EQ227()*self.__EQ78()**2
        elif constituent == "K2":
            return self.__EQ235()
        elif constituent == "M8":
            return self.__EQ78()**4
        elif constituent == "MS4":
            return self.__EQ78()
        else:
            raise RuntimeError('Unrecognized constituent '
                               + '{}.'.format(constituent))

    def __get_greenwich_term(self, constituent):
        if constituent == "M2":
            return 2.*(self.DT-self.DS+self.DH)+2.*(self.DXI-self.DNU)
        elif constituent == "S2":
            return 2.*self.DT
        elif constituent == "N2":
            return 2.*(self.DT+self.DH)-3.*self.DS+self.DP+2. \
                * (self.DXI-self.DNU)
        elif constituent == "K1":
            return self.DT+self.DH-90.-self.DNUP
        elif constituent == "M4":
            return 4.*(self.DT-self.DS+self.DH)+4.*(self.DXI-self.DNU)
        elif constituent == "O1":
            return self.DT-2.*self.DS+self.DH+90.+2.*self.DXI-self.DNU
        elif constituent == "M6":
            return 6.*(self.DT-self.DS+self.DH)+6.*(self.DXI-self.DNU)
        elif constituent == "MK3":
            return 3.*(self.DT+self.DH)-2.*self.DS-90.+2.*(self.DXI-self.DNU) \
                - self.DNUP
        elif constituent == "S4":
            return 4.*self.DT
        elif constituent == "MN4":
            return 4.*(self.DT+self.DH)-5.*self.DS+self.DP+4.\
                * (self.DXI-self.DNU)
        elif constituent == "Nu2":
            return 2.*self.DT-3.*self.DS+4.*self.DH-self.DP+2. \
                * (self.DXI-self.DNU)
        elif constituent == "S6":
            return 6.*self.DT
        elif constituent == "MU2":
            return 2.*(self.DT+2.*(self.DH-self.DS))+2.*(self.DXI-self.DNU)
        elif constituent == "2N2":
            return 2.*(self.DT-2.*self.DS+self.DH+self.DP)+2. \
                * (self.DXI-self.DNU)
        elif constituent == "OO1":
            return self.DT+2.*self.DS+self.DH-90.-2.*self.DXI-self.DNU
        elif constituent == "lambda2":
            return 2.*self.DT-self.DS+self.DP+180.+2.*(self.DXI-self.DNU)
        elif constituent == "S1":
            return self.DT
        elif constituent == "M1":
            return self.DT-self.DS+self.DH-90.+self.DXI-self.DNU+self.DQ
        elif constituent == "J1":
            return self.DT+self.DS+self.DH-self.DP-90.-self.DNU
        elif constituent == "Mm":
            return self.DS-self.DP
        elif constituent == "Ssa":
            return 2.*self.DH
        elif constituent == "Sa":
            return self.DH
        elif constituent == "Msf":
            return 2.*(self.DS-self.DH)
        elif constituent == "Mf":
            return 2.*self.DS-2.*self.DXI
        elif constituent == "RHO":
            return self.DT+3.*(self.DH-self.DS)-self.DP+90.+2.\
                * self.DXI-self.DNU
        elif constituent == "Q1":
            return self.DT-3.*self.DS+self.DH+self.DP+90.+2.*self.DXI-self.DNU
        elif constituent == "T2":
            return 2.*self.DT-self.DH+self.DP1
        elif constituent == "R2":
            return 2.*self.DT+self.DH-self.DP1+180.
        elif constituent == "2Q1":
            return self.DT-4.*self.DS+self.DH+2.*self.DP+90.+2.*self.DXI \
                - self.DNU
        elif constituent == "P1":
            return self.DT-self.DH+90.
        elif constituent == "2SM2":
            return 2.*(self.DT+self.DS-self.DH)+2.*(self.DNU-self.DXI)
        elif constituent == "M3":
            return 3.*(self.DT-self.DS+self.DH)+3.*(self.DXI-self.DNU)
        elif constituent == "L2":
            return 2.*(self.DT+self.DH)-self.DS-self.DP+180.+2.\
                * (self.DXI-self.DNU)-self.DR
        elif constituent == "2MK3":
            return 3.*(self.DT+self.DH)-4.*self.DS+90.+4.*(self.DXI-self.DNU) \
                + self.DNUP
        elif constituent == "K2":
            return 2.*(self.DT+self.DH)-2.*self.DNUP2
        elif constituent == "M8":
            return 8.*(self.DT-self.DS+self.DH)+8.*(self.DXI-self.DNU)
        elif constituent == "MS4":
            return 2.*(2.*self.DT-self.DS+self.DH)+2.*(self.DXI-self.DNU)
        else:
            raise RuntimeError('Unrecognized constituent '
                               + '{}.'.format(constituent))

    def __get_lunar_node(self, hours):
        return (259.1560564 - 19.328185764 * self.DYR - .0529539336 * self.DDAY
                - .0022064139 * hours)

    def __get_lunar_perigee(self, hours):
        return (334.3837214 + 40.66246584 * self.DYR + .111404016 * self.DDAY
                + .004641834 * hours)

    def __get_lunar_mean_longitude(self, hours):
        return (277.0256206 + 129.38482032 * self.DYR + 13.176396768
                * self.DDAY + .549016532 * hours)

    def __get_solar_perigee(self, hours):
        return (281.2208569 + .01717836 * self.DYR + .000047064 * self.DDAY
                + .000001961 * hours)

    def __get_solar_mean_longitude(self, hours):
        return (280.1895014 - .238724988 * self.DYR + .9856473288 * self.DDAY
                + .0410686387 * hours)

    def __EQ73(self):
        return (2./3.-np.sin(self.I)**2)/.5021

    def __EQ74(self):
        return np.sin(self.I)**2/.1578

    def __EQ75(self):
        return np.sin(self.I)*np.cos(self.I/2.)**2/.37988

    def __EQ76(self):
        return np.sin(2.*self.I)/.7214

    def __EQ77(self):
        return np.sin(self.I)*np.sin(self.I/2.)**2/.0164

    def __EQ78(self):
        return (np.cos(self.I/2)**4)/.91544

    def __EQ149(self):
        return np.cos(self.I/2.)**6/.8758

    def __EQ197(self):
        return np.sqrt(2.310+1.435*np.cos(2.*(self.P - self.XI)))

    def __EQ207(self):
        return self.__EQ75()*self.__EQ197()

    def __EQ213(self):
        return np.sqrt(1.-12.*np.tan(self.I/2.)**2*np.cos(2.*self.P)
                       + 36.*np.tan(self.I/2.)**4)

    def __EQ215(self):
        return self.__EQ78()*self.__EQ213()

    def __EQ227(self):
        return np.sqrt(.8965*np.sin(2.*self.I)**2+.6001*np.sin(2.*self.I)
                       * np.cos(self.NU)+.1006)

    def __EQ235(self):
        return .001+np.sqrt(19.0444*np.sin(self.I)**4+2.7702*np.sin(self.I)**2
                            * np.cos(2.*self.NU)+.0981)

    @property
    def units(self):
        return 'rad/sec'

    @property
    def nc(self):
        return self.__nc

    @property
    def start_date(self):
        return self.__start_date

    @property
    def end_date(self):
        return self.__end_date

    @property
    def spinup_days(self):
        return self.__spinup_days

    @property
    def forcing_start_date(self):
        return self.__forcing_start_date

    @property
    def forcing_end_date(self):
        return self.__forcing_end_date

    @property
    def constituents(self):
        return sorted(self.__constituents)

    @property
    def tidal_potential_amplitude(self):
        return CaseInsensitiveDict(self.__tidal_potential_amplitude)

    @property
    def orbital_frequency(self):
        return CaseInsensitiveDict(self.__orbital_frequency)

    @property
    def earth_tidal_potential(self):
        return CaseInsensitiveDict(self.__earth_tidal_potential)

    @property
    def nodal_factor(self):
        return self.__nodal_factor

    @property
    def TPXO_interp(self):
        return self.__TPXO_interp

    @property
    def tpxo_constituents(self):
        return self.__tpxo_constituents

    @property
    def major8(self):
        return self.__major8

    @property
    def greenwich_term(self):
        return self.__greenwich_term


class TidalForcing(unittest.TestCase):

    def test_control_Sandy(self):
        """
        K1    0.94843     298.0830551604721
        O1    0.91613     167.70941789973676
        P1    1.0         70.01137799920002
        Q1    0.91613     266.0233347217413
        N2    1.0203      208.1293472169341
        M2    1.0203      109.81543039492954
        S2    1.0         0.0
        K2    0.86465     55.559822189650276
        M4    1.041       219.63086078985907
        MS4   1.0203      109.81543039492954
        """
        start_date = datetime(2012, 10, 26, 00)
        end_date = start_date + timedelta(days=4.25)
        t = _TidalForcing(start_date, end_date, 15)
        ordered = ["K1", "O1", "P1", "Q1", "N2", "M2", "S2", "K2", "MF", "MM",
                   "M4", "MS4"]
        for constituent in ordered:
            if constituent in t.constituents:
                s = '{:<6}'.format(constituent)
                s += '{:<7.5}'.format(t.nodal_factor[constituent])
                s += 5*' '
                s += '{} '.format(t.greenwich_term[constituent])
                print(s)

    def _test_TPXO_interp(self):
        mesh = Model.AdcircMesh(os.getenv('ShinnecockInlet'), 4326, 'LMSL')
        now = datetime.now()
        three_days = now + timedelta(minutes=3.*24.*60.)
        TF = _TidalForcing(now, three_days, 7)
        TF.set_TPXO_interp(mesh)
