from datetime import datetime, timedelta
import numpy as np
from ordered_set import OrderedSet
# from scipy.interpolate import griddata
# from netCDF4 import Dataset
# from requests.structures import CaseInsensitiveDict


class TidalForcing:
    __orbital_frequencies = {'M4':      0.0002810378050173,
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

    __tidal_potential_amplitudes = {'M2': 0.242334,
                                    'S2': 0.112841,
                                    'N2': 0.046398,
                                    'K2': 0.030704,
                                    'K1': 0.141565,
                                    'O1': 0.100514,
                                    'P1': 0.046843,
                                    'Q1': 0.019256}
    __earth_tidal_potentials = {'M2': 0.693,
                                'S2': 0.693,
                                'N2': 0.693,
                                'K2': 0.693,
                                'K1': 0.736,
                                'O1': 0.695,
                                'P1': 0.706,
                                'Q1': 0.695}
    __major_constituents = ['Q1',
                            'O1',
                            'P1',
                            'K1',
                            'N2',
                            'M2',
                            'S2',
                            'K2']

    __constituents = [*__major_constituents,
                      'Mm',
                      'Mf',
                      'M4',
                      'MN4',
                      'MS4',
                      '2N2',
                      'S1']

    def __init__(self):
        self.__active_constituents = OrderedSet()
        self.__start_date = None
        self.__end_date = None
        self.__spinup_time = timedelta(0.)

    def __call__(self, constituent):
        return self.get_tidal_constituent(constituent)

    def __iter__(self):
        for constituent in self.__active_constituents:
            yield (constituent, self.get_tidal_constituent(constituent))

    def __len__(self):
        return len(self.__active_constituents)

    def use_all(self):
        for constituent in self.__constituents:
            self.__active_constituents.append(constituent)

    def use_major(self):
        for constituent in self.__major_constituents:
            self.use_constituent(constituent)

    def use_constituent(self, constituent):
        msg = "Constituent must be one of "
        msg += "{}".format(self.__constituents)
        assert constituent in self.__constituents, msg
        self.__active_constituents.add(constituent)

    def drop_constituent(self, constituent):
        msg = "constituent must be one of: "
        msg += "{}".format(self.active_constituents)
        assert constituent in self.active_constituents, msg
        self.__active_constituents.pop(
            self.active_constituents.index(constituent))

    def get_active_constituents(self):
        return self.active_constituents

    def get_tidal_potential_amplitude(self, constituent):
        try:
            return self.tidal_potential_amplitudes[constituent]
        except KeyError:
            pass

    def get_earth_tidal_potential(self, constituent):
        try:
            return self.earth_tidal_potentials[constituent]
        except KeyError:
            pass

    def get_orbital_frequency(self, constituent):
        return self.orbital_frequencies[constituent]

    def get_tidal_constituent(self, constituent):
        return (self.get_tidal_potential_amplitude(constituent),
                self.get_orbital_frequency(constituent),
                self.get_earth_tidal_potential(constituent),
                self.get_nodal_factor(constituent),
                self.get_greenwich_term(constituent))

    def get_nodal_factor(self, constituent):
        if constituent == "M2":
            return self.EQ78
        elif constituent == "S2":
            return 1.0
        elif constituent == "N2":
            return self.EQ78
        elif constituent == "K1":
            return self.EQ227
        elif constituent == "M4":
            return (self.EQ78)**2.
        elif constituent == "O1":
            return self.EQ75
        elif constituent == "M6":
            return (self.EQ78)**3.
        elif constituent == "MK3":
            return self.EQ78*self.EQ227
        elif constituent == "S4":
            return 1.0
        elif constituent == "MN4":
            return (self.EQ78)**2.
        elif constituent == "Nu2":
            return self.EQ78
        elif constituent == "S6":
            return 1.0
        elif constituent == "MU2":
            return self.EQ78
        elif constituent == "2N2":
            return self.EQ78
        elif constituent == "OO1":
            return self.EQ77
        elif constituent == "lambda2":
            return self.EQ78
        elif constituent == "S1":
            return 1.0
        elif constituent == "M1":
            return self.EQ207
        elif constituent == "J1":
            return self.EQ76
        elif constituent == "Mm":
            return self.EQ73
        elif constituent == "Ssa":
            return 1.0
        elif constituent == "Sa":
            return 1.0
        elif constituent == "Msf":
            return self.EQ78
        elif constituent == "Mf":
            return self.EQ74
        elif constituent == "RHO":
            return self.EQ75
        elif constituent == "Q1":
            return self.EQ75
        elif constituent == "T2":
            return 1.0
        elif constituent == "R2":
            return 1.0
        elif constituent == "2Q1":
            return self.EQ75
        elif constituent == "P1":
            return 1.0
        elif constituent == "2SM2":
            return self.EQ78
        elif constituent == "M3":
            return self.EQ149
        elif constituent == "L2":
            return self.EQ215
        elif constituent == "2MK3":
            return self.EQ227*self.EQ78**2
        elif constituent == "K2":
            return self.EQ235
        elif constituent == "M8":
            return self.EQ78**4
        elif constituent == "MS4":
            return self.EQ78
        else:
            raise RuntimeError('Unrecognized constituent '
                               + '{}.'.format(constituent))

    def normalize_to_360(f):
        def decorator(self, constituent):
            return f(self, constituent) % 360.
        return decorator

    @normalize_to_360
    def get_greenwich_term(self, constituent):
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

    def get_lunar_node(self):
        return (259.1560564 - 19.328185764 * self.DYR - .0529539336 * self.DDAY
                - .0022064139 * self.hour_middle)

    def get_lunar_perigee(self):
        return (334.3837214 + 40.66246584 * self.DYR + .111404016 * self.DDAY
                + .004641834 * self.hour_middle)

    def get_lunar_mean_longitude(self):
        return (277.0256206 + 129.38482032 * self.DYR + 13.176396768
                * self.DDAY + .549016532 * self.forcing_start_date.hour)

    def get_solar_perigee(self):
        return (281.2208569 + .01717836 * self.DYR + .000047064 * self.DDAY
                + .000001961 * self.start_date.hour)

    def get_solar_mean_longitude(self):
        return (280.1895014 - .238724988 * self.DYR + .9856473288 * self.DDAY
                + .0410686387 * self.start_date.hour)

    @property
    def EQ73(self):
        """ """
        return (2./3.-np.sin(self.I)**2)/.5021

    @property
    def EQ74(self):
        """ """
        return np.sin(self.I)**2/.1578

    @property
    def EQ75(self):
        """ """
        return np.sin(self.I)*np.cos(self.I/2.)**2/.37988

    @property
    def EQ76(self):
        """ """
        return np.sin(2.*self.I)/.7214

    @property
    def EQ77(self):
        """ """
        return np.sin(self.I)*np.sin(self.I/2.)**2/.0164

    @property
    def EQ78(self):
        """ """
        return (np.cos(self.I/2)**4)/.91544

    @property
    def EQ149(self):
        """ """
        return np.cos(self.I/2.)**6/.8758

    @property
    def EQ197(self):
        """ """
        return np.sqrt(2.310+1.435*np.cos(2.*(self.P - self.XI)))

    @property
    def EQ207(self):
        """ """
        return self.EQ75*self.EQ197

    @property
    def EQ213(self):
        """ """
        return np.sqrt(1.-12.*np.tan(self.I/2.)**2*np.cos(2.*self.P)
                       + 36.*np.tan(self.I/2.)**4)

    @property
    def EQ215(self):
        """ """
        return self.EQ78*self.EQ213

    @property
    def EQ227(self):
        """ """
        return np.sqrt(.8965*np.sin(2.*self.I)**2+.6001*np.sin(2.*self.I)
                       * np.cos(self.NU)+.1006)

    @property
    def EQ235(self):
        """ """
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
        try:
            return self.__start_date
        except AttributeError:
            raise AttributeError(
                "Must set start_date before operation can take place.")

    @property
    def end_date(self):
        try:
            return self.__end_date
        except AttributeError:
            raise AttributeError(
                "Must set end_date before operation can take place.")

    @property
    def forcing_start_date(self):
        return self.start_date - self.spinup_time

    @property
    def spinup_time(self):
        try:
            return self.__spinup_time
        except AttributeError:
            return timedelta(0.)

    @property
    def active_constituents(self):
        return self.__active_constituents

    @property
    def major_constituents(self):
        return self.__major_constituents

    @property
    def all_constituents(self):
        return self.__constituents

    @property
    def orbital_frequencies(self):
        return self.__orbital_frequencies

    @property
    def tidal_potential_amplitudes(self):
        return self.__tidal_potential_amplitudes

    @property
    def earth_tidal_potentials(self):
        return self.__earth_tidal_potentials

    @property
    def hour_middle(self):
        return self.forcing_start_date.hour + (
            (self.end_date - self.forcing_start_date).total_seconds()
            / 3600 / 2)

    @property
    def I(self):  # noqa:E743
        return np.arccos(.9136949-.0356926*np.cos(self.N))

    @property
    def N(self):
        return np.deg2rad(self.DN)

    @property
    def DN(self):
        return self.get_lunar_node()

    @property
    def DYR(self):
        return self.forcing_start_date.year - 1900.

    @property
    def DDAY(self):
        return (self.forcing_start_date.timetuple().tm_yday
                + int((self.forcing_start_date.year-1901.)/4.) - 1)

    @property
    def NU(self):
        return np.arcsin(.0897056*np.sin(self.N)/np.sin(self.I))

    @property
    def DT(self):
        return (180.+self.start_date.hour*(360./24))

    @property
    def DS(self):
        return self.get_lunar_mean_longitude()

    @property
    def DP(self):
        return self.get_lunar_perigee()

    @property
    def P(self):
        return np.deg2rad(self.DP)

    @property
    def DH(self):
        return self.get_solar_mean_longitude()

    @property
    def H(self):
        return np.deg2rad(self.DH)

    @property
    def S(self):
        return np.deg2rad(self.DS)

    @property
    def DP1(self):
        return self.get_solar_perigee()  # HR

    @property
    def P1(self):
        return np.deg2rad(self.DP1)

    @property
    def DI(self):
        return np.rad2deg(self.I)

    @property
    def DNU(self):
        return np.rad2deg(self.NU)

    @property
    def XI(self):
        return self.N-2.*np.arctan(.64412*np.tan(self.N/2)) - self.NU

    @property
    def DXI(self):
        return np.rad2deg(self.XI)

    @property
    def T(self):
        return np.deg2rad(self.DT)

    @property
    def NUP(self):
        return np.arctan(np.sin(self.NU) / (
            np.cos(self.NU)+.334766/np.sin(2.*self.I)))

    @property
    def DNUP(self):
        return np.rad2deg(self.NUP)

    @property
    def DPC(self):
        return self.DP - self.DXI

    @property
    def PC(self):
        return np.deg2rad(self.DPC)

    @property
    def R(self):
        return (np.arctan(np.sin(2.*self.PC)
                / ((1./6.)*(1./np.tan(.5*self.I))**2 - np.cos(2.*self.PC))))

    @property
    def DR(self):
        return np.rad2deg(self.R)

    @property
    def NUP2(self):
        return (np.arctan(np.sin(2.*self.NU)
                / (np.cos(2.*self.NU)+.0726184 / np.sin(self.I)**2))/2.)

    @property
    def DNUP2(self):
        return np.rad2deg(self.NUP2)

    @property
    def Q(self):
        return np.arctan2((5.*np.cos(self.I)-1.)*np.sin(self.PC),
                          (7.*np.cos(self.I)+1.)*np.cos(self.PC))

    @property
    def DQ(self):
        return np.rad2deg(self.Q)

    @start_date.setter
    def start_date(self, start_date):
        assert isinstance(start_date, (datetime, type(None))), \
            "start_date must be a datetime instance."
        if self.__end_date is not None and start_date is not None:
            assert start_date < self.end_date
        self.__start_date = start_date

    @end_date.setter
    def end_date(self, end_date):
        assert isinstance(end_date, (datetime, type(None))), \
            "end_date must be a datetime instance."
        if self.__start_date is not None and end_date is not None:
            assert end_date > self.start_date
        self.__end_date = end_date

    @spinup_time.setter
    def spinup_time(self, spinup_time):
        msg = "spinup_time must be timedelta object"
        assert isinstance(spinup_time, timedelta), msg
        self.__spinup_time = np.abs(spinup_time)