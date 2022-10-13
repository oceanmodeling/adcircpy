from datetime import datetime, timedelta
from enum import Enum
from os import PathLike
from typing import Union

import numpy as np

from adcircpy.forcing import bctypes
from adcircpy.forcing.tides.hamtide import HAMTIDE
from adcircpy.forcing.tides.tpxo import TPXO


class TidalSource(Enum):
    HAMTIDE = HAMTIDE
    TPXO = TPXO


class Tides(bctypes.EtaBc):
    def __init__(
        self, tidal_source: Union[str, TidalSource] = None, resource: PathLike = None
    ):
        if tidal_source is None:
            tidal_source = TidalSource.HAMTIDE
        elif isinstance(tidal_source, str):
            try:
                tidal_source = TidalSource[tidal_source.upper()]
            except:
                raise NotImplementedError(
                    f'tidal source {tidal_source} not recognized; '
                    f'must be one of {[entry.__name__ for entry in TidalSource]}'
                )

        self._active_constituents = {}

        self.tidal_source = tidal_source
        self.tidal_dataset = tidal_source.value(resource)

    def __call__(self, constituent: str) -> (float, float, float, float, float):
        return self.get_tidal_constituent(constituent)

    def __iter__(self):
        for constituent in self.active_constituents:
            yield constituent, self.get_tidal_constituent(constituent)

    def __len__(self) -> int:
        return len(self.active_constituents)

    def use_all(self, potential=True, forcing=True):
        for constituent in self.constituents:
            if constituent not in self.major_constituents:
                potential = False
            self._active_constituents[constituent] = {
                'potential': potential,
                'forcing': forcing,
            }

    def use_major(self, potential=True, forcing=True):
        for constituent in self.major_constituents:
            self.use_constituent(constituent, potential, forcing)

    def use_constituent(self, constituent, potential=True, forcing=True):
        msg = 'Constituent must be one of '
        msg += f'{self.constituents}'
        assert constituent in self.constituents, msg
        if constituent not in self.major_constituents:
            potential = False
        self._active_constituents[constituent] = {
            'potential': potential,
            'forcing': forcing,
        }

    def drop_constituent(self, constituent):
        msg = 'constituent must be one of: '
        msg += f'{self.active_constituents}'
        assert constituent in self.active_constituents, msg
        self._active_constituents.pop(constituent)

    def get_active_constituents(self) -> [str]:
        return list(self.active_constituents.keys())

    def get_active_forcing_constituents(self) -> [str]:
        fc = []
        for c, d in self.active_constituents.items():
            if d['forcing']:
                fc.append(c)
        return fc

    def get_active_potential_constituents(self) -> [str]:
        fc = []
        for c, d in self.active_constituents.items():
            if d['potential']:
                fc.append(c)
        return fc

    def get_tidal_potential_amplitude(self, constituent) -> float:
        if constituent in self.tidal_potential_amplitudes:
            return self.tidal_potential_amplitudes[constituent]

    def get_tidal_species_type(self, constituent) -> int:
        if constituent in self.tidal_species_type:
            return self.tidal_species_type[constituent]
        return 0

    def get_orbital_frequency(self, constituent) -> float:
        return self.orbital_frequencies[constituent]

    def get_tidal_constituent(self, constituent):
        return (
            self.get_tidal_potential_amplitude(constituent),
            self.get_orbital_frequency(constituent),
            self.get_earth_tidal_potential(constituent),
            self.get_nodal_factor(constituent),
            self.get_greenwich_factor(constituent),
        )

    def get_earth_tidal_potential(self, constituent) -> float:
        try:
            return self.earth_tidal_potentials[constituent]
        except KeyError:
            pass

    def get_nodal_factor(self, constituent) -> float:
        if constituent == 'M2':
            return self.EQ78
        elif constituent == 'S2':
            return 1.0
        elif constituent == 'N2':
            return self.EQ78
        elif constituent == 'K1':
            return self.EQ227
        elif constituent == 'M4':
            return (self.EQ78) ** 2.0
        elif constituent == 'O1':
            return self.EQ75
        elif constituent == 'M6':
            return (self.EQ78) ** 3.0
        elif constituent == 'MK3':
            return self.EQ78 * self.EQ227
        elif constituent == 'S4':
            return 1.0
        elif constituent == 'MN4':
            return (self.EQ78) ** 2.0
        elif constituent == 'Nu2':
            return self.EQ78
        elif constituent == 'S6':
            return 1.0
        elif constituent == 'MU2':
            return self.EQ78
        elif constituent == '2N2':
            return self.EQ78
        elif constituent == 'OO1':
            return self.EQ77
        elif constituent == 'lambda2':
            return self.EQ78
        elif constituent == 'S1':
            return 1.0
        elif constituent == 'M1':
            return self.EQ207
        elif constituent == 'J1':
            return self.EQ76
        elif constituent == 'Mm':
            return self.EQ73
        elif constituent == 'Ssa':
            return 1.0
        elif constituent == 'Sa':
            return 1.0
        elif constituent == 'Msf':
            return self.EQ78
        elif constituent == 'Mf':
            return self.EQ74
        elif constituent == 'RHO':
            return self.EQ75
        elif constituent == 'Q1':
            return self.EQ75
        elif constituent == 'T2':
            return 1.0
        elif constituent == 'R2':
            return 1.0
        elif constituent == '2Q1':
            return self.EQ75
        elif constituent == 'P1':
            return 1.0
        elif constituent == '2SM2':
            return self.EQ78
        elif constituent == 'M3':
            return self.EQ149
        elif constituent == 'L2':
            return self.EQ215
        elif constituent == '2MK3':
            return self.EQ227 * self.EQ78 ** 2
        elif constituent == 'K2':
            return self.EQ235
        elif constituent == 'M8':
            return self.EQ78 ** 4
        elif constituent == 'MS4':
            return self.EQ78
        else:
            msg = f'Unrecognized constituent {constituent}'
            raise TypeError(msg)

    def _normalize_to_360(f):
        def decorator(self, constituent):
            return f(self, constituent) % 360.0

        return decorator

    @_normalize_to_360
    def get_greenwich_factor(self, constituent) -> float:
        if constituent == 'M2':
            return 2.0 * (self.DT - self.DS + self.DH) + 2.0 * (self.DXI - self.DNU)
        elif constituent == 'S2':
            return 2.0 * self.DT
        elif constituent == 'N2':
            return (
                2.0 * (self.DT + self.DH)
                - 3.0 * self.DS
                + self.DP
                + 2.0 * (self.DXI - self.DNU)
            )
        elif constituent == 'K1':
            return self.DT + self.DH - 90.0 - self.DNUP
        elif constituent == 'M4':
            return 4.0 * (self.DT - self.DS + self.DH) + 4.0 * (self.DXI - self.DNU)
        elif constituent == 'O1':
            return self.DT - 2.0 * self.DS + self.DH + 90.0 + 2.0 * self.DXI - self.DNU
        elif constituent == 'M6':
            return 6.0 * (self.DT - self.DS + self.DH) + 6.0 * (self.DXI - self.DNU)
        elif constituent == 'MK3':
            return (
                3.0 * (self.DT + self.DH)
                - 2.0 * self.DS
                - 90.0
                + 2.0 * (self.DXI - self.DNU)
                - self.DNUP
            )
        elif constituent == 'S4':
            return 4.0 * self.DT
        elif constituent == 'MN4':
            return (
                4.0 * (self.DT + self.DH)
                - 5.0 * self.DS
                + self.DP
                + 4.0 * (self.DXI - self.DNU)
            )
        elif constituent == 'Nu2':
            return (
                2.0 * self.DT
                - 3.0 * self.DS
                + 4.0 * self.DH
                - self.DP
                + 2.0 * (self.DXI - self.DNU)
            )
        elif constituent == 'S6':
            return 6.0 * self.DT
        elif constituent == 'MU2':
            return 2.0 * (self.DT + 2.0 * (self.DH - self.DS)) + 2.0 * (self.DXI - self.DNU)
        elif constituent == '2N2':
            return 2.0 * (self.DT - 2.0 * self.DS + self.DH + self.DP) + 2.0 * (
                self.DXI - self.DNU
            )
        elif constituent == 'OO1':
            return self.DT + 2.0 * self.DS + self.DH - 90.0 - 2.0 * self.DXI - self.DNU
        elif constituent == 'lambda2':
            return 2.0 * self.DT - self.DS + self.DP + 180.0 + 2.0 * (self.DXI - self.DNU)
        elif constituent == 'S1':
            return self.DT
        elif constituent == 'M1':
            return self.DT - self.DS + self.DH - 90.0 + self.DXI - self.DNU + self.DQ
        elif constituent == 'J1':
            return self.DT + self.DS + self.DH - self.DP - 90.0 - self.DNU
        elif constituent == 'Mm':
            return self.DS - self.DP
        elif constituent == 'Ssa':
            return 2.0 * self.DH
        elif constituent == 'Sa':
            return self.DH
        elif constituent == 'Msf':
            return 2.0 * (self.DS - self.DH)
        elif constituent == 'Mf':
            return 2.0 * self.DS - 2.0 * self.DXI
        elif constituent == 'RHO':
            return (
                self.DT
                + 3.0 * (self.DH - self.DS)
                - self.DP
                + 90.0
                + 2.0 * self.DXI
                - self.DNU
            )
        elif constituent == 'Q1':
            return (
                self.DT - 3.0 * self.DS + self.DH + self.DP + 90.0 + 2.0 * self.DXI - self.DNU
            )
        elif constituent == 'T2':
            return 2.0 * self.DT - self.DH + self.DP1
        elif constituent == 'R2':
            return 2.0 * self.DT + self.DH - self.DP1 + 180.0
        elif constituent == '2Q1':
            return (
                self.DT
                - 4.0 * self.DS
                + self.DH
                + 2.0 * self.DP
                + 90.0
                + 2.0 * self.DXI
                - self.DNU
            )
        elif constituent == 'P1':
            return self.DT - self.DH + 90.0
        elif constituent == '2SM2':
            return 2.0 * (self.DT + self.DS - self.DH) + 2.0 * (self.DNU - self.DXI)
        elif constituent == 'M3':
            return 3.0 * (self.DT - self.DS + self.DH) + 3.0 * (self.DXI - self.DNU)
        elif constituent == 'L2':
            return (
                2.0 * (self.DT + self.DH)
                - self.DS
                - self.DP
                + 180.0
                + 2.0 * (self.DXI - self.DNU)
                - self.DR
            )
        elif constituent == '2MK3':
            return (
                3.0 * (self.DT + self.DH)
                - 4.0 * self.DS
                + 90.0
                + 4.0 * (self.DXI - self.DNU)
                + self.DNUP
            )
        elif constituent == 'K2':
            return 2.0 * (self.DT + self.DH) - 2.0 * self.DNUP2
        elif constituent == 'M8':
            return 8.0 * (self.DT - self.DS + self.DH) + 8.0 * (self.DXI - self.DNU)
        elif constituent == 'MS4':
            return 2.0 * (2.0 * self.DT - self.DS + self.DH) + 2.0 * (self.DXI - self.DNU)
        else:
            msg = f'Unrecognized constituent {constituent}'
            raise TypeError(msg)

    def get_lunar_node(self) -> float:
        return (
            259.1560564
            - 19.328185764 * self.DYR
            - 0.0529539336 * self.DDAY
            - 0.0022064139 * self.hour_middle
        )

    def get_lunar_perigee(self) -> float:
        return (
            334.3837214
            + 40.66246584 * self.DYR
            + 0.111404016 * self.DDAY
            + 0.004641834 * self.hour_middle
        )

    def get_lunar_mean_longitude(self) -> float:
        return (
            277.0256206
            + 129.38482032 * self.DYR
            + 13.176396768 * self.DDAY
            + 0.549016532 * self.forcing_start_date.hour
        )

    def get_solar_perigee(self) -> float:
        return (
            281.2208569
            + 0.01717836 * self.DYR
            + 0.000047064 * self.DDAY
            + 0.000001961 * self.start_date.hour
        )

    def get_solar_mean_longitude(self) -> float:
        return (
            280.1895014
            - 0.238724988 * self.DYR
            + 0.9856473288 * self.DDAY
            + 0.0410686387 * self.start_date.hour
        )

    @property
    def EQ73(self) -> float:
        """ """
        return (2.0 / 3.0 - np.sin(self.I) ** 2) / 0.5021

    @property
    def EQ74(self) -> float:
        """ """
        return np.sin(self.I) ** 2 / 0.1578

    @property
    def EQ75(self) -> float:
        """ """
        return np.sin(self.I) * np.cos(self.I / 2.0) ** 2 / 0.37988

    @property
    def EQ76(self) -> float:
        """ """
        return np.sin(2.0 * self.I) / 0.7214

    @property
    def EQ77(self) -> float:
        """ """
        return np.sin(self.I) * np.sin(self.I / 2.0) ** 2 / 0.0164

    @property
    def EQ78(self) -> float:
        """ """
        return (np.cos(self.I / 2) ** 4) / 0.91544

    @property
    def EQ149(self) -> float:
        """ """
        return np.cos(self.I / 2.0) ** 6 / 0.8758

    @property
    def EQ197(self) -> float:
        """ """
        return np.sqrt(2.310 + 1.435 * np.cos(2.0 * (self.P - self.XI)))

    @property
    def EQ207(self) -> float:
        """ """
        return self.EQ75 * self.EQ197

    @property
    def EQ213(self) -> float:
        """ """
        return np.sqrt(
            1.0
            - 12.0 * np.tan(self.I / 2.0) ** 2 * np.cos(2.0 * self.P)
            + 36.0 * np.tan(self.I / 2.0) ** 4
        )

    @property
    def EQ215(self) -> float:
        """ """
        return self.EQ78 * self.EQ213

    @property
    def EQ227(self) -> float:
        return np.sqrt(
            0.8965 * np.sin(2.0 * self.I) ** 2
            + 0.6001 * np.sin(2.0 * self.I) * np.cos(self.NU)
            + 0.1006
        )

    @property
    def EQ235(self) -> float:
        return 0.001 + np.sqrt(
            19.0444 * np.sin(self.I) ** 4
            + 2.7702 * np.sin(self.I) ** 2 * np.cos(2.0 * self.NU)
            + 0.0981
        )

    @property
    def start_date(self) -> datetime:
        try:
            return self.__start_date
        except AttributeError:
            msg = 'Must set start_date attribute.'
            raise AttributeError(msg)

    @start_date.setter
    def start_date(self, start_date: datetime):
        if start_date is None:
            del self.start_date
            return
        msg = f'start_date must be an instance of type {datetime}.'
        assert isinstance(start_date, datetime), msg
        try:
            msg = 'start_date must be smaller than end_date.'
            assert start_date < self.end_date, msg
        except AttributeError:
            pass
        self.__start_date = start_date

    @start_date.deleter
    def start_date(self):
        try:
            del self.__start_date
        except AttributeError:
            pass

    @property
    def end_date(self) -> datetime:
        try:
            return self.__end_date
        except AttributeError:
            msg = 'Must set end_date attribute.'
            raise AttributeError(msg)

    @end_date.setter
    def end_date(self, end_date: datetime):
        if end_date is None:
            del self.end_date
            return
        msg = f'end_date must be an instance of type {datetime}.'
        assert isinstance(end_date, datetime), msg
        try:
            msg = f'end_date ({end_date}) must be larger than '
            msg += f'start_date ({self.start_date}).'
            assert end_date > self.start_date, msg
        except AttributeError:
            pass
        self.__end_date = end_date

    @end_date.deleter
    def end_date(self):
        try:
            del self.__end_date
        except AttributeError:
            pass

    @property
    def forcing_start_date(self) -> datetime:
        return self.start_date - self.spinup_time

    @property
    def spinup_time(self) -> timedelta:
        try:
            return self.__spinup_time
        except AttributeError:
            return timedelta(0.0)

    @spinup_time.setter
    def spinup_time(self, spinup_time: timedelta):
        if spinup_time is None:
            del self.spinup_time
            return
        msg = f'spinup_time must be of and instance of type {timedelta}.'
        assert isinstance(spinup_time, timedelta), msg
        self.__spinup_time = np.abs(spinup_time)

    @spinup_time.deleter
    def spinup_time(self):
        try:
            del self.__spinup_time
        except AttributeError:
            pass

    @property
    def active_constituents(self) -> [str]:
        return self._active_constituents.copy()

    @property
    def major_constituents(self) -> [str]:
        return ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']

    @property
    def constituents(self) -> [str]:
        constituents = self.major_constituents
        if self.tidal_source == TidalSource.TPXO:
            constituents.extend(['Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1'])
        return constituents

    @property
    def orbital_frequencies(self) -> {str: float}:
        return {
            'M4': 0.0002810378050173,
            'M6': 0.0004215567080107,
            'MK3': 0.0002134400613513,
            'S4': 0.0002908882086657,
            'MN4': 0.0002783986019952,
            'S6': 0.0004363323129986,
            'M3': 0.0002107783537630,
            '2MK3': 0.0002081166466594,
            'M8': 0.0005620756090649,
            'MS4': 0.0002859630068415,
            'M2': 0.0001405189025086,
            'S2': 0.0001454441043329,
            'N2': 0.0001378796994865,
            'Nu2': 0.0001382329037065,
            'MU2': 0.0001355937006844,
            '2N2': 0.0001352404964644,
            'lambda2': 0.0001428049013108,
            'T2': 0.0001452450073529,
            'R2': 0.0001456432013128,
            '2SM2': 0.0001503693061571,
            'L2': 0.0001431581055307,
            'K2': 0.0001458423172006,
            'K1': 0.0000729211583579,
            'O1': 0.0000675977441508,
            'OO1': 0.0000782445730498,
            'S1': 0.0000727220521664,
            'M1': 0.0000702594512543,
            'J1': 0.0000755603613800,
            'RHO': 0.0000653117453487,
            'Q1': 0.0000649585411287,
            '2Q1': 0.0000623193381066,
            'P1': 0.0000725229459750,
            'Mm': 0.0000026392030221,
            'Ssa': 0.0000003982128677,
            'Sa': 0.0000001991061914,
            'Msf': 0.0000049252018242,
            'Mf': 0.0000053234146919,
        }

    @property
    def tidal_potential_amplitudes(self) -> {str: float}:
        return {
            'M2': 0.242334,
            'S2': 0.112841,
            'N2': 0.046398,
            'K2': 0.030704,
            'K1': 0.141565,
            'O1': 0.100514,
            'P1': 0.046843,
            'Q1': 0.019256,
        }

    @property
    def tidal_species_type(self) -> {str: float}:
        return {'M2': 2, 'S2': 2, 'N2': 2, 'K2': 2, 'K1': 1, 'O1': 1, 'P1': 1, 'Q1': 1}

    @property
    def earth_tidal_potentials(self) -> {str: float}:
        return {
            'M2': 0.693,
            'S2': 0.693,
            'N2': 0.693,
            'K2': 0.693,
            'K1': 0.736,
            'O1': 0.695,
            'P1': 0.706,
            'Q1': 0.695,
        }

    @property
    def hour_middle(self) -> float:
        return self.forcing_start_date.hour + (
            (self.end_date - self.forcing_start_date) / timedelta(hours=1) / 2
        )

    @property
    def I(self) -> float:  # noqa:E743
        return np.arccos(0.9136949 - 0.0356926 * np.cos(self.N))

    @property
    def N(self) -> float:
        return np.deg2rad(self.DN)

    @property
    def DN(self) -> float:
        return self.get_lunar_node()

    @property
    def DYR(self) -> float:
        return self.forcing_start_date.year - 1900.0

    @property
    def DDAY(self) -> float:
        return (
            self.forcing_start_date.timetuple().tm_yday
            + int((self.forcing_start_date.year - 1901.0) / 4.0)
            - 1
        )

    @property
    def NU(self) -> float:
        return np.arcsin(0.0897056 * np.sin(self.N) / np.sin(self.I))

    @property
    def DT(self) -> float:
        return 180.0 + self.start_date.hour * (360.0 / 24)

    @property
    def DS(self) -> float:
        return self.get_lunar_mean_longitude()

    @property
    def DP(self) -> float:
        return self.get_lunar_perigee()

    @property
    def P(self) -> float:
        return np.deg2rad(self.DP)

    @property
    def DH(self) -> float:
        return self.get_solar_mean_longitude()

    @property
    def DP1(self) -> float:
        return self.get_solar_perigee()  # HR

    @property
    def DNU(self) -> float:
        return np.rad2deg(self.NU)

    @property
    def XI(self) -> float:
        return self.N - 2.0 * np.arctan(0.64412 * np.tan(self.N / 2)) - self.NU

    @property
    def DXI(self) -> float:
        return np.rad2deg(self.XI)

    @property
    def NUP(self) -> float:
        return np.arctan(np.sin(self.NU) / (np.cos(self.NU) + 0.334766 / np.sin(2.0 * self.I)))

    @property
    def DNUP(self) -> float:
        return np.rad2deg(self.NUP)

    @property
    def DPC(self) -> float:
        return self.DP - self.DXI

    @property
    def PC(self) -> float:
        return np.deg2rad(self.DPC)

    @property
    def R(self) -> float:
        return np.arctan(
            np.sin(2.0 * self.PC)
            / ((1.0 / 6.0) * (1.0 / np.tan(0.5 * self.I)) ** 2 - np.cos(2.0 * self.PC))
        )

    @property
    def DR(self) -> float:
        return np.rad2deg(self.R)

    @property
    def NUP2(self) -> float:
        return (
            np.arctan(
                np.sin(2.0 * self.NU)
                / (np.cos(2.0 * self.NU) + 0.0726184 / np.sin(self.I) ** 2)
            )
            / 2.0
        )

    @property
    def DNUP2(self) -> float:
        return np.rad2deg(self.NUP2)

    @property
    def Q(self) -> float:
        return np.arctan2(
            (5.0 * np.cos(self.I) - 1.0) * np.sin(self.PC),
            (7.0 * np.cos(self.I) + 1.0) * np.cos(self.PC),
        )

    @property
    def DQ(self) -> float:
        return np.rad2deg(self.Q)

    @property
    def ntip(self) -> int:
        return len(self.get_active_potential_constituents())

    @property
    def cutoff_depth(self) -> float:
        if self.ntip == 0:
            return 0
        try:
            return self.__cutoff_depth
        except AttributeError:
            return 40.0

    @cutoff_depth.setter
    def cutoff_depth(self, cutoff_depth: float):
        assert isinstance(cutoff_depth, (int, float))
        self.__cutoff_depth = cutoff_depth

    @property
    def nbfr(self) -> int:
        return len(self.get_active_forcing_constituents())

    @property
    def btype(self) -> str:
        return 'iettype'

    @property
    def iettype(self) -> int:
        return 3

    def __eq__(self, other: 'Tides') -> bool:
        return self.tidal_dataset == other.tidal_dataset
