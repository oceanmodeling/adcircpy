from datetime import datetime, timedelta
from functools import lru_cache
import math
from os import PathLike
import pathlib

import numpy as np

from adcircpy.forcing.tides.tpxo import TPXO
from adcircpy.mesh.mesh import AdcircMesh


class EnumeratedValue:
    def __init__(self, values: []):
        self.values = values
        self.value = None

    def __set__(self, instance, value):
        if value not in self.values:
            raise ValueError(f'value must be one of {self.values}')
        self.value = value


class DefaultedValue:
    def __init__(self, default=None):
        self.default = default
        self.value = None

    def __get__(self, instance, owner):
        if self.value is not None:
            return self.value
        elif isinstance(self.default, Callable):
            return self.default()
        else:
            return self.default

    def __set__(self, instance, value):
        self.value = value


class DefaultedTypedValue:
    def __init__(self, value_type: type, default=None):
        self.type = value_type
        if default is not None and not isinstance(default, self.type):
            default = self.type(default)
        super().__init__(default)

    def __set__(self, instance, value):
        if value is not None and not isinstance(value, self.type):
            value = self.type(value)
        super().__set__(instance, value)


class DefaultedString(DefaultedTypedValue):
    def __init__(self, default: str = None):
        super().__init__(str, default)

    def __get__(self, instance, owner) -> str:
        return super().__get__(instance, owner)

    def __set__(self, instance, value: str):
        super().__set__(instance, value)


class DefaultedBoolean(DefaultedTypedValue):
    def __init__(self, default: bool = False):
        super().__init__(bool, default)

    def __get__(self, instance, owner) -> bool:
        return super().__get__(instance, owner)

    def __set__(self, instance, value: bool):
        super().__set__(instance, value)


class DefaultedInteger(DefaultedTypedValue):
    def __init__(self, default: int = None):
        super().__init__(int, default)

    def __get__(self, instance, owner) -> int:
        return super().__get__(instance, owner)

    def __set__(self, instance, value: int):
        super().__set__(instance, value)


class DefaultedFloat(DefaultedTypedValue):
    def __init__(self, default: float = None):
        super().__init__(float, default)

    def __get__(self, instance, owner) -> float:
        return super().__get__(instance, owner)

    def __set__(self, instance, value: float):
        super().__set__(instance, value)


class AbsoluteDefaultedFloat(DefaultedFloat):
    def __get__(self, instance, owner) -> float:
        return abs(super().__get__(instance, owner))


class AbsoluteDefaultedInteger(DefaultedInteger):
    def __get__(self, instance, owner) -> float:
        return abs(super().__get__(instance, owner))


class EnumeratedString(EnumeratedValue, DefaultedString):
    def __init__(self, values: [str], default: str = None):
        self.values = values
        EnumeratedValue.__init__(self, values)
        DefaultedString.__init__(self, default)

    def __set__(self, instance, value: str):
        EnumeratedValue.__set__(self, instance, value)
        DefaultedString.__set__(self, instance, value)


class EnumeratedFloat(EnumeratedValue, DefaultedFloat):
    def __init__(self, values: [float], default: float = None):
        self.values = values
        EnumeratedValue.__init__(self, values)
        DefaultedFloat.__init__(self, default)

    def __set__(self, instance, value: float):
        EnumeratedValue.__set__(self, instance, value)
        DefaultedFloat.__set__(self, instance, value)


class EnumeratedInteger(EnumeratedValue, DefaultedInteger):
    def __init__(self, values: [int], default: int = None):
        if not isinstance(values, list):
            values = list(values)
        self.values = values
        EnumeratedValue.__init__(self, values)
        DefaultedInteger.__init__(self, default)

    def __set__(self, instance, value: int):
        EnumeratedValue.__set__(self, instance, value)
        DefaultedInteger.__set__(self, instance, value)


class LateralStressGWCE(EnumeratedString):
    def __init__(self, values: [str], default: str = None):
        super().__init__(values, default)

    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if not instance.smagorinsky:
                return 'kolar_grey'
            else:
                return 'velocity_based'


class LateralStressGWCEisSymmetric(DefaultedBoolean):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if not instance.smagorinsky:
                return False
            else:
                return True


class RUNDES(DefaultedString):
    @property
    def default(self):
        return f'created on {datetime.now():%Y-%m-%d %H:%M}'


class IHOT(EnumeratedFloat):
    def __init__(self):
        super().__init__([0, 567, 568], 0)

    def __get__(self, instance, owner) -> float:
        if self.value is not None:
            return super().__get__(instance, owner)
        elif self._runtype == 'coldstart':
            return 0
        elif self._runtype == 'hotstart':
            return 567
        else:
            return self.default


class WarnElev(DefaultedFloat):
    def __get__(self, instance, owner) -> float:
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            raise NotImplementedError()


class NOLIBF(EnumeratedInteger):
    def __init__(self):
        super().__init__([0, 1, 2], 2)

    def __get__(self, instance, owner) -> int:
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            mesh_attributes = self.mesh.get_nodal_attribute_names()
            for attribute in [
                'quadratic_friction_coefficient_at_sea_floor',
                'mannings_n_at_sea_floor',
                'chezy_friction_coefficient_at_sea_floor',
            ]:
                if attribute in mesh_attributes:
                    attr = self.mesh.get_nodal_attribute(attribute)
                    if self._runtype == 'coldstart':
                        if attr['coldstart'] is True:
                            return 1
                    else:
                        if attr['hotstart'] is True:
                            return 1
            else:
                return self.default


class NTIP(EnumeratedInteger):
    def __init__(self):
        super().__init__([0, 1, 2], 1)

    def __get__(self, instance, owner) -> int:
        value = super().__get__(instance, owner)
        if value == 2:
            try:
                instance.fort24
            except AttributeError:
                raise RuntimeError('Must generate fort.24 file.')
        return value


class DTDP(AbsoluteDefaultedFloat):
    def __get__(self, instance, owner) -> float:
        if self.value is not None:
            value = super().__get__(instance, owner)
        else:
            value = self.mesh.critical_timestep(instance.CFL)
            if not instance.predictor_corrector:
                value = -value
            self.value = value
        return value

    def __set__(self, instance, value: float):
        if value == 0:
            raise ValueError('timestep cannot be 0')
        super().__set__(instance, value)


class TAUO(DefaultedFloat):
    def __get__(self, instance, owner) -> float:
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.mesh.has_nodal_attribute(
                'primitive_weighting_in_continuity_equation',
                instance._runtype):
                return -3
            elif instance.NOLIBF != 2:
                return instance.CF
            else:
                return 0.005


class FFACTOR(DefaultedString):
    def __get__(self, instance, owner) -> str:
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            value = f'{instance.CF:G}'
            if instance.NOLIBF == 2:
                value += f' {instance.HBREAK:G} {instance.FTHETA:G} {instance.FGAMMA:G}'
            return value

    def __set__(self, instance, value: str):
        try:
            value = f'{float(value):G}'
        except:
            pass
        super().__set__(instance, value)


class CF(FFACTOR):
    def __init__(self):
        super().__init__(0.0025)

    def __set__(self, instance, value: str):
        instance.FFACTOR = value
        super().__set__(instance, value)


class ESLM(AbsoluteDefaultedFloat):
    def __get__(self, instance, owner) -> float:
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.smagorinsky:
                return -instance.smagorinsky_coefficient
            else:
                return instance.horizontal_mixing_coefficient


class DRAMP(DefaultedString):
    def __get__(self, instance, owner):
        if self.value is not None:
            return f'{self.value:<.16G}{10 * " "}'
        else:
            DRAMP = instance.spinup_factor * \
                    ((instance.start_date -
                      instance.forcing_start_date).total_seconds()
                     / (60.0 * 60.0 * 24.0))
            if self.NRAMP in [0, 1]:
                return f'{DRAMP:<.16G}{10 * " "}'
            else:
                return f'{DRAMP:<.3f} ' \
                       f'{instance.DRAMPExtFlux:<.3f} ' \
                       f'{instance.FluxSettlingTime:<.3f} ' \
                       f'{instance.DRAMPIntFlux:<.3f} ' \
                       f'{instance.DRAMPElev:<.3f} ' \
                       f'{instance.DRAMPTip:<.3f} ' \
                       f'{instance.DRAMPMete:<.3f} ' \
                       f'{instance.DRAMPWRad:<.3f} ' \
                       f'{instance.DUnRampMete:<.3f} '

    def __set__(self, instance, value: str):
        try:
            value = f'{float(value):G}'
        except:
            pass
        super().__set__(instance, value)


class SLAM0(DefaultedFloat):
    def __get__(self, instance, owner) -> float:
        try:
            return super().__get__(instance, owner)
        except AttributeError:
            return np.median(instance.mesh.x)


class SFEA0(DefaultedFloat):
    def __get__(self, instance, owner) -> float:
        try:
            return super().__get__(instance, owner)
        except AttributeError:
            return np.median(instance.mesh.y)


class CORI(DefaultedFloat):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.NCOR == 0:
                raise NotImplementedError('CORI without NCOR not implemented')
            else:
                return 0.0

    def __set__(self, instance, value):
        if value is None:
            if instance.NCOR == 0:
                raise RuntimeError('Must pass CORI when NCOR=0')
            else:
                value = 0
        super.__set__(instance, value)


class NHSTAR(EnumeratedInteger):
    def __init__(self):
        super().__init__([0, 1, 2, 3, 5])

    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.forcing_start_date == instance.start_date:
                return 0
            if instance._runtype == 'coldstart':
                if instance.netcdf is True:
                    return 5
                else:
                    return 3
            else:
                return 0


class NHSINC(DefaultedInteger):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.NHSTAR == 0:
                return 0
            else:
                dt = instance.start_date - instance.forcing_start_date
                if dt.total_seconds() == 0:
                    dt = instance.end_date - instance.forcing_start_date
                    return int(
                        dt.total_seconds() / np.around(instance.DTDP, 6))
                else:
                    return int(
                        dt.total_seconds() / np.around(instance.DTDP, 6))


class THAS(DefaultedFloat):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.NFREQ > 0:
                try:
                    if instance._runtype == 'coldstart':
                        return instance.STATIM + float(instance.DRAMP)
                    else:
                        dt = instance.start_date - instance.forcing_start_date
                        return (instance.STATIM + dt.total_seconds()) / (
                            24.0 * 60.0 * 60.0)
                except TypeError:
                    #  if self.DRAMP is not castable to float()
                    raise
            else:
                return 0

    def __set__(self, instance, value):
        if value < 0:
            raise ValueError('value must be greater than or equal to 0')
        super().__set__(instance, value)


class THAF(DefaultedFloat):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            if instance.NFREQ == 0:
                return 0
            dt = instance.start_date - instance.forcing_start_date
            if instance._runtype == 'coldstart':
                if dt.total_seconds() == 0:
                    dt = instance.end_date - instance.start_date
                    return dt.total_seconds() / (24.0 * 60.0 * 60.0)
                else:
                    return dt.days
            else:
                dt = instance.start_date - instance.forcing_start_date
                return (instance.STATIM + dt.total_seconds()) / (
                    24.0 * 60.0 * 60.0)

    def __set__(self, instance, value):
        if value < 0:
            raise ValueError('value must be greater than or equal to 0')
        super().__set__(instance, value)


class NHAINC(DefaultedInteger):
    def __get__(self, instance, owner):
        if self.value is not None:
            return super().__get__(instance, owner)
        else:
            value = float('inf')
            for _output in instance._outputs:
                if _output['harmonic_analysis']:
                    if instance._runtype == 'coldstart':
                        if _output['spinup']:
                            fs = _output['sampling_rate']
                            value = np.min([value, fs.total_seconds()])
                    else:  # consider a "metonly" run?
                        fs = _output['sampling_rate']
                        value = np.min([value, fs.total_seconds()])
            if value == float('inf'):
                value = 0
            return int(value / instance.DTDP)

    def __set__(self, instance, value):
        if value < 0:
            raise ValueError('value must be greater than or equal to 0')
        super().__set__(instance, value)


class FMV(DefaultedFloat):
    def __init__(self):
        super().__init__(0)

    def __set__(self, instance, value: float):
        if value < 0 or value > 1:
            raise ValueError('value must be between 0 and 1')
        super().__set__(instance, value)


class Fort15:
    def __init__(self, mesh: AdcircMesh = None):
        self._mesh = mesh
        self._runtype = None

    @property
    def mesh(self):
        return self._mesh

    def fort15(self, runtype: str):
        self._runtype = runtype
        # ----------------
        # model options
        # ----------------
        f = []
        f.extend([
            f'{self.RUNDES:<63} ! RUNDES',
            f'{self.RUNID:<63} ! RUNID',
            f'{self.NFOVER:<63d} ! NFOVER',
            f'{self.NABOUT:<63d} ! NABOUT',
            f'{self.NSCREEN:<63d} ! NSCREEN',
            f'{self.IHOT:<63d} ! IHOT',
            f'{self.ICS:<63d} ! ICS',
            f'{self.IM:<63d} ! IM',
        ])
        if self.IM in [21, 611113]:
            f.append(f'{self.IDEN:<63d} ! IDEN')
        f.extend([
            f'{self.NOLIBF:<63G} ! NOLIBF',
            f'{self.NOLIFA:<63d} ! NOLIFA',
            f'{self.NOLICA:<63d} ! NOLICA',
            f'{self.NOLICAT:<63d} ! NOLICAT',
            f'{self.NWP:<63d} ! NWP',
        ])
        if self._runtype == 'coldstart':
            attributes = self.mesh.get_coldstart_attributes()
        elif self._runtype == 'hotstart':
            attributes = self.mesh.get_hotstart_attributes()
        f.extend(f'{attribute:<63}' for attribute in attributes)
        f.extend([
            f'{self.NCOR:<63d} ! NCOR',
            f'{self.NTIP:<63d} ! NTIP',
            f'{int((self.NRS * 100) + self.NWS):<63d} ! NWS',
            f'{self.NRAMP:<63d} ! NRAMP',
            f'{self.G:<63G} ! gravitational acceleration',
            f'{self.TAU0:<63G} ! TAU0',
        ])
        if self.TAU0 == -5:
            f.append(
                f'{self.Tau0FullDomainMin:G} '
                f'{self.Tau0FullDomainMax:G}'.ljust(63)
                + ' ! Tau0FullDomainMin Tau0FullDomainMax'
            )
        f.extend([
            f'{self.DTDP:<63.6f} ! DTDP',
            f'{self.STATIM:<63G} ! STATIM',
            f'{self.REFTIM:<63G} ! REFTIM',
        ])
        if self.NWS not in [0, 1, 9, 11]:
            interval = f'{self.WTIMINC}'
            description = 'WTIMINC - meteorological data time increment'
            if self.NRS in [1, 3, 4, 5]:
                interval += f' {self.RSTIMINC}'
                description += ', RSTIMINC wave forcing increment'
            f.append(f'{interval:<63} ! {description}')
        f.extend([
            f'{self.RNDAY:<63G} ! RNDAY',
            f'{self.DRAMP:<63} ! DRAMP',
            f'{self.A00:G} {self.B00:G} {self.C00:G}'.ljust(63)
            + ' ! A00 B00 C00',
            f'{self.H0:G} 0 0 {self.VELMIN:G}'.ljust(63) + ' ! H0 ? ? VELMIN',
            f'{self.SLAM0:G} {self.SFEA0:G}'.ljust(63) + ' ! SLAM0 SFEA0',
            f'{self.FFACTOR:<63} ! {"CF HBREAK FTHETA FGAMMA" if self.NOLIBF == 2 else "FFACTOR"}',
            f'{self.ESLM:<63G} ! {"ESL - LATERAL EDDY VISCOSITY COEFFICIENT" if not self.smagorinsky else "smagorinsky coefficient"}',
            f'{self.CORI:<63G} ! CORI',
        ])
        # ----------------
        # tidal forcings
        # ----------------
        f.append(f'{self.NTIF:<63d} ! NTIF')
        active = self._get_active_tidal_potential_constituents()
        for constituent in active:
            forcing = self.tidal_forcing(constituent)
            f.extend([
                f'{constituent}',
                f'{forcing[0]:G} {forcing[1]:G} {forcing[2]:G} {forcing[3]:G} {forcing[4]:G}',
            ])
        f.append(f'{self.NBFR:d}')
        active = self._get_active_tidal_forcing_constituents()
        for constituent in active:
            forcing = self.tidal_forcing(constituent)
            f.extend(
                [
                    f'{constituent}',
                    f'{forcing[1]:G} ' f'{forcing[3]:G} ' f'{forcing[4]:G}',
                    # f'{len(self.mesh.open_boundaries)}',
                ]
            )
        for id, bnd in self.mesh.open_boundaries.items():
            # f.append(f'{bnd["neta"]}')
            # elevation
            if bnd['iettype'] in [0, 1, 4]:
                pass
            elif bnd['iettype'] in [3, 5]:
                for constituent in self.tidal_forcing.get_active_constituents():
                    f.append(f'{constituent}')
                    vertices = self.mesh.get_xy(crs='EPSG:4326')[
                               bnd['indexes'], :]
                    amp, phase = self.tidal_forcing.tpxo(constituent, vertices)
                    f.extend(f'{amp[i]:.8e} {phase[i]:.8e}' for i in
                             range(len(vertices)))
            elif bnd['iettype'] in 2:
                bnd['iettype']['obj'].ethconst
        f.append(f'{self.ANGINN:<63G} ! ANGINN')
        # ----------------
        # other boundary forcings go here.
        # (e.g. river boundary forcing)
        # ----------------
        for id, bnd in self.mesh.open_boundaries.items():
            # f.append(f'{bnd["neta"]}')
            # velocity
            if bnd['ifltype'] in [0, 1, 4]:
                pass
            else:
                raise NotImplementedError(f'bctides generation not implemented'
                                          f' for "ifltype={bnd["ifltype"]}"')
            # temperature
            if bnd['itetype'] == 0:
                pass
            else:
                raise NotImplementedError(f'bctides generation not implemented'
                                          f' for "itetype={bnd["itetype"]}"')
            # salinity
            if bnd['isatype'] == 0:
                pass
            else:
                raise NotImplementedError(f'bctides generation not implemented'
                                          f' for "isatype={bnd["isatype"]}"')
            # tracers
            if bnd['itrtype'] == 0:
                pass
            else:
                raise NotImplementedError(f'bctides generation not implemented'
                                          f' for "itrtype={bnd["itrtype"]}"')
        # ----------------
        # output requests
        # ----------------
        # elevation out stations
        f.extend([
            f'{self.NOUTE:G} {self.TOUTSE:G} '
            f'{self.TOUTFE:G} {self.NSPOOLE:G}'.ljust(63)
            + f' ! NOUTE TOUTSE TOUTFE NSPOOLE',
            f'{self.NSTAE:<63d} ! NSTAE',
        ])
        stations = self.elevation_stations_output
        if stations['sampling_rate'] is not None:
            if self._runtype == 'coldstart':
                if stations['spinup']:
                    f.extend(
                        f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                        for station_id, (x, y) in
                        stations['collection'].items()
                    )
            else:
                f.extend(
                    f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                    for station_id, (x, y) in stations['collection'].items()
                )
        # velocity out stations
        f.extend([
            (f'{self.NOUTV:G} {self.TOUTSV:G} '
             + f'{self.TOUTFV:G} {self.NSPOOLV:G}').ljust(63)
            + ' ! NOUTV TOUTSV TOUTFV NSPOOLV'
              f'{self.NSTAV:<63G} ! NSTAV'
        ])
        stations = self.velocity_stations_output
        if stations['sampling_rate'] is not None:
            if self._runtype == 'coldstart':
                if stations['spinup']:
                    f.extend(
                        f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                        for station_id, (x, y) in
                        stations['collection'].items()
                    )
            else:
                f.extend(
                    f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                    for station_id, (x, y) in stations['collection'].items()
                )
        if self.IM == 10:
            # concentration out stations
            f.extend([
                (f'{self.NOUTC:G} {self.TOUTSC:G} '
                 + f'{self.TOUTFC:G} {self.NSPOOLC:G}').ljust(63)
                + ' ! NOUTC TOUTSC TOUTFC NSPOOLC\n',
                f'{self.NSTAC:<63d} ! NSTAC\n',
            ])
            stations = self.concentration_stations_output
            if stations['sampling_rate'] is not None:
                if self._runtype == 'coldstart':
                    if stations['spinup']:
                        f.extend(
                            f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}\n'
                            for station_id, (x, y) in
                            stations['collection'].items()
                        )
                else:
                    f.extend(
                        f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                        for station_id, (x, y) in
                        stations['collection'].items()
                    )
        if self.NWS > 0:
            # meteorological out stations
            f.extend([
                (f'{self.NOUTM:G} {self.TOUTSM:G} '
                 + f'{self.TOUTFM:G} {self.NSPOOLM:G}').ljust(63)
                + ' ! NOUTM TOUTSM TOUTFM NSPOOLM',
                f'{self.NSTAM:<63d} ! NSTAM',
            ])
            stations = self.meteorological_stations_output
            if stations['sampling_rate'] is not None:
                if stations['sampling_rate'] is not None:
                    if self._runtype == 'coldstart':
                        if stations['spinup']:
                            f.extend(
                                f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                                for station_id, (x, y) in
                                stations['collection'].items()
                            )
                    else:
                        f.extend(
                            f'{x:G} {y:G}'.ljust(63) + f' ! {station_id}'
                            for station_id, (x, y) in
                            stations['collection'].items()
                        )
        # elevation global outputs
        f.append(
            (f'{self.NOUTGE:d} {self.TOUTSGE:f} '
             + f'{self.TOUTFGE:f} {self.NSPOOLGE:d}').ljust(63)
            + ' ! NOUTGE TOUTSGE TOUTFGE NSPOOLGE'
        )
        # velocity global otuputs
        f.append(
            (f'{self.NOUTGV:d} {self.TOUTSGV:f} '
             + f'{self.TOUTFGV:f} {self.NSPOOLGV:d}').ljust(63)
            + ' ! NOUTGV TOUTSGV TOUTFGV NSPOOLGV'
        )
        if self.IM == 10:
            f.append(
                (f'{self.NOUTGC:d} {self.TOUTSGC:f} '
                 + f'{self.TOUTFGC:f} {self.NSPOOLGC:d}').ljust(63)
                + ' ! NOUTSGC TOUTGC TOUTFGC NSPOOLGC'
            )
        if self.NWS != 0:
            f.append(
                (f'{self.NOUTGM:d} {self.TOUTSGM:f} '
                 + f'{self.TOUTFGM:f} {self.NSPOOLGM:d}').ljust(63)
                + ' ! NOUTGM TOUTSGM TOUTFGM NSPOOLGM'
            )
        # harmonic analysis requests
        harmonic_analysis = False
        self._outputs = [
            self.elevation_surface_output,
            self.velocity_surface_output,
            self.elevation_stations_output,
            self.velocity_stations_output,
        ]
        for _output in self._outputs:
            if _output['harmonic_analysis']:
                if self._runtype == 'coldstart':
                    if _output['spinup']:
                        harmonic_analysis = True
                        break
                else:
                    harmonic_analysis = True
                    break
        f.append(f'{self.NFREQ:<63d} ! NFREQ')
        if harmonic_analysis:
            for constituent, forcing in self.tidal_forcing:
                f.extend([
                    f'{constituent:<63} ',
                    (f'{forcing[1]:<.16G} {forcing[3]:<.16G}'
                     + f'{forcing[4]:<.16G}').ljust(63),
                ])
        f.extend([
            (f'{self.THAS:G} {self.THAF:G} '
             + f'{self.NHAINC} {self.FMV}').ljust(63)
            + ' ! THAS THAF NHAINC FMV',
            (f'{self.NHASE:G} {self.NHASV:G} '
             + f'{self.NHAGE:G} {self.NHAGV:G}').ljust(63)
            + ' ! NHASE NHASV NHAGE NHAGV',
        ])
        # ----------------
        # hostart file generation
        # ----------------
        f.extend([
            f'{self.NHSTAR:d} {self.NHSINC:d}'.ljust(63) + ' ! NHSTAR NHSINC',
            (f'{self.ITITER:<1d} {self.ISLDIA:<1d} '
             + f'{self.CONVCR:<.15G} {self.ITMAX:<4d}').ljust(63)
            + ' ! ITITER ISLDIA CONVCR ITMAX',
        ])
        if self.vertical_mode == '3D':
            raise NotImplementedError('3D runs not yet implemented')
        f.extend([
            f'{self.NCPROJ:<63} ! NCPROJ',
            f'{self.NCINST:<63} ! NCINST',
            f'{self.NCSOUR:<63} ! NCSOUR',
            f'{self.NCHIST:<63} ! NCHIST',
            f'{self.NCREF:<63} ! NCREF',
            f'{self.NCCOM:<63} ! NCCOM',
            f'{self.NCHOST:<63} ! NCHOST',
            f'{self.NCCONV:<63} ! NCONV',
            f'{self.NCCONT:<63} ! NCCONT',
            f'{self.NCDATE:<63} ! Forcing start date / NCDATE',
        ])
        del self._outputs

        for name, namelist in self.namelists.items():
            f.append(f'&{name} ' +
                     ', '.join([f'{key}={value}'
                                for key, value in namelist.items()]) +
                     ' \\')

        return '\n'.join(f)

    def write(self, runtype: str, path: PathLike, overwrite: bool = False):
        assert runtype in ['coldstart', 'hotstart']
        fort15 = pathlib.Path(path)
        if fort15.exists() and not overwrite:
            raise Exception(f'{fort15} exists. '
                            f'Pass `overwrite=True` to overwrite.')
        with open(fort15, 'w', newline='\n') as f:
            f.write(self.fort15(runtype))

    @property
    def namelists(self) -> {str: {str: str}}:
        namelists = {}
        if self.NRS in [1, 3, 4, 5]:
            namelists['SWANOutputControl'] = {
                'SWAN_OutputHS': 'False',
                'SWAN_OutputDIR': 'False',
                'SWAN_OutputTM01': 'False',
                'SWAN_OutputTPS': 'False',
                'SWAN_OutputWIND': 'False',
                'SWAN_OutputTM02': 'False',
                'SWAN_OutputTMM10': 'False'
            }
        return namelists

    def set_time_weighting_factors_in_gcwe(self, A00, B00, C00):
        A00 = float(A00)
        B00 = float(B00)
        C00 = float(C00)
        assert A00 + B00 + C00 == 1.0, '"A00 + B00 + C00" must equal 1'
        self.__A00 = A00
        self.__B00 = B00
        self.__C00 = C00

    @staticmethod
    def parse_stations(path, station_type):
        stations = dict()
        with open(path, 'r') as f:
            for line in f:
                if station_type in line:
                    line = f.readline().split('!')[0]
                    num = int(line)
                    for i in range(num):
                        line = f.readline().split('!')
                        if len(line) > 0:
                            station_name = line[1].strip(' \n')
                        else:
                            station_name = str(i)
                        vertices = line[0].split(' ')
                        vertices = [float(x) for x in vertices if x != '']
                        x = float(vertices[0])
                        y = float(vertices[1])
                        try:
                            z = float(vertices[2])
                            stations[station_name] = (x, y, z)
                        except IndexError:
                            stations[station_name] = (x, y)
        return stations

    vertical_mode = EnumeratedString(['2D', '3D'], '2D')
    lateral_stress_in_gwce = LateralStressGWCE()
    lateral_stress_in_gwce_is_symmetrical = LateralStressGWCEisSymmetric()
    advection_in_gwce = EnumeratedString(
        ['non_conservative', 'form_1', 'form_2'], 'non_conservative')
    lateral_stress_in_momentum = EnumeratedString(
        ['velocity_based', 'flux_based'], 'velocity_based')
    lateral_stress_in_momentum_is_symmetrical = DefaultedBoolean()
    lateral_stress_in_momentum_method = EnumeratedString(
        ['2_part', 'integration_by_parts'], 'integration_by_parts')
    advection_in_momentum = EnumeratedString(
        ['non_conservative', 'form_1', 'form_2'], 'non_conservative')
    area_integration_in_momentum = EnumeratedString(
        ['corrected', 'original'], 'corrected')
    baroclinicity = DefaultedBoolean()
    smagorinsky = DefaultedBoolean()
    smagorinsky_coefficient = AbsoluteDefaultedFloat(0.2)
    gwce_solution_scheme = EnumeratedString(
        ['semi-implicit', 'explicit'], 'semi-implicit')
    horizontal_mixing_coefficient = AbsoluteDefaultedFloat(10)
    passive_scalar_transport = DefaultedBoolean()
    stress_based_3D = DefaultedBoolean()
    predictor_corrector = DefaultedBoolean(True)

    @property
    @lru_cache(maxsize=None)
    def TPXO(self):
        return TPXO()

    @property
    def timestep(self):
        return np.abs(self.DTDP)

    @timestep.setter
    def timestep(self, timestep):
        self.DTDP = timestep

    RUNDES = RUNDES()
    RUNID = DefaultedString(self.mesh.description)

    @property
    def IHOT(self):
        return self._IHOT

    _IHOT = IHOT()
    NFOVER = EnumeratedInteger([0, 1], 1)
    WarnElev = WarnElev()

    @property
    def iWarnElevDump(self):
        try:
            return self.__iWarnElevDump
        except AttributeError:
            raise NotImplementedError

    @iWarnElevDump.setter
    def iWarnElevDump(self, iWarnElevDump):
        if iWarnElevDump is not None:
            iWarnElevDump = int(iWarnElevDump)
            if iWarnElevDump not in [0, 1]:
                raise TypeError('iWarnElevDump must be 0 or 1')
            self.__iWarnElevDump = int(iWarnElevDump)
        else:
            if self.WarnElev is not None:
                raise RuntimeError(
                    'Must set iWarnElevDump if WarnElev is not ' + 'None')

    @property
    def WarnElevDumpLimit(self):
        try:
            return self.__WarnElevDumpLimit
        except AttributeError:
            raise NotImplementedError

    @WarnElevDumpLimit.setter
    def WarnElevDumpLimit(self, WarnElevDumpLimit):
        if WarnElevDumpLimit is not None:
            assert isinstance(WarnElevDumpLimit, int)
            assert WarnElevDumpLimit > 0
            self.__WarnElevDumpLimit = WarnElevDumpLimit
        else:
            if self.WarnElev is not None:
                raise RuntimeError(
                    'Must set WarnElevDumpLimit if WarnElev is ' + 'not None')

    @property
    def ErrorElev(self):
        try:
            return self.__ErrorElev
        except AttributeError:
            raise NotImplementedError

    @ErrorElev.setter
    def ErrorElev(self, ErrorElev):
        if ErrorElev is not None:
            self.__ErrorElev = float(ErrorElev)
        else:
            if self.WarnElev is not None:
                raise RuntimeError(
                    'Must set iWarnElevDump if WarnElev is not ' + 'None')

    NABOUT = EnumeratedInteger(range(-1, 4), 1)
    NSCREEN = DefaultedInteger(100)

    @property
    def NWS(self) -> int:
        """
        wind stress number
        http://adcirc.org/home/documentation/users-manual-v50/input-file-descriptions/nws-values-table/
        """

        if self._runtype == 'coldstart':
            nws = 0
        elif self.mesh is not None:
            wind_forcing = self.mesh._surface_forcing['imetype']
            if wind_forcing is not None:
                # check for wave forcing here as well.
                nws = int(wind_forcing.NWS % 100)
        else:
            nws = 0

        return nws

    @property
    def NRS(self) -> int:
        """
        radiative (wave) stress number

        100 - fort.23
        300 - SWAN
        400 - STWAVE
        500 - WW3
        """

        if self._runtype == 'coldstart':
            nrs = 0
        elif self.mesh is not None:
            wind_forcing = self.mesh._surface_forcing['imetype']
            wave_forcing = self.mesh._surface_forcing['iwrtype']
            if wave_forcing is not None:
                nrs = wave_forcing.NRS
            elif wind_forcing is not None:
                nrs = int(math.floor(wind_forcing.NWS / 100))
            else:
                nrs = 0
        else:
            nrs = 0

        return nrs

    @property
    def ICS(self):
        if self.mesh.crs.is_geographic:
            return 2
        else:
            return 1

    @property
    def IM(self):
        if self.stress_based_3D:
            return 2

        elif self.passive_scalar_transport:
            if self.vertical_mode == '2D':
                if not self.baroclinicity:
                    return 10
                else:
                    return 30
            else:
                if not self.baroclinicity:
                    return 11
                else:
                    return 31

        else:
            def get_digit_1():
                if self.vertical_mode == '2D':
                    if self.lateral_stress_in_gwce == 'kolar_grey':
                        return 1
                    elif self.lateral_stress_in_gwce == 'flux_based':
                        if self.lateral_stress_in_gwce_is_symmetrical:
                            return 4
                        else:
                            return 2
                    elif self.lateral_stress_in_gwce == 'velocity_based':
                        if self.lateral_stress_in_gwce_is_symmetrical:
                            return 5
                        else:
                            return 3
                else:
                    return 6

            def get_digit_2():
                if self.advection_in_gwce == 'non_conservative':
                    return 1
                elif self.advection_in_gwce == 'form_1':
                    return 2
                elif self.advection_in_gwce == 'form_2':
                    return 3

            def get_digit_3():
                if self.lateral_stress_in_momentum == 'velocity_based':
                    if self.lateral_stress_in_momentum_method == 'integration_by_parts':
                        if self.lateral_stress_in_momentum_is_symmetrical:
                            return 3
                        else:
                            return 1
                    else:
                        raise NotImplementedError('not implemented in adcirc')
                        return 5
                elif self.lateral_stress_in_momentum == 'flux_based':
                    if self.lateral_stress_in_momentum_method == 'integration_by_parts':
                        if self.lateral_stress_in_momentum_is_symmetrical:
                            return 4
                        else:
                            return 2
                    else:
                        raise NotImplementedError('not implemented in adcirc')
                        return 6

            def get_digit_4():
                if self.advection_in_momentum == 'non_conservative':
                    return 1
                elif self.advection_in_momentum == 'form_1':
                    return 2
                elif self.advection_in_momentum == 'form_2':
                    return 3

            def get_digit_5():
                if self.area_integration_in_momentum == 'corrected':
                    return 1
                elif self.area_integration_in_momentum == 'original':
                    return 2

            def get_digit_6():
                if (not self.baroclinicity and
                    self.gwce_solution_scheme == 'semi-implicit-legacy'):
                    return 1

                elif (not self.baroclinicity and
                      self.gwce_solution_scheme == 'explicit'):
                    return 2

                elif (not self.baroclinicity and
                      self.gwce_solution_scheme == 'semi-implicit'):
                    return 3

                else:
                    raise Exception(
                        f'No IM digit 6 for {self.baroclinicity}, '
                        f'{self.gwce_solution_scheme}')

                # elif (self.baroclinicity and
                #       self.gwce_solution_scheme == 'semi-implicit'):
                #     return 3
                # elif (self.baroclinicity and
                #       self.gwce_solution_scheme == 'explicit'):
                #     return 4

            return int(f'{get_digit_1():d}' \
                       f'{get_digit_2():d}' \
                       f'{get_digit_3():d}' \
                       f'{get_digit_4():d}' \
                       f'{get_digit_5():d}' \
                       f'{get_digit_6():d}')

    @property
    def IDEN(self):
        raise NotImplementedError('3D runs not yet supported.')
        # return self.__IDEN

    @IDEN.setter
    def IDEN(self, IDEN):
        if IDEN is not None:
            raise NotImplementedError('3D runs not yet supported.')

    NOLIBF = NOLIBF()
    NOLIFA = EnumeratedInteger([0, 1, 2], 2)
    NOLICA = EnumeratedInteger([0, 1], 1)
    NOLICAT = EnumeratedInteger([0, 1], 1)

    @property
    def NWP(self):
        if self._runtype == 'coldstart':
            return len(self.mesh.get_coldstart_attributes())
        else:
            return len(self.mesh.get_hotstart_attributes())

    @property
    def NRAMP(self):
        if self.spinup_time.total_seconds() == 0:
            return 1
        if self._runtype == 'coldstart':
            return 1
        else:
            return 8

    NCOR = EnumeratedInteger([0, 1], 1)
    NTIP = NTIP()
    CFL = DefaultedFloat(0.7)
    G = DefaultedFloat(9.81)
    DTDP = DTDP()
    TAU0 = TAUO()
    FFACTOR = FFACTOR()
    CF = CF()
    ESLM = ESLM()

    # Looks like this has always to be zero!
    # the following makes adcirc crash with not enough time
    # in meteorological inputs.
    # return (
    #     (self.start_date - self.forcing_start_date).total_seconds()
    #     / (60.*60.*24))
    STATIM = DefaultedFloat(0)
    REFTIM = DefaultedFloat(0)

    @property
    def WTIMINC(self):
        if self.NWS not in [0, 1, 9, 11]:
            return self.wind_forcing.WTIMINC
        else:
            return 0

    @property
    def RSTIMINC(self):
        if self.NRS in [1, 3, 4, 5]:
            if self.wave_forcing is not None:
                return self.wave_forcing.RSTIMINC
            else:
                return self.WTIMINC
        else:
            return 0

    @property
    def RNDAY(self):
        if self._runtype == 'coldstart':
            if self.spinup_time.total_seconds() > 0.0:
                RNDAY = self.start_date - self.forcing_start_date
            else:
                RNDAY = self.end_date - self.start_date
            return RNDAY.total_seconds() / (60.0 * 60.0 * 24.0)
        else:
            RNDAY = self.end_date - self.forcing_start_date
            return RNDAY.total_seconds() / (60.0 * 60.0 * 24.0)

    DRAMP = DRAMP()
    DRAMPExtFlux = DefaultedFloat(0)
    FluxSettlingTime = DefaultedFloat(0)
    DRAMPIntFlux = DefaultedFloat(0)

    @property
    def DRAMPElev(self):
        try:
            return self.__DRAMPElev
        except AttributeError:
            return self.spinup_factor * (
                (self.start_date - self.forcing_start_date).total_seconds()
                / (60.0 * 60.0 * 24.0)
            )

    @DRAMPElev.setter
    def DRAMPElev(self, DRAMPElev):
        self.__DRAMPElev = float(DRAMPElev)

    @property
    def DRAMPTip(self):
        try:
            return self.__DRAMPTip
        except AttributeError:
            return self.spinup_factor * (
                (self.start_date - self.forcing_start_date).total_seconds()
                / (60.0 * 60.0 * 24.0)
            )

    @DRAMPTip.setter
    def DRAMPTip(self, DRAMPTip):
        self.__DRAMPTip = float(DRAMPTip)

    DRAMPMete = DefaultedFloat(1)
    DRAMPWRad = DefaultedFloat(0)

    @property
    def DUnRampMete(self):
        try:
            return self.__DUnRampMete
        except AttributeError:
            dt = self.start_date - self.forcing_start_date
            return (self.STATIM + dt.total_seconds()) / (24.0 * 60.0 * 60.0)

    @DUnRampMete.setter
    def DUnRampMete(self, DUnRampMete):
        if DUnRampMete is None:
            DUnRampMete = self.DRAMP
        self.__DUnRampMete = float(DUnRampMete)

    @property
    def A00(self):
        try:
            return self.__A00
        except AttributeError:
            if self.gwce_solution_scheme == 'explicit':
                return 0
            if self.gwce_solution_scheme == 'semi-implicit-legacy':
                return 0.35
            if self.gwce_solution_scheme == 'semi-implicit':
                return 0.5

    @property
    def B00(self):
        try:
            return self.__B00
        except AttributeError:
            if self.gwce_solution_scheme == 'explicit':
                return 0
            if self.gwce_solution_scheme == 'semi-implicit-legacy':
                return 0.30
            if self.gwce_solution_scheme == 'semi-implicit':
                return 0.5

    @property
    def C00(self):
        try:
            return self.__C00
        except AttributeError:
            if self.gwce_solution_scheme == 'explicit':
                return 0
            if self.gwce_solution_scheme == 'semi-implicit-legacy':
                return 0.35
            if self.gwce_solution_scheme == 'semi-implicit':
                return 0

    H0 = DefaultedFloat(0.01)
    NODEDRYMIN = DefaultedInteger(0)
    NODEWETRMP = DefaultedInteger(0)
    VELMIN = DefaultedFloat(0.01)
    SLAM0 = SLAM0()
    SFEA0 = SFEA0()
    HBREAK = DefaultedFloat(1)
    FTHETA = DefaultedFloat(10)
    FGAMMA = DefaultedFloat(1 / 3)
    CORI = CORI()

    @property
    def NTIF(self):
        NTIF = 0
        if self.tidal_forcing is not None:
            for constituent in self.tidal_forcing.get_active_constituents():
                if constituent in self.tidal_forcing.major_constituents:
                    NTIF += 1
        return NTIF

    @property
    def NBFR(self):
        if self.iettype in [3, 5]:
            return self.tidal_forcing.nbfr
        return 0

    ANGINN = DefaultedFloat(110)
    NOUTE = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'stations', 'elevation'))
    TOUTSE = DefaultedInteger(
        partial(Fort15._get_TOUTS__, self, 'stations', 'elevation'))
    TOUTFE = DefaultedInteger(
        partial(Fort15._get_TOUTF__, self, 'stations', 'elevation'))
    NSPOOLE = DefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'stations', 'elevation'))
    NSTAE = DefaultedInteger(
        partial(Fort15._get_NSTA_, self, 'elevation'))
    NOUTV = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'stations', 'velocity'))
    TOUTSV = DefaultedInteger(
        partial(Fort15._get_TOUTS__, self, 'stations', 'velocity'))
    TOUTFV = DefaultedInteger(
        partial(Fort15._get_TOUTF__, self, 'stations', 'velocity'))
    NSPOOLV = DefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'stations', 'velocity'))
    NSTAV = DefaultedInteger(
        partial(Fort15._get_NSTA_, self, 'velocity'))
    NOUTM = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'stations', 'meteorological'))
    TOUTSM = DefaultedInteger(
        partial(Fort15._get_TOUTS__, self, 'stations', 'meteorological'))
    TOUTFM = DefaultedInteger(
        partial(Fort15._get_TOUTF__, self, 'stations', 'meteorological'))
    NSPOOLM = DefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'stations', 'meteorological'))
    NSTAM = DefaultedInteger(
        partial(Fort15._get_NSTA_, self, 'meteorological'))
    NOUTC = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'stations', 'concentration'))
    TOUTSC = DefaultedInteger(
        partial(Fort15._get_TOUTS__, self, 'stations', 'concentration'))
    TOUTFC = DefaultedInteger(
        partial(Fort15._get_TOUTF__, self, 'stations', 'concentration'))
    NSPOOLC = DefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'stations', 'concentration'))
    NSTAC = DefaultedInteger(
        partial(Fort15._get_NSTA__, self, 'concentration'))
    NOUTGE = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'surface', 'elevation'))
    TOUTSGE = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTS__, self, 'surface', 'elevation'))
    TOUTFGE = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTF__, self, 'surface', 'elevation'))
    NSPOOLGE = AbsoluteDefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'surface', 'elevation'))
    NOUTGV = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'surface', 'velocity'))
    TOUTSGV = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTS__, self, 'surface', 'velocity'))
    TOUTFGV = AbsoluteDefaultedFloat(
        partial(self._get_TOUTF__, self, 'surface', 'velocity'))
    NSPOOLGV = AbsoluteDefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'surface', 'velocity'))
    NOUTGM = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'surface', 'meteorological'))
    TOUTSGM = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTS__, self, 'surface', 'meteorological'))
    TOUTFGM = AbsoluteDefaultedFloat(
        partial(self._get_TOUTF__, self, 'surface', 'meteorological'))
    NSPOOLGM = AbsoluteDefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'surface', 'meteorological'))
    NOUTGC = DefaultedInteger(
        partial(Fort15._get_NOUT__, self, 'surface', 'concentration'))
    TOUTSGC = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTS__, self, 'surface', 'concentration'))
    TOUTFGC = AbsoluteDefaultedFloat(
        partial(Fort15._get_TOUTF__, self, 'surface', 'concentration'))
    NSPOOLGC = AbsoluteDefaultedInteger(
        partial(Fort15._get_NSPOOL__, self, 'surface', 'concentration'))

    @property
    def NFREQ(self):
        if self._runtype == 'coldstart':
            if np.any([_['spinup'] for _ in self._outputs]):
                if np.any([_['sampling_rate'] for _ in self._outputs]):
                    if np.any([_['harmonic_analysis'] for _ in self._outputs]):
                        return len(
                            self.tidal_forcing.get_active_constituents())
        else:
            if np.any([_['sampling_rate'] for _ in self._outputs]):
                if np.any([_['harmonic_analysis'] for _ in self._outputs]):
                    return len(self.tidal_forcing.get_active_constituents())
        return 0

    THAS = THAS()
    THAF = THAF()
    NHAINC = NHAINC()
    FMV = FMV()

    @property
    def NHASE(self):
        try:
            return self.__NHASE
        except AttributeError:
            return self._get_harmonic_analysis_state(
                self.elevation_stations_output)

    @property
    def NHASV(self):
        try:
            return self.__NHASV
        except AttributeError:
            return self._get_harmonic_analysis_state(
                self.velocity_stations_output)

    @property
    def NHAGE(self):
        try:
            return self.__NHAGE
        except AttributeError:
            return self._get_harmonic_analysis_state(
                self.elevation_surface_output)

    @property
    def NHAGV(self):
        try:
            return self.__NHAGV
        except AttributeError:
            return self._get_harmonic_analysis_state(
                self.velocity_surface_output)

    NHSTAR = NHSTAR()
    NHSINC = NHSINC()
    ITITER = EnumeratedInteger([-1, 1], 1)
    ISLDIA = EnumeratedInteger(range(6), 0)

    # https://stackoverflow.com/questions/19141432/python-numpy-machine-epsilon
    # return 500*(7./3 - 4./3 - 1)
    CONVCR = DefaultedFloat(1.0e-8)

    ITMAX = DefaultedInteger(25)
    NCPROJ = DefaultedString('')
    NCINST = DefaultedString('')
    NCSOUR = DefaultedString('')
    NCHIST = DefaultedString('')
    NCREF = DefaultedString('')
    NCCOM = DefaultedString('')
    NCHOST = DefaultedString('')
    NCCONV = DefaultedString('')
    NCCONT = DefaultedString('')

    @property
    def NCDATE(self):
        return self.forcing_start_date.strftime('%Y-%m-%d %H:%M')

    @property
    def FortranNamelists(self):
        return self.__FortranNamelists

    def _get_active_tidal_potential_constituents(self):
        if self.iettype in [3, 5]:
            return self.tidal_forcing.get_active_potential_constituents()
        else:
            return []

    def _get_active_tidal_forcing_constituents(self):
        if self.iettype in [3, 5]:
            return self.tidal_forcing.get_active_forcing_constituents()
        else:
            return []

    @property
    def elevbc(self):
        return self.mesh._boundary_forcing['iettype']['obj']

    @property
    def iettype(self):
        if self.elevbc is not None:
            return self.elevbc.iettype
        else:
            return 0

    def _get_NSTA_(self, physical_var):
        stations = self._container['stations'][physical_var]
        if self._runtype == 'coldstart':
            if stations['spinup'] is not None:
                return len(stations['collection'].keys())
            else:
                return 0
        else:
            if np.abs(self._get_NOUT__('stations', physical_var)) > 0:
                return len(stations['collection'].keys())
            else:
                return 0

    def _get_NOUT__(self, output_type, physical_var):
        output = self._container[output_type][physical_var]
        if self._runtype == 'coldstart':
            if output['spinup'] is not None:
                if output['netcdf'] is True:
                    return -5
                else:
                    return -1
            else:
                return 0

        elif self._runtype == 'hotstart':
            if output['sampling_rate'] is not None:
                if output['netcdf'] is True:
                    return -5
                else:
                    return -1
            else:
                return 0

    def _get_TOUTS__(self, output_type, physical_var):
        output = self._container[output_type][physical_var]
        # coldstart
        if self._runtype == 'coldstart':
            if output['spinup'] is not None:
                start = output['spinup_start']
            else:
                return 0

        # hotstart
        else:
            if output['sampling_rate'] is not None:
                start = output['start']
            else:
                return 0

        # typecast
        if isinstance(start, int):  # int interpreted as literal timestep
            start = timedelta(seconds=(start * self.timestep))
            if self._runtype == 'coldstart':
                start = start - self.forcing_start_date
            else:
                start = start - self.start_date

        elif isinstance(start, type(None)):
            if self._runtype == 'hotstart':
                dt = self.start_date - self.forcing_start_date
                return dt.total_seconds() / (60 * 60 * 24)
            else:
                return 0

        return start.total_seconds() / (60.0 * 60.0 * 24.0)

    def _get_TOUTF__(self, output_type, physical_var):
        output = self._container[output_type][physical_var]
        # coldstart
        if self._runtype == 'coldstart':
            if output['spinup'] is not None:
                if output['spinup_end'] is None:
                    if self.NOUTGE != 0:
                        time = self.spinup_time.total_seconds() / (
                            60.0 * 60.0 * 24.0)
                        if time > 0:
                            return time
                        else:
                            dt = self.end_date - self.start_date
                            return dt.total_seconds() / (60.0 * 60.0 * 24.0)
                    else:
                        return 0
                else:
                    raise NotImplementedError
            else:
                return 0

        # hotstart
        elif self._runtype == 'hotstart':
            if output['sampling_rate'] is not None:
                if output['end'] is None:
                    if self._runtype == 'hotstart':
                        dt = self.end_date - self.forcing_start_date
                        return dt.total_seconds() / (60 * 60 * 24)
                    # if self.NOUTGE != 0:
                    #     time = self.spinup_time.total_seconds()/(60.*60.*24.)
                    #     if time > 0:
                    #         return time
                    #     else:
                    #         dt = self.end_date - self.forcing_start_date
                    #         return dt.total_seconds()/(60.*60.*24.)
                    else:
                        dt = self.start_date - self.forcing_start_date
                        return dt.total_seconds() / (60 * 60 * 24)
                else:
                    raise NotImplementedError
            else:
                return 0

    def _get_NSPOOL__(self, output_type, physical_var):
        output = self._container[output_type][physical_var]
        if self._runtype == 'coldstart':
            if output['spinup']:
                return int(round(output['spinup'].total_seconds() / self.DTDP))
            else:
                return 0
        else:
            if output['sampling_rate'] is not None:
                if output_type == 'surface' and output[
                    'sampling_rate'].total_seconds() == 0:
                    return int((
                                   self.end_date - self.start_date).total_seconds() / self.DTDP)
                return int(round(
                    (output['sampling_rate'].total_seconds() / self.DTDP)))
            else:
                return 0

    def _get_harmonic_analysis_state(self, output):
        state = 0
        if self._runtype == 'coldstart':
            if output['spinup'] and output['harmonic_analysis']:
                if self.netcdf:
                    return 5
                else:
                    return 1
        else:
            if output['harmonic_analysis']:
                if self.netcdf:
                    return 5
                else:
                    return 1
        return state
