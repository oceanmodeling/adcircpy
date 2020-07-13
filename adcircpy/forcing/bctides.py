from functools import lru_cache
from datetime import datetime, timedelta


class Bctides:

    def __init__(self, mesh, start_date, run_time, spinup_time):
        self._mesh = mesh
        self._start_date = start_date
        self._run_time = run_time
        self._spinup_time = spinup_time

    @property
    def mesh(self):
        return self._mesh

    @property
    def bctides(self):
        f = ""
        active = self._get_active_tidal_potential_constituents()
        for constituent in active:
            forcing = self._tidal_forcing(constituent)
            f += f'{constituent} \n'
            f += f'{forcing[0]:G} '
            f += f"{forcing[1]:G} "
            f += f'{forcing[2]:G} '
            f += f'{forcing[3]:G} '
            f += f'{forcing[4]:G}'
            f += '\n'
        f += f'{self.nbfr:d}\n'
        active = self._get_active_tidal_forcing_constituents()
        for constituent in active:
            forcing = self._tidal_forcing(constituent)
            f += f'{constituent} \n'
            f += f"{forcing[2]:G} "
            f += f'{forcing[3]:G} '
            f += f'{forcing[4]:G}'
            f += '\n'
        f += f"{len(self.mesh.open_boundaries)}\n"
        for id, bnd in self.mesh.open_boundaries.items():
            f += f"{bnd['neta']} "
            f += f"{bnd['iettype']} "
            f += f"{bnd['ifltype']} "
            f += f"{bnd['itetype']} "
            f += f"{bnd['isatype']}"
            itrtype = bnd['itrtype']
            if itrtype != 0:
                f += f" {itrtype}"
            f += "\n"
            # elevation
            if bnd['iettype'] in [0, 1, 4]:
                pass
            elif bnd['iettype'] in [3, 5]:
                for constituent in self._tidal_forcing.get_active_constituents():
                    f += f'{constituent}\n'
                    vertices = self.mesh.hgrid.get_xy(
                        crs='EPSG:4326')[bnd['indexes'], :]
                    amp, phase = self._tidal_forcing.tpxo(
                        constituent, vertices)
                    for i in range(len(vertices)):
                        f += f'{amp[i]:.8e} {phase[i]:.8e}\n'
            elif bnd['iettype'] in 2:
                f += bnd['iettype']['class'].ethconst
            # velocity
            if bnd['ifltype'] in [0, 1, 4]:
                pass
            else:
                msg = "bctides generation not implemented for "
                msg += f"ifltype={bnd['ifltype']}"
                raise NotImplementedError(msg)
            # temperature
            if bnd['itetype'] == 0:
                pass
            else:
                msg = "bctides generation not implemented for "
                msg += f"itetype={bnd['itetype']}"
                raise NotImplementedError(msg)
            # salinity
            if bnd['isatype'] == 0:
                pass
            else:
                msg = "bctides generation not implemented for "
                msg += f"isatype={bnd['isatype']}"
                raise NotImplementedError(msg)
            # tracers
            if bnd['itrtype'] == 0:
                pass
            else:
                msg = "bctides generation not implemented for "
                msg += f"itrtype={bnd['itrtype']}"
                raise NotImplementedError(msg)
        return f

    @property
    @lru_cache
    def elevbc(self):
        # the key itself is the forcing object
        return self.mesh.get_boundary_forcing('iettype')["class"]

    @property
    def start_date(self):
        return self._start_date

    @property
    def spinup_time(self):
        return self._spinup_time

    @property
    def run_time(self):
        return self._run_time

    @property
    def forcing_start_date(self):
        return self.start_date - self.spinup_time

    @property
    def end_date(self):
        return self.start_date + self.run_time

    @property
    def iettype(self):
        return self.elevbc.iettype

    @property
    def ntip(self):
        if self.iettype in [3, 5]:
            return self._tidal_forcing.ntip
        return 0

    @property
    def ntif(self):
        NTIF = 0
        for constituent in self._tidal_forcing.get_active_constituents():
            if constituent in self._tidal_forcing.major_constituents:
                NTIF += 1
        return NTIF

    @property
    def nbfr(self):
        if self.iettype in [3, 5]:
            return self._tidal_forcing.nbfr
        return 0

    @property
    def cutoff_depth(self):
        if self.iettype in [3, 5]:
            return self._tidal_forcing.cutoff_depth
        return 0

    def _get_active_tidal_potential_constituents(self):
        if self.iettype in [3, 5]:
            return self._tidal_forcing.get_active_potential_constituents()
        else:
            return []

    def _get_active_tidal_forcing_constituents(self):
        if self.iettype in [3, 5]:
            return self._tidal_forcing.get_active_forcing_constituents()
        else:
            return []

    @property
    def _mesh(self):
        return self.__mesh

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _run_time(self):
        return self.__run_time

    @property
    def _spinup_time(self):
        return self.__spinup_time

    @property
    @lru_cache
    def _tidal_forcing(self):
        elevbc = self.elevbc
        elevbc.start_date = self.start_date
        elevbc.end_date = self.end_date
        elevbc.spinup_time = self.spinup_time
        return elevbc

    @_mesh.setter
    def _mesh(self, mesh):
        from adcircpy import AdcircMesh
        assert isinstance(mesh, AdcircMesh)
        self.__mesh = mesh

    @_start_date.setter
    def _start_date(self, start_date):
        assert isinstance(start_date, datetime)
        self.__start_date = start_date

    @_run_time.setter
    def _run_time(self, run_time):
        assert isinstance(run_time, timedelta)
        self.__run_time = run_time

    @_spinup_time.setter
    def _spinup_time(self, spinup_time):
        assert isinstance(spinup_time, timedelta)
        self.__spinup_time = spinup_time
