import numpy as np
from datetime import datetime
from adcircpy.lib._TPXO import _TPXO


class _Fort15(object):

    def __init__(self):
        self.__TPXO = _TPXO()

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

    def __get_output_stations(self, output_type):
        if output_type == 'elevation':
            return self.get_elevation_output_stations()
        elif output_type == 'velocity':
            return self.get_velocity_output_stations()
        elif output_type == 'concentration':
            return self.get_concentration_output_stations()
        elif output_type == 'meteorological':
            return self.get_meteorological_output_stations()

    def __get_global_output_request(self, output_type):
        if output_type == 'elevation':
            return self.get_elevation_global_output()
        elif output_type == 'velocity':
            return self.get_velocity_global_output()
        elif output_type == 'concentration':
            return self.get_concentration_global_output()
        elif output_type == 'meteorological':
            return self.get_meteorological_global_output()

    def __get_NOUT_(self, output_type):
        stations = self.__get_output_stations(output_type)
        if stations['sampling_frequency'].total_seconds() > 0.:
            if self.runtype == 'coldstart':
                if stations['spinup']:
                    if stations['netcdf'] is True:
                        return -3
                    else:
                        return 3
                else:
                    return 0
            else:
                if stations['netcdf'] is True:
                    return -3
                else:
                    return 3
        return 0

    def __get_TOUTS_(self, _):
        if self.runtype == 'coldstart':
            return 0
        else:
            return self.spinup_time.total_seconds()/(60.*60.*24.)

    def __get_TOUTF_(self, output_type):
        stations = self.__get_output_stations(output_type)
        if stations['sampling_frequency'].total_seconds() > 0.:
            if self.runtype == 'coldstart':
                if stations['spinup']:
                    return self.spinup_time.total_seconds()/(60.*60.*24.)
            else:
                return (
                    (self.end_date - self.forcing_start_date).total_seconds()
                    / (60.*60.*24.))
        return 0

    def __get_NSPOOL_(self, output_type):
        stations = self.__get_output_stations(output_type)
        if self.runtype == 'coldstart':
            if stations['spinup'] is True:
                return int((stations['sampling_frequency'].total_seconds()
                            / self.DTDP))
            else:
                return 0
        else:
            if stations['sampling_frequency'].total_seconds() > 0:
                return int((stations['sampling_frequency'].total_seconds()
                            / self.DTDP))
            else:
                return 0

    def __get_NSTA_(self, output_type):
        stations = self.__get_output_stations(output_type)
        if self.runtype == 'coldstart':
            if stations['spinup'] is True:
                return len(stations['station_id'].keys())
            else:
                return 0
        else:
            return len(stations['station_id'].keys())

    def __get_NOUTG_(self, output_type):
        stations = self.__get_global_output_request(output_type)
        if stations['sampling_frequency'].total_seconds() > 0.:
            if self.runtype == 'coldstart':
                if stations['spinup']:
                    if stations['netcdf'] is True:
                        return -3
                    else:
                        return 3
                else:
                    return 0
            else:
                if stations['netcdf'] is True:
                    return -3
                else:
                    return 3
        return 0

    def __get_TOUTSG_(self, output_type):
        if self.runtype == 'coldstart':
            return 0
        else:
            return self.spinup_time.total_seconds()/(60.*60.*24.)

    def __get_TOUTFG_(self, output_type):
        if self.runtype == 'coldstart':
            return self.spinup_time.total_seconds()/(60.*60.*24.)
        else:
            return (
                (self.end_date - self.forcing_start_date).total_seconds()
                / (60.*60.*24.))

    def __get_NSPOOLG_(self, output_type):
        request = self.__get_global_output_request(output_type)
        if self.runtype == 'coldstart':
            if request['spinup'] is True:
                return int((request['sampling_frequency'].total_seconds()
                            / self.DTDP))
            else:
                return 0
        else:
            if request['sampling_frequency'].total_seconds() > 0:
                return int((request['sampling_frequency'].total_seconds()
                            / self.DTDP))
            else:
                return 0

    def __get_harmonic_analysis_state(self, output):
        state = 0
        if self.runtype == 'coldstart':
            if output['spinup'] and output['harmonic_analysis']:
                return 1
        else:
            if output['harmonic_analysis']:
                return 1
        return state

    @property
    def fort15(self):
        # ----------------
        # model options
        # ----------------
        f = '{:<32}! RUNDES\n'.format(self.RUNDES)
        f += '{:<32}! RUNID\n'.format(self.RUNID)
        f += '{:<32d}! NFOVER\n'.format(self.NFOVER)
        f += '{:<32d}! NABOUT\n'.format(self.NABOUT)
        f += '{:<32d}! NSCREEN\n'.format(self.NSCREEN)
        f += '{:<32d}! IHOT\n'.format(self.IHOT)
        f += '{:<32d}! ICS\n'.format(self.ICS)
        f += '{:<32d}! IM\n'.format(self.IM)
        if self.IM == 21:
            f += '{}\n'.format(self.IDEN)
        f += '{:<32d}! NOLIBF\n'.format(self.NOLIBF)
        f += '{:<32d}! NOLIFA\n'.format(self.NOLIFA)
        f += '{:<32d}! NOLICA\n'.format(self.NOLICA)
        f += '{:<32d}! NOLICAT\n'.format(self.NOLICAT)
        f += '{:<32d}! NWP\n'.format(self.NWP)
        if self.runtype == 'coldstart':
            attributes = self.AdcircMesh.get_coldstart_attributes()
        elif self.runtype == 'hotstart':
            attributes = self.AdcircMesh.get_hotstart_attributes()
        for attribute in attributes.keys():
            f += '{:<32}\n'.format(attribute)
        f += '{:<32d}! NCOR\n'.format(self.NCOR)
        f += '{:<32d}! NTIP\n'.format(self.NTIP)
        f += '{:<32d}! NWS\n'.format(self.NWS)
        f += '{:<32d}! NRAMP\n'.format(self.NRAMP)
        f += '{:<32.16E}! G\n'.format(self.G)
        if isinstance(self.TAU0, int):
            f += '{:<32d}! TAU0\n'.format(self.TAU0)
        else:
            f += '{:<32.16E}! TAU0\n'.format(self.TAU0)
        if self.TAU0 == -5:
            f += '{} {}\n'.format(self.Tau0FullDomainMin,
                                  self.Tau0FullDomainMax)
        f += '{:<32.6E}! DTDP\n'.format(self.DTDP)
        f += '{:<32.6E}! STATIM\n'.format(self.STATIM)
        f += '{:<32.6E}! REFTIM\n'.format(self.REFTIM)
        if self.NWS not in [0, 1, 9, 11]:
            f += '{}'.format(self.WTIMINC)
            f += 12*' '
            f += '! WTIMINC\n'
        f += '{:<32.6E}! RNDAY\n'.format(self.RNDAY)
        f += '{:<32}! DRAMP\n'.format(self.DRAMP)
        f += '{:<.3E} {:<.3E} {:<.3E}'.format(self.A00, self.B00, self.C00)
        f += '   ! A00 B00 C00\n'
        f += '{:<.6E} 0 0 {:<.6E}'.format(self.H0, self.VELMIN)
        f += '   ! H0 NODEDRYMIN NODEWETRMP VELMIN\n'
        f += '{:<.6E} {:<.6E}'.format(self.SLAM0, self.SFEA0)
        f += '      ! SLAM0 SFEA0\n'
        if self.NOLIBF in [0, 1]:
            f += '{:<32.4f}! FFACTOR\n'.format(self.FFACTOR)
        elif self.NOLIBF == 2:
            raise NotImplementedError('Need to setup for NOLIBF==2')
        if self.IM in [0, 1, 2]:
            f += '{:<32.4f}! ESLM (eddy viscosity)\n'.format(self.ESLM)
        elif self.IM == 10:
            raise NotImplementedError('IM==10')
        f += '{:<32.4f}! CORI\n'.format(self.CORI)
        # ----------------
        # tidal forcings
        # ----------------
        f += '{:<32d}! NTIF\n'.format(self.NTIF)
        for constituent, forcing in self.TidalForcing:
            if constituent in self.TidalForcing.major_constituents:
                f += '{:<32}\n'.format(constituent)
                f += '{:<.8E} '.format(forcing[0])
                f += '{:<.8E} '.format(forcing[1])
                f += '{:<.8E} '.format(forcing[2])
                f += '{:<.8E} '.format(forcing[3])
                f += '{:<.8E} '.format(forcing[4])
                f += '\n'
        f += '{:<32d}! NBFR\n'.format(len(self.TidalForcing))
        for constituent, forcing in self.TidalForcing:
            f += '{:<32}\n'.format(constituent)
            f += '{:<.8E} '.format(forcing[1])
            f += '{:<.8E} '.format(forcing[3])
            f += '{:<.8E} '.format(forcing[4])
            f += '\n'
        # NOTE:  This part is written as one-constituent then all boundaries
        # as opposed to one-boundary then all constituents for that boundary.
        # Not exactly sure how ADCIRC handles multiple open boundaries.
        for constituent in self.TidalForcing.get_active_constituents():
            vertices = self.AdcircMesh.xy[self.AdcircMesh.ocean_boundary, :]
            amp, phase = self.TPXO(constituent, vertices)
            f += '{}\n'.format(constituent)
            for i in range(len(vertices)):
                f += '{:<.16E} '.format(amp[i])
                f += '{:<.16E}\n'.format(phase[i])
        f += "{:<32.2f}! ANGINN\n".format(self.ANGINN)
        # ----------------
        # other boundary type forcing go here.
        # (e.g. river boundary forcing)
        # ----------------
        # ----------------
        # output requests
        # ----------------
        # elevation out stations
        f += "{:<2d} ".format(self.NOUTE)
        f += "{:<.4E} ".format(self.TOUTSE)
        f += "{:<.4E} ".format(self.TOUTFE)
        f += "{:<7d}".format(self.NSPOOLE)
        f += "! NOUTE TOUTSE TOUTFE NSPOOLE\n"
        f += "{:<32d}! NSTAE\n".format(self.NSTAE)
        if self.NSTAE > 0:
            stations = self.get_elevation_output_stations()
            for station_id, (x, y) in stations['station_id'].items():
                f += "{:<10.4f} ".format(x)
                f += "{:<10.4f} ".format(y)
                f += "{:10}".format('')
                f += "! {}\n".format(station_id)
        # velocity out stations
        f += "{:<2d} ".format(self.NOUTV)
        f += "{:<.4E} ".format(self.TOUTSV)
        f += "{:<.4E} ".format(self.TOUTFV)
        f += "{:<7d}".format(self.NSPOOLV)
        f += "! NOUTV TOUTSV TOUTFV NSPOOLV\n"
        f += "{:<32d}! NSTAV\n".format(self.NSTAV)
        if self.NSTAV > 0:
            stations = self.get_velocity_output_stations()
            for station_id, (x, y) in stations['station_id'].items():
                f += "{:<10.4f} ".format(x)
                f += "{:<10.4f} ".format(y)
                f += "{:20}".format('')
                f += "! {}\n".format(station_id)
        if self.IM == 10:
            # concentration out stations
            f += "{:<2d} ".format(self.NOUTC)
            f += "{:<.4E} ".format(self.TOUTSC)
            f += "{:<.4E} ".format(self.TOUTFC)
            f += "{:<7d}\n".format(self.NSPOOLC)
            f += "{:<32d}! NSTAV\n".format(self.NSTAC)
            if self.NSTAC > 0:
                stations = self.get_concentration_output_stations()
                for station_id, (x, y) in stations['station_id'].items():
                    f += "{:<10.4f} ".format(x)
                    f += "{:<10.4f} ".format(y)
                    f += "{:20}".format('')
                    f += "! {}\n".format(station_id)
        if self.NWS > 0:
            # meteorological out stations
            f += "{:<2d} ".format(self.NOUTM)
            f += "{:<.4E} ".format(self.TOUTSM)
            f += "{:<.4E} ".format(self.TOUTFM)
            f += "{:<7d}\n".format(self.NSPOOLM)
            f += "{:<32d}! NSTAM\n".format(self.NSTAM)
            if self.NSTAM > 0:
                stations = self.get_meteorological_output_stations()
                for station_id, (x, y) in stations['station_id'].items():
                    f += "{:<10.4f} ".format(x)
                    f += "{:<10.4f} ".format(y)
                    f += "{:20}".format('')
                    f += "! {}\n".format(station_id)
        # elevation global outputs
        f += "{:<2d} ".format(self.NOUTGE)
        f += "{:<.4E} ".format(self.TOUTSGE)
        f += "{:<.4E} ".format(self.TOUTFGE)
        f += "{:<7d}".format(self.NSPOOLGE)
        f += "! NOUTGE TOUTSGE TOUTFGE NSPOOLGE\n"
        # velocity global otuputs
        f += "{:<2d} ".format(self.NOUTGV)
        f += "{:<.4E} ".format(self.TOUTSGV)
        f += "{:<.4E} ".format(self.TOUTFGV)
        f += "{:<7d}".format(self.NSPOOLGV)
        f += "! NOUTGV TOUTSGV TOUTFGV NSPOOLGV\n"
        if self.IM == 10:
            f += "{:<2d} ".format(self.NOUTGC)
            f += "{:<.4E} ".format(self.TOUTSGC)
            f += "{:<.4E} ".format(self.TOUTFGC)
            f += "{:<7d}".format(self.NSPOOLGC)
            f += "! NOUTGC TOUTSGC TOUTFGC NSPOOLGC\n"
        if self.NWS != 0:
            f += "{:<2d} ".format(self.NOUTGM)
            f += "{:<.4E} ".format(self.TOUTSGM)
            f += "{:<.4E} ".format(self.TOUTFGM)
            f += "{:<7d}".format(self.NSPOOLGM)
            f += "! NOUTGM TOUTSGM TOUTFGM NSPOOLGM\n"
        # harmonic analysis requests
        harmonic_analysis = False
        self._outputs = [self.get_elevation_global_output(),
                         self.get_velocity_global_output(),
                         self.get_elevation_output_stations(),
                         self.get_velocity_output_stations()]
        for _output in self._outputs:
            if _output['harmonic_analysis']:
                if self.runtype == 'coldstart':
                    if _output['spinup']:
                        harmonic_analysis = True
                        break
                else:
                    harmonic_analysis = True
                    break
        if harmonic_analysis:
            f += '{:<32d}! NFREQ\n'.format(self.NFREQ)
            for constituent, forcing in self.TidalForcing:
                f += '{:<32}\n'.format(constituent)
                f += '{:<.16E} '.format(forcing[1])
                f += '{:<.16E} '.format(forcing[3])
                f += '{:<.16E} '.format(forcing[4])
                f += '\n'
        else:
            f += '{:<32d}! NFREQ\n'.format(self.NFREQ)
        f += '{:<.3f} '.format(self.THAS)
        f += '{:<.3f} '.format(self.THAF)
        f += '{:<10d} '.format(self.NHAINC)
        f += '{:<.3f} '.format(self.FMV)
        f += 3*' '
        f += '! THAS THAF NHAINC FMV\n'
        f += "{:<2d} ".format(self.NHASE)
        f += "{:<2d} ".format(self.NHASV)
        f += "{:<2d} ".format(self.NHAGE)
        f += "{:<2d}".format(self.NHAGV)
        f += 21*' '
        f += "! NHASE NHASV NHAGE NHAGV\n"
        # ----------------
        # hostart file generation
        # ----------------
        f += "{:<2d} ".format(self.NHSTAR)
        f += "{:<10d} ".format(self.NHSINC)
        f += 18*' '
        f += "! NHSTAR NHSINC\n"
        f += "{:<1d} ".format(self.ITITER)
        f += "{:<1d} ".format(self.ISLDIA)
        f += "{:<.15E}  ".format(self.CONVCR)
        f += "{:<4d} ".format(self.ITMAX)
        f += "! ITITER ISLDIA CONVCR ITMAX\n"
        if self.IM not in [0, 111112]:
            raise NotImplementedError('3d runs not yet implemented')
        f += "{}\n".format(self.NCPROJ)
        f += "{}\n".format(self.NCINST)
        f += "{}\n".format(self.NCSOUR)
        f += "{}\n".format(self.NCHIST)
        f += "{}\n".format(self.NCREF)
        f += "{}\n".format(self.NCCOM)
        f += "{}\n".format(self.NCHOST)
        f += "{}\n".format(self.NCCONV)
        f += "{}\n".format(self.NCCONT)
        f += "{}".format(self.NCDATE)
        f += 16*" "
        f += "! NCDATE (forcing start date)\n"
        del self._outputs
        return f

    @property
    def runtype(self):
        try:
            return self.__runtype
        except AttributeError:
            raise AttributeError('Must set runtype attribute.')

    @property
    def RUNDES(self):
        try:
            self.__RUNDES
        except AttributeError:
            return datetime.now().strftime('created on %Y-%m-%d %H:%M')

    @property
    def RUNID(self):
        try:
            self.__RUNID
        except AttributeError:
            return self.AdcircMesh.description

    @property
    def IHOT(self):
        try:
            return self.__IHOT
        except AttributeError:
            if self.runtype == 'coldstart':
                return 0
            elif self.runtype == 'hotstart':
                return 567

    @property
    def NFOVER(self):
        try:
            self.__NFOVER
        except AttributeError:
            return 1

    @property
    def WarnElev(self):
        try:
            return self.__WarnElev
        except AttributeError:
            raise NotImplementedError

    @property
    def iWarnElevDump(self):
        try:
            return self.__iWarnElevDump
        except AttributeError:
            raise NotImplementedError

    @property
    def WarnElevDumpLimit(self):
        try:
            return self.__WarnElevDumpLimit
        except AttributeError:
            raise NotImplementedError

    @property
    def ErrorElev(self):
        try:
            return self.__ErrorElev
        except AttributeError:
            raise NotImplementedError

    @property
    def NABOUT(self):
        try:
            self.__NABOUT
        except AttributeError:
            return 1

    @property
    def NSCREEN(self):
        try:
            return self.__NSCREEN
        except AttributeError:
            return 100

    @property
    def NWS(self):
        if self.runtype == 'coldstart':
            return 0
        else:
            if self.WindForcing is not None:
                # check for wave forcing here as well.
                return self.WindForcing.NWS
            else:
                return 0

    @property
    def ICS(self):
        if self.AdcircMesh.SpatialReference.IsGeographic():
            return 2
        else:
            return 1

    @property
    def IM(self):
        try:
            return self.__IM
        except AttributeError:
            return 0

    @property
    def IDEN(self):
        return self.__IDEN

    @property
    def NOLIBF(self):
        try:
            return self.__NOLIBF
        except AttributeError:
            NOLIBF = 0
            for attribute in [
                    'quadratic_friction_coefficient_at_sea_floor',
                    'mannings_n_at_sea_floor',
                    'chezy_friction_coefficient_at_sea_floor']:
                try:
                    attr = self.AdcircMesh.get_nodal_attribute(attribute)
                    if self.runtype == 'coldstart':
                        if attr['coldstart'] is True:
                            NOLIBF = 1
                    else:
                        if attr['hotstart'] is True:
                            NOLIBF = 1
                except AssertionError:
                    pass
            return NOLIBF

    @property
    def NOLIFA(self):
        try:
            return self.__NOLIFA
        except AttributeError:
            return 2

    @property
    def NOLICA(self):
        try:
            return self.__NOLICA
        except AttributeError:
            return 1

    @property
    def NOLICAT(self):
        try:
            return self.__NOLICAT
        except AttributeError:
            return 1

    @property
    def NWP(self):
        if self.runtype == 'coldstart':
            return len(self.AdcircMesh.get_coldstart_attributes())
        else:
            return len(self.AdcircMesh.get_hotstart_attributes())

    @property
    def NRAMP(self):
        if self.runtype == 'coldstart':
            return 1
        else:
            return 8

    @property
    def NCOR(self):
        try:
            return self.__NCOR
        except AttributeError:
            return 1

    @property
    def NTIP(self):
        try:
            return self.__NTIP
        except AttributeError:
            return 1

    @property
    def G(self):
        try:
            return self.__G
        except AttributeError:
            return 9.81

    @property
    def DTDP(self):
        try:
            return self.__DTDP
        except AttributeError:
            raise Exception('Must set DTDP attribute.')

    @property
    def TAU0(self):
        try:
            return self.__TAU0
        except AttributeError:
            TAU0 = self.AdcircMesh.has_primitive_weighting(self.runtype)
            if TAU0:
                return -3
            else:
                return self.FFACTOR

    @property
    def FFACTOR(self):
        try:
            return self.__FFACTOR
        except AttributeError:
            return 0.0025

    @property
    def ESLM(self):
        try:
            return self.__ESLM
        except AttributeError:
            return 10.

    @property
    def STATIM(self):
        try:
            return self.__STATIM
        except AttributeError:
            if self.runtype == 'coldstart':
                return 0
            else:
                return (
                    (self.start_date - self.forcing_start_date).total_seconds()
                    / (60.*60.*24))

    @property
    def REFTIM(self):
        try:
            return self.__REFTIM
        except AttributeError:
            return 0.

    @property
    def WTIMINC(self):
        if self.NWS not in [0, 1, 9, 11]:
            return self.WindForcing.WTIMINC
        else:
            return 0

    @property
    def RNDAY(self):
        if self.spinup_time.total_seconds() > 0.:
            RNDAY = self.start_date - self.forcing_start_date
        else:
            RNDAY = self.end_date - self.forcing_start_date
        return RNDAY.total_seconds()/(60.*60.*24.)

    @property
    def DRAMP(self):
        DRAMP = str(self.spinup_factor * (
                (self.start_date - self.forcing_start_date).total_seconds()
                / (60.*60.*24.)))
        if self.NRAMP in [0, 1]:
            return DRAMP
        else:
            DRAMP += '{} '.format(self.DRAMPExtFlux)
            DRAMP += '{} '.format(self.FluxSettlingTime)
            DRAMP += '{} '.format(self.DRAMPIntFlux)
            DRAMP += '{} '.format(self.DRAMPElev)
            DRAMP += '{} '.format(self.DRAMPTip)
            DRAMP += '{} '.format(self.DRAMPMete)
            DRAMP += '{} '.format(self.DRAMPWRad)
            DRAMP += '{}'.format(self.DUnRampMete)
            return DRAMP

    @property
    def DRAMPExtFlux(self):
        try:
            return self.__DRAMPExtFlux
        except AttributeError:
            return 0.

    @property
    def FluxSettlingTime(self):
        try:
            return self.__FluxSettlingTime
        except AttributeError:
            return 0.

    @property
    def DRAMPIntFlux(self):
        try:
            return self.__DRAMPIntFlux
        except AttributeError:
            return 0.

    @property
    def DRAMPElev(self):
        try:
            return self.__DRAMPElev
        except AttributeError:
            return self.spinup_factor * (
                (self.start_date - self.forcing_start_date).total_seconds()
                / (60.*60.*24.))

    @property
    def DRAMPTip(self):
        try:
            return self.__DRAMPTip
        except AttributeError:
            return self.spinup_factor * (
                (self.start_date - self.forcing_start_date).total_seconds()
                / (60.*60.*24.))

    @property
    def DRAMPMete(self):
        try:
            return self.__DRAMPMete
        except AttributeError:
            return 1.

    @property
    def DRAMPWRad(self):
        try:
            return self.__DRAMPWRad
        except AttributeError:
            return 0.

    @property
    def DUnRampMete(self):
        try:
            return self.__DUnRampMete
        except AttributeError:
            dt = self.start_date - self.forcing_start_date
            return (self.STATIM + dt.total_seconds()) / (24.*60.*60.)

    @property
    def A00(self):
        try:
            return self.__A00
        except AttributeError:
            return 0.35

    @property
    def B00(self):
        try:
            return self.__B00
        except AttributeError:
            return 0.30

    @property
    def C00(self):
        try:
            return self.__C00
        except AttributeError:
            return 0.35

    @property
    def H0(self):
        try:
            return self.__H0
        except AttributeError:
            return 0.05

    @property
    def NODEDRYMIN(self):
        try:
            return self.__NODEDRYMIN
        except AttributeError:
            return 0

    @property
    def NODEWETRMP(self):
        try:
            return self.__NODEWETRMP
        except AttributeError:
            return 0

    @property
    def VELMIN(self):
        try:
            return self.__VELMIN
        except AttributeError:
            return 0.05

    @property
    def SLAM0(self):
        try:
            return self.__SLAM0
        except AttributeError:
            # The meaning of SLAM0 is not clear from documentation.
            # approximating the number as the median of the x component.
            return np.median(self.AdcircMesh.x)

    @property
    def SFEA0(self):
        try:
            return self.__SFEA0
        except AttributeError:
            # The meaning of SFEA0 is not clear from documentation.
            # approximating the number as the median of the y component.
            return np.median(self.AdcircMesh.y)

    @property
    def TAU(self):
        return self.__TAU

    @property
    def HBREAK(self):
        return self.__HBREAK

    @property
    def FTHETA(self):
        return self.__FTHETA

    @property
    def FGAMMA(self):
        return self.__FGAMMA

    @property
    def ESLC(self):
        return self.__ESLC

    @property
    def CORI(self):
        try:
            return self.__CORI
        except AttributeError:
            if self.NCOR == 0:
                raise NotImplementedError
            else:
                return 0.

    @property
    def NTIF(self):
        if self.TidalForcing is not None:
            NTIF = 0
            for constituent in self.TidalForcing.get_active_constituents():
                if constituent in self.TidalForcing.major_constituents:
                    NTIF += 1
            return NTIF
        else:
            return 0

    @property
    def ANGINN(self):
        try:
            return self.__ANGINN
        except AttributeError:
            return 110.

    @property
    def NOUTE(self):
        # netcdf, ascii or nothing
        try:
            self.__NOUTE
        except AttributeError:
            return self.__get_NOUT_('elevation')

    @property
    def TOUTSE(self):
        # number of days relative to coldstart after which data is recorded
        try:
            return self.__TOUTSE
        except AttributeError:
            return self.__get_TOUTS_('elevation')

    @property
    def TOUTFE(self):
        try:
            return self.__TOUTFE
        except AttributeError:
            return self.__get_TOUTF_('elevation')

    @property
    def NSPOOLE(self):
        # sampling frequency
        try:
            return self.__NSPOOLE
        except AttributeError:
            return self.__get_NSPOOL_('elevation')

    @property
    def NSTAE(self):
        try:
            return self.__NSTAE
        except AttributeError:
            return self.__get_NSTA_('elevation')

    @property
    def NOUTV(self):
        # netcdf, ascii or nothing
        try:
            self.__NOUTV
        except AttributeError:
            return self.__get_NOUT_('velocity')

    @property
    def TOUTSV(self):
        # number of days relative to coldstart after which data is recorded
        try:
            return self.__TOUTSV
        except AttributeError:
            return self.__get_TOUTS_('velocity')

    @property
    def TOUTFV(self):
        try:
            return self.__TOUTFV
        except AttributeError:
            return self.__get_TOUTF_('velocity')

    @property
    def NSPOOLV(self):
        # sampling frequency
        try:
            return self.__NSPOOLV
        except AttributeError:
            return self.__get_NSPOOL_('velocity')

    @property
    def NSTAV(self):
        try:
            return self.__NSTAV
        except AttributeError:
            return self.__get_NSTA_('velocity')

    @property
    def NOUTM(self):
        # netcdf, ascii or nothing
        try:
            self.__NOUTM
        except AttributeError:
            return self.__get_NOUT_('meteorological')

    @property
    def TOUTSM(self):
        # number of days relative to coldstart after which data is recorded
        try:
            return self.__TOUTSM
        except AttributeError:
            return self.__get_TOUTS_('meteorological')

    @property
    def TOUTFM(self):
        try:
            return self.__TOUTFM
        except AttributeError:
            return self.__get_TOUTF_('meteorological')

    @property
    def NSPOOLM(self):
        # sampling frequency
        try:
            return self.__NSPOOLM
        except AttributeError:
            return self.__get_NSPOOL_('meteorological')

    @property
    def NSTAM(self):
        try:
            return self.__NSTAM
        except AttributeError:
            return self.__get_NSTA_('meteorological')

    @property
    def NOUTC(self):
        # netcdf, ascii or nothing
        try:
            self.__NOUTC
        except AttributeError:
            return self.__get_NOUT_('concentration')

    @property
    def TOUTSC(self):
        # number of days relative to coldstart after which data is recorded
        try:
            return self.__TOUTSC
        except AttributeError:
            return self.__get_TOUTS_('concentration')

    @property
    def TOUTFC(self):
        try:
            return self.__TOUTFC
        except AttributeError:
            return self.__get_TOUTF_('concentration')

    @property
    def NSPOOLC(self):
        # sampling frequency
        try:
            return self.__NSPOOLC
        except AttributeError:
            return self.__get_NSPOOL_('concentration')

    @property
    def NSTAC(self):
        try:
            return self.__NSTAC
        except AttributeError:
            return self.__get_NSTA_('concentration')

    @property
    def NOUTGE(self):
        try:
            return self.__NOUTGE
        except AttributeError:
            return self.__get_NOUTG_('elevation')

    @property
    def TOUTSGE(self):
        try:
            return self.__TOUTSGE
        except AttributeError:
            return self.__get_TOUTSG_('elevation')

    @property
    def TOUTFGE(self):
        try:
            return self.__TOUTFGE
        except AttributeError:
            if self.NOUTGE == 0:
                return 0.
            else:
                return self.__get_TOUTFG_('elevation')

    @property
    def NSPOOLGE(self):
        try:
            return self.__NSPOOLGE
        except AttributeError:
            return self.__get_NSPOOLG_('elevation')

    @property
    def NOUTGV(self):
        try:
            return self.__NOUTGV
        except AttributeError:
            return self.__get_NOUTG_('velocity')

    @property
    def TOUTSGV(self):
        try:
            return self.__TOUTSGV
        except AttributeError:
            return self.__get_TOUTSG_('velocity')

    @property
    def TOUTFGV(self):
        try:
            return self.__TOUTFGV
        except AttributeError:
            if self.NOUTGV == 0:
                return 0.
            else:
                return self.__get_TOUTFG_('velocity')

    @property
    def NSPOOLGV(self):
        try:
            return self.__NSPOOLGV
        except AttributeError:
            return self.__get_NSPOOLG_('velocity')

    @property
    def NOUTGM(self):
        try:
            return self.__NOUTGM
        except AttributeError:
            return self.__get_NOUTG_('meteorological')

    @property
    def TOUTSGM(self):
        try:
            return self.__TOUTSGM
        except AttributeError:
            return self.__get_TOUTSG_('meteorological')

    @property
    def TOUTFGM(self):
        try:
            return self.__TOUTFGM
        except AttributeError:
            return self.__get_TOUTFG_('meteorological')

    @property
    def NSPOOLGM(self):
        try:
            return self.__NSPOOLGM
        except AttributeError:
            return self.__get_NSPOOLG_('meteorological')

    @property
    def NOUTGC(self):
        try:
            return self.__NOUTGC
        except AttributeError:
            return self.__get_NOUTG_('concentration')

    @property
    def TOUTSGC(self):
        try:
            return self.__TOUTSGC
        except AttributeError:
            return self.__get_TOUTSG_('concentration')

    @property
    def TOUTFGC(self):
        try:
            return self.__TOUTFGC
        except AttributeError:
            return self.__get_TOUTFG_('concentration')

    @property
    def NSPOOLGC(self):
        try:
            return self.__NSPOOLGC
        except AttributeError:
            return self.__get_NSPOOLG_('concentration')

    @property
    def NFREQ(self):
        if self.runtype == 'coldstart':
            if np.any([_['spinup'] for _ in self._outputs]):
                if np.any([_['sampling_frequency'] for _ in self._outputs]):
                    if np.any([_['harmonic_analysis'] for _ in self._outputs]):
                        return len(self.TidalForcing.get_active_constituents())
        else:
            if np.any([_['sampling_frequency'] for _ in self._outputs]):
                if np.any([_['harmonic_analysis'] for _ in self._outputs]):
                    return len(self.TidalForcing.get_active_constituents())
        return 0

    @property
    def THAS(self):
        try:
            return self.__THAS
        except AttributeError:
            for _output in self._outputs:
                if _output['harmonic_analysis']:
                    if self.runtype == 'coldstart':
                        if _output['spinup']:
                            return self.STATIM
                    else:
                        dt = self.start_date - self.forcing_start_date
                        return (self.STATIM + dt.total_seconds())/(24.*60.*60.)
                else:
                    return self.STATIM

    @property
    def THAF(self):
        try:
            return self.__THAF
        except AttributeError:
            if self.NFREQ == 0:
                return 0
            dt = self.start_date - self.forcing_start_date
            return (self.STATIM + dt.total_seconds())/(24.*60.*60.)

    @property
    def NHAINC(self):
        try:
            return self.__NHAINC
        except AttributeError:
            NHAINC = float('inf')
            for _output in self._outputs:
                if _output['harmonic_analysis']:
                    if self.runtype == 'coldstart':
                        if _output['spinup']:
                            fs = _output['sampling_frequency']
                            NHAINC = np.min([NHAINC, int(fs.total_seconds())])
                    else:  # consider a "metonly" run?
                        fs = _output['sampling_frequency']
                        NHAINC = np.min([NHAINC, int(fs.total_seconds())])
            if NHAINC == float('inf'):
                NHAINC = 0
            return NHAINC

    @property
    def FMV(self):
        try:
            return self.__FMV
        except AttributeError:
            return 0

    @property
    def NHASE(self):
        try:
            return self.__NHASE
        except AttributeError:
            return self.__get_harmonic_analysis_state(
                self.get_elevation_output_stations())

    @property
    def NHASV(self):
        try:
            return self.__NHASV
        except AttributeError:
            return self.__get_harmonic_analysis_state(
                self.get_velocity_output_stations())

    @property
    def NHAGE(self):
        try:
            return self.__NHAGE
        except AttributeError:
            return self.__get_harmonic_analysis_state(
                self.get_elevation_global_output())

    @property
    def NHAGV(self):
        try:
            return self.__NHAGV
        except AttributeError:
            return self.__get_harmonic_analysis_state(
                self.get_velocity_global_output())

    @property
    def NHSTAR(self):
        try:
            return self.__NHSTAR
        except AttributeError:
            if self.runtype == 'coldstart':
                if self.netcdf is True:
                    return 5
                else:
                    return 3
            else:
                return 0

    @property
    def NHSINC(self):
        try:
            return self.__NHSINC
        except AttributeError:
            if self.NHSTAR == 0:
                return 0
            else:
                dt = self.start_date - self.forcing_start_date
                return int((dt.total_seconds())/self.DTDP)

    @property
    def ITITER(self):
        try:
            return self.__ITITER
        except AttributeError:
            return 1

    @property
    def ISLDIA(self):
        try:
            return self.__ISLDIA
        except AttributeError:
            return 0

    @property
    def CONVCR(self):
        try:
            return self.__CONVCR
        except AttributeError:
            # https://stackoverflow.com/questions/19141432/python-numpy-machine-epsilon
            # return 500*(7./3 - 4./3 - 1)
            return 1.E-6

    @property
    def ITMAX(self):
        try:
            return self.__ITMAX
        except AttributeError:
            return 25

    @property
    def NCPROJ(self):
        try:
            return self.__NCPROJ
        except AttributeError:
            return '{:<32}! NCPROJ'.format(' ')

    @property
    def NCINST(self):
        try:
            return self.__NCINST
        except AttributeError:
            return '{:<32}! NCINST'.format(' ')

    @property
    def NCSOUR(self):
        try:
            return self.__NCSOUR
        except AttributeError:
            return '{:<32}! NCSOUR'.format(' ')

    @property
    def NCHIST(self):
        try:
            return self.__NCHIST
        except AttributeError:
            return '{:<32}! NCHIST'.format(' ')

    @property
    def NCREF(self):
        try:
            return self.__NCREF
        except AttributeError:
            return '{:<32}! NCREF'.format(' ')

    @property
    def NCCOM(self):
        try:
            return self.__NCCOM
        except AttributeError:
            return '{:<32}! NCCOM'.format(' ')

    @property
    def NCHOST(self):
        try:
            return self.__NCHOST
        except AttributeError:
            return '{:<32}! NCHOST'.format(' ')

    @property
    def NCCONV(self):
        try:
            return self.__NCCONV
        except AttributeError:
            return '{:<32}! NCCONV'.format(' ')

    @property
    def NCCONT(self):
        try:
            return self.__NCCONT
        except AttributeError:
            return '{:<32}! NCCONT'.format(' ')

    @property
    def NCDATE(self):
        return self.forcing_start_date.strftime('%Y-%m-%d %H:%M')

    @property
    def TPXO(self):
        return self.__TPXO

    @property
    def FortranNamelists(self):
        return self.__FortranNamelists

    @runtype.setter
    def runtype(self, runtype):
        assert runtype in ['coldstart', 'hotstart']
        self.__runtype = runtype

    @RUNDES.setter
    def RUNDES(self, RUNDES):
        self.__RUNDES = str(RUNDES)

    @IHOT.setter
    def IHOT(self, IHOT):
        assert IHOT in [0]
        self.__IHOT = IHOT

    @RUNID.setter
    def RUNID(self, RUNID):
        self.__RUNID = str(RUNID)

    @NFOVER.setter
    def NFOVER(self, NFOVER):
        assert NFOVER in [0, 1]
        self.__NFOVER = NFOVER

    @WarnElev.setter
    def WarnElev(self, WarnElev):
        if WarnElev is not None:
            self.__WarnElev = float(WarnElev)
        else:
            self.__WarnElev = None

    @iWarnElevDump.setter
    def iWarnElevDump(self, iWarnElevDump):
        if iWarnElevDump is not None:
            iWarnElevDump = int(iWarnElevDump)
            if iWarnElevDump not in [0, 1]:
                raise TypeError('iWarnElevDump must be 0 or 1')
            self.__iWarnElevDump = int(iWarnElevDump)
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set iWarnElevDump if WarnElev is not '
                                   + 'None')

    @WarnElevDumpLimit.setter
    def WarnElevDumpLimit(self, WarnElevDumpLimit):
        if WarnElevDumpLimit is not None:
            assert isinstance(WarnElevDumpLimit, int)
            assert WarnElevDumpLimit > 0
            self.__WarnElevDumpLimit = WarnElevDumpLimit
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set WarnElevDumpLimit if WarnElev is '
                                   + 'not None')

    @ErrorElev.setter
    def ErrorElev(self, ErrorElev):
        if ErrorElev is not None:
            self.__ErrorElev = float(ErrorElev)
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set iWarnElevDump if WarnElev is not '
                                   + 'None')

    @NABOUT.setter
    def NABOUT(self, NABOUT):
        assert isinstance(NABOUT, int)
        assert NABOUT in [-1, 0, 1, 2, 3]
        self.__NABOUT = NABOUT

    @NSCREEN.setter
    def NSCREEN(self, NSCREEN):
        assert isinstance(NSCREEN, int)
        self.__NSCREEN = NSCREEN

    @IM.setter
    def IM(self, IM):
        assert isinstance(IM, int)
        if IM in [0, 111112]:
            self.__IM = IM
        elif IM in [1, 2, 611112]:
            raise NotImplementedError('3D runs not yet supported.')
        else:
            raise TypeError('IM must be 0, 1, 21, 111112 or 611112.')

    @IDEN.setter
    def IDEN(self, IDEN):
        if IDEN is not None:
            raise NotImplementedError('3D runs not yet supported.')

    @NOLIBF.setter
    def NOLIBF(self, NOLIBF):
        assert NOLIBF in [0, 1, 2]
        self.__NOLIBF = NOLIBF

    @NOLIFA.setter
    def NOLIFA(self, NOLIFA):
        NOLIFA = int(NOLIFA)
        assert NOLIFA in [0, 1, 2]
        self.__NOLIFA = NOLIFA

    @NOLICA.setter
    def NOLICA(self, NOLICA):
        NOLICA = int(NOLICA)
        assert NOLICA in [0, 1]
        self.__NOLICA = NOLICA

    @NOLICAT.setter
    def NOLICAT(self, NOLICAT):
        NOLICAT = int(NOLICAT)
        assert NOLICAT in [0, 1]
        self.__NOLICAT = NOLICAT

    @NCOR.setter
    def NCOR(self, NCOR):
        assert NCOR in [0, 1]
        self.__NCOR = NCOR

    @NTIP.setter
    def NTIP(self, NTIP):
        NTIP = int(NTIP)
        assert NTIP in [0, 1, 2]
        self.__NTIP = NTIP

    @G.setter
    def G(self, G):
        self.__G = float(G)

    @DTDP.setter
    def DTDP(self, DTDP):
        DTDP = float(DTDP)
        DTDP = np.abs(DTDP)
        assert DTDP != 0.
        self.__DTDP = DTDP

    @STATIM.setter
    def STATIM(self, STATIM):
        self.__STATIM = float(STATIM)

    @REFTIM.setter
    def REFTIM(self, REFTIM):
        self.__REFTIM = float(REFTIM)

    @DRAMP.setter
    def DRAMP(self, DRAMP):
        self.__DRAMP = float(DRAMP)

    @DRAMPExtFlux.setter
    def DRAMPExtFlux(self, DRAMPExtFlux):
        self.__DRAMPExtFlux = float(DRAMPExtFlux)

    @FluxSettlingTime.setter
    def FluxSettlingTime(self, FluxSettlingTime):
        self.__FluxSettlingTime = float(FluxSettlingTime)

    @DRAMPIntFlux.setter
    def DRAMPIntFlux(self, DRAMPIntFlux):
        self.__DRAMPIntFlux = float(DRAMPIntFlux)

    @DRAMPElev.setter
    def DRAMPElev(self, DRAMPElev):
        self.__DRAMPElev = float(DRAMPElev)

    @DRAMPTip.setter
    def DRAMPTip(self, DRAMPTip):
        self.__DRAMPTip = float(DRAMPTip)

    @DRAMPMete.setter
    def DRAMPMete(self, DRAMPMete):
        self.__DRAMPMete = float(DRAMPMete)

    @DRAMPWRad.setter
    def DRAMPWRad(self, DRAMPWRad):
        self.__DRAMPWRad = float(DRAMPWRad)

    @DUnRampMete.setter
    def DUnRampMete(self, DUnRampMete):
        if DUnRampMete is None:
            DUnRampMete = self.DRAMP
        self.__DUnRampMete = float(DUnRampMete)

    @A00.setter
    def A00(self, A00):
        self.__A00 = float(A00)

    @B00.setter
    def B00(self, B00):
        self.__B00 = float(B00)

    @C00.setter
    def C00(self, C00):
        self.__C00 = float(C00)

    @H0.setter
    def H0(self, H0):
        self.__H0 = float(H0)

    @NODEDRYMIN.setter
    def NODEDRYMIN(self, NODEDRYMIN):
        self.__NODEDRYMIN = int(NODEDRYMIN)

    @NODEWETRMP.setter
    def NODEWETRMP(self, NODEWETRMP):
        self.__NODEWETRMP = int(NODEWETRMP)

    @VELMIN.setter
    def VELMIN(self, VELMIN):
        self.__VELMIN = float(VELMIN)

    @SLAM0.setter
    def SLAM0(self, SLAM0):
        self.__SLAM0 = float(SLAM0)

    @SFEA0.setter
    def SFEA0(self, SFEA0):
        self.__SFEA0 = float(SFEA0)

    @FFACTOR.setter
    def FFACTOR(self, FFACTOR):
        self.__FFACTOR = float(FFACTOR)

    @ESLM.setter
    def ESLM(self, ESLM):
        self.__ESLM = float(ESLM)

    @HBREAK.setter
    def HBREAK(self, HBREAK):
        self.__HBREAK = float(HBREAK)

    @FTHETA.setter
    def FTHETA(self, FTHETA):
        self.__FTHETA = float(FTHETA)

    @FGAMMA.setter
    def FGAMMA(self, FGAMMA):
        self.__FGAMMA = float(FGAMMA)

    @ESLM.setter
    def ESLM(self, ESLM):
        self.__ESLM = float(ESLM)

    @ESLC.setter
    def ESLC(self, ESLC):
        self.__ESLC = ESLC

    @CORI.setter
    def CORI(self, CORI):
        if CORI is None:
            if self.NCOR == 0:
                raise Exception('Must pass CORI when NCOR=0')
            else:
                CORI = 0.
        else:
            CORI = float(CORI)
        self.__CORI = CORI

    # @NTIF.setter
    # def NTIF(self, NTIF):
    #     if self.TidalForcing is None:
    #         NTIF = []
    #     else:
    #         if NTIF == 'all':
    #             NTIF = []
    #             for constituent in self.TidalForcing.constituents:
    #                 if constituent in self.TidalForcing.major8:
    #                     NTIF.append(constituent)
    #         else:
    #             NTIF = list(NTIF)
    #     self.__NTIF = NTIF

    @ANGINN.setter
    def ANGINN(self, ANGINN):
        self.__ANGINN = float(ANGINN)

    @THAS.setter
    def THAS(self, THAS):
        THAS = float(THAS)
        assert THAS >= 0.
        self.__THAS = THAS

    @THAF.setter
    def THAF(self, THAF):
        THAF = float(THAF)
        assert THAF >= 0.
        self.__THAF = THAF

    @NHAINC.setter
    def NHAINC(self, NHAINC):
        NHAINC = int(NHAINC)
        assert NHAINC >= 0
        self.__NHAINC = NHAINC

    @FMV.setter
    def FMV(self, FMV):
        FMV = float(FMV)
        assert FMV >= 0. and FMV <= 1.
        self.__FMV = FMV

    @NHSTAR.setter
    def NHSTAR(self, NHSTAR):
        assert NHSTAR in [0, 1, 2, 3, 5]
        self.__NHSINC = NHSTAR

    @NHSINC.setter
    def NHSINC(self, NHSINC):
        self.__NHSINC = int(NHSINC)

    @ITITER.setter
    def ITITER(self, ITITER):
        ITITER = int(ITITER)
        assert ITITER in [1, -1]
        self.__ITITER = ITITER

    @ISLDIA.setter
    def ISLDIA(self, ISLDIA):
        ISLDIA = int(ISLDIA)
        assert ISLDIA in [0, 1, 2, 3, 4, 5]
        self.__ISLDIA = ISLDIA

    @CONVCR.setter
    def CONVCR(self, CONVCR):
        self.__CONVCR = float(CONVCR)

    @ITMAX.setter
    def ITMAX(self, ITMAX):
        self.__ITMAX = int(ITMAX)

    @NCPROJ.setter
    def NCPROJ(self, NCPROJ):
        self.__NCPROJ = str(NCPROJ)

    @NCINST.setter
    def NCINST(self, NCINST):
        self.__NCINST = str(NCINST)

    @NCSOUR.setter
    def NCSOUR(self, NCSOUR):
        self.__NCSOUR = str(NCSOUR)

    @NCHIST.setter
    def NCHIST(self, NCHIST):
        self.__NCHIST = str(NCHIST)

    @NCREF.setter
    def NCREF(self, NCREF):
        self.__NCREF = str(NCREF)

    @NCCOM.setter
    def NCCOM(self, NCCOM):
        self.__NCCOM = str(NCCOM)

    @NCHOST.setter
    def NCHOST(self, NCHOST):
        self.__NCHOST = str(NCHOST)

    @NCCONV.setter
    def NCCONV(self, NCCONV):
        self.__NCCONV = str(NCCONV)

    @NCCONT.setter
    def NCCONT(self, NCCONT):
        self.__NCCONT = str(NCCONT)
