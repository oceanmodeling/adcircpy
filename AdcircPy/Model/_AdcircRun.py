# global imports
import numpy as np
from pathlib import Path

# local imports
from AdcircPy import Model
from AdcircPy.Model._Fort15 import _Fort15


class _AdcircRun(_Fort15):

    def __init__(self,
                 AdcircMesh,
                 TidalForcing=None,
                 WindForcing=None,
                 netcdf=True,
                 ElevationGlobalOutput=None,
                 VelocityGlobalOutput=None,
                 MeteorologicalGlobalOutput=None,
                 ElevationStationsOutput=None,
                 VelocityStationsOutput=None,
                 MeteorologicalStationsOutput=None, **fort15):
        self.AdcircMesh = AdcircMesh
        self.TidalForcing = TidalForcing
        self.WindForcing = WindForcing
        if not self.TidalForcing and not self.WindForcing:
            raise RuntimeError('Must pass at least one of TidalForcing or '
                               + 'WindForcing objects.')
        self.__set_start_date()
        self.__set_end_date()
        self.__set_forcing_start_date()
        self.netcdf = netcdf
        self.ElevationGlobalOutput = ElevationGlobalOutput
        self.VelocityGlobalOutput = VelocityGlobalOutput
        self.MeteorologicalGlobalOutput = MeteorologicalGlobalOutput
        self.ElevationStationsOutput = ElevationStationsOutput
        self.VelocityStationsOutput = VelocityStationsOutput
        self.MeteorologicalStationsOutput = MeteorologicalStationsOutput
        super(_AdcircRun, self).__init__(**fort15)

    def dump(self, output_dir=None, coldstart=True, hotstart=True,
             filename_coldstart='fort.15.coldstart',
             filename_hotstart='fort.15.hotstart',
             fort13_filename='fort.13'):
        if output_dir is not None:
            output_dir = Path(output_dir)
        self.AdcircMesh.NodalAttributes.dump(
                                        output_dir, filename=fort13_filename)
        if coldstart:
            fort15 = self.get_fort15('coldstart')
            if output_dir is None:
                print(fort15)
            else:
                _ = Path(str(output_dir)+'/'+filename_coldstart)
                with open(str(_), 'w') as f:
                    f.write(fort15)
        if hotstart:
            fort15 = self.get_fort15('hotstart')
            if output_dir is None:
                print(fort15)
            else:
                _ = Path(str(output_dir)+'/'+filename_hotstart)
                with open(_, 'w') as f:
                    f.write(fort15)

    def get_fort15(self, runtype):
        assert runtype in ['coldstart', 'hotstart', 'metonly']
        f = '{:<32}! RUNDES\n'.format(self.RUNDES)
        f += '{:<32}! RUNID\n'.format(self.RUNID, '')
        f += '{:<32d}! NFOVER\n'.format(self.NFOVER)
        f += '{:<32d}! NABOUT\n'.format(self.NABOUT)
        f += '{:<32d}! NSCREEN\n'.format(self.NSCREEN)
        if runtype == 'coldstart':
            f += '{:<32d}! IHOT\n'.format(0)
        elif runtype == 'hotstart':
            dt = self.start_date - self.forcing_start_date
            NHSINC = int(dt.total_seconds()/self.DTDP)/self.NHSINC
            if NHSINC == 1:
                f += '{:<32d}! IHOT\n'.format(567)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError
        f += '{:<32d}! ICS\n'.format(self.ICS)
        f += '{:<32d}! IM\n'.format(self.IM)
        if self.IM == 21:
            f += '{}\n'.format(self.IDEN)
        f += '{:<32d}! NOLIBF\n'.format(self.NOLIBF)
        f += '{:<32d}! NOLIFA\n'.format(self.NOLIFA)
        f += '{:<32d}! NOLICA\n'.format(self.NOLICA)
        f += '{:<32d}! NOLICAT\n'.format(self.NOLICAT)
        if runtype == 'coldstart':
            NWP = len(self.AdcircMesh.NodalAttributes.coldstart)
            f += '{:<32d}! NWP\n'.format(NWP)
            for attribute in self.AdcircMesh.NodalAttributes.coldstart:
                f += '{:<32}\n'.format(attribute)
        elif runtype == 'hotstart':
            NWP = len(self.AdcircMesh.NodalAttributes.hotstart)
            f += '{:<32d}! NWP\n'.format(NWP)
            for attribute in self.AdcircMesh.NodalAttributes.hotstart:
                f += '{:<32}\n'.format(attribute)
        else:
            raise NotImplementedError('Need to sort this part for metonly')
        f += '{:<32d}! NCOR\n'.format(self.NCOR)
        f += '{:<32d}! NTIP\n'.format(self.NTIP)
        if runtype == 'coldstart':
            NWS = 0
            f += '{:<32d}! NWS\n'.format(NWS)
            f += '{:<32d}! NRAMP\n'.format(1)
        elif runtype == 'hotstart':
            if self.WindForcing:
                NWS = self.WindForcing.NWS
                f += '{:<32d}! NWS\n'.format(NWS)
                f += '{:<32d}! NRAMP\n'.format(8)  # NRAMP
            else:
                # tidal only hostart
                NWS = 0
                f += '{:<32d}! NWS\n'.format(NWS)
                f += '{:<32d}! NRAMP\n'.format(0)  # NRAMP
        f += '{:<32.2f}! G\n'.format(self.G)
        TAU0 = self.AdcircMesh.NodalAttributes.TAU0
        if TAU0 is None:
            f += '{:<32.4f}! TAU0\n'.format(self.FFACTOR)
        elif TAU0 in [3, -3]:
            TAU0 = -3
            f += '{:<32d}! TAU0\n'.format(TAU0)
        elif TAU0 == -5:
            PW = self.NodalAttributes[
                'primitive_weighting_in_continuity_equation']
            f += '{:<32d}! TAU0\n'.format(TAU0)
            f += '{} {}\n'.format(PW.Tau0FullDomainMin, PW.Tau0FullDomainMax)
        else:
            f += '{:<32d}! TAU0\n'.format(TAU0)
        f += '{:<32.2f}! DTDP\n'.format(self.DTDP)
        f += '{:<32.2f}! STATIM\n'.format(self.STATIM)
        f += '{:<32.2f}! REFTIM\n'.format(self.REFTIM)
        if runtype in ['hotstart', 'metonly']:
            if self.WindForcing:
                f += '{:<32d} ! WTIMINC\n'.format(self.WindForcing.WTIMINC)
        if runtype == 'coldstart':
            RNDAY = self.start_date - self.forcing_start_date
        else:
            RNDAY = self.end_date - self.forcing_start_date
        f += '{:<32.2f}! RNDAY\n'.format(RNDAY.total_seconds()/(60.*60.*24.))
        if runtype == 'coldstart':
            DRAMP = self.DRAMP
        elif runtype == 'hotstart':
            if self.TidalForcing and not self.WindForcing:
                DRAMP = self.DRAMP
            else:
                raise NotImplementedError('Need to add metramp info')
        f += '{:<32.2f}! DRAMP\n'.format(DRAMP)
        f += '{:<4.2f} {:<4.2f} {:<4.2f}{:<18}! A00 B00 C00\n'.format(
                self.A00, self.B00, self.C00, '')
        f += '{:<5.3f} {:<2d} {:<2d} {:<5.3f}{:15}'.format(
            self.H0, self.NODEDRYMIN, self.NODEWETRMP, self.VELMIN, '')
        f += '! H0 NODEDRYMIN NODEWETRMP VELMIN\n'
        f += '{:<4.1f} {:<4.1f}{:22}'.format(self.SLAM0, self.SFEA0, '')
        f += '! SLAM0 SFEA0\n'
        if self.NOLIBF in [0, 1]:
            f += '{:<32.4f}! FFACTOR\n'.format(self.FFACTOR)
        elif self.NOLIBF == 2:
            raise NotImplementedError('Need to setup for NOLIBF==2')
        if self.IM in [0, 1, 2]:
            f += '{:<32.4f}! ESLM (eddy viscosity)\n'.format(self.ESLM)
        elif self.IM == 10:
            raise NotImplementedError('')
        f += '{:<32.4f}! CORI\n'.format(self.CORI)
        f += '{:<32d}! NTIF\n'.format(len(self.NTIF))
        for constituent in self.NTIF:
            forcing = self.TidalForcing(constituent)
            f += '{:<32}\n'.format(constituent)
            f += '{:<.4E} '.format(forcing[0])
            f += '{:<.6E} '.format(forcing[1])
            f += '{:<4.2f} '.format(forcing[2])
            f += '{:<.5E} '.format(forcing[3])
            f += '{:<.5E} '.format(forcing[4])
            f += '\n'
        NBFR = len(self.TidalForcing.constituents)
        f += '{:<32d}! NBFR\n'.format(NBFR)
        for constituent in self.TidalForcing.constituents:
            if constituent in self.NTIF:
                forcing = self.TidalForcing(constituent)
                f += '{:<32}\n'.format(constituent)
                f += '{:<.6E} '.format(forcing[1])
                f += '{:<.5E} '.format(forcing[3])
                f += '{:<.5E} '.format(forcing[4])
                f += '\n'
        for constituent in self.TidalForcing.constituents:
            if constituent not in self.NTIF:
                forcing = self.TidalForcing(constituent)
                f += '{:<32}\n'.format(constituent)
                f += '{:<.6E} '.format(forcing[1])
                f += '{:<.5E} '.format(forcing[3])
                f += '{:<.5E} '.format(forcing[4])
                f += '\n'
        # NOTE:  This part is written as one-constituent then all boundaries
        # as opposed to one-boundary then all constituents for that boundary.
        # Not exactly sure how ADCIRC handles multiple open boundaries.
        for constituent in self.TidalForcing.constituents:
            if constituent in self.NTIF:
                for i, boundary in enumerate(self.TidalForcing.TPXO_interp):
                    if i == 0:
                        f += '{:<32}\n'.format(constituent)
                        tpxo = boundary[constituent]
                        for ha, hp in tpxo:
                            f += '{:<.6E} '.format(ha)
                            f += '{:>8.3f} '.format(hp)
                            f += '\n'
        for constituent in self.TidalForcing.constituents:
            if constituent not in self.NTIF:
                for i, boundary in enumerate(self.TidalForcing.TPXO_interp):
                    if i == 0:
                        f += '{:<32}\n'.format(constituent)
                        tpxo = boundary[constituent]
                        for ha, hp in tpxo:
                            f += '{:<.6E} '.format(ha)
                            f += '{:>8.3f} '.format(hp)
                            f += '\n'
        f += "{:<32.2f}! ANGINN\n".format(self.ANGINN)
        # other boundary type forcings go here as NFFR (e.g. river boundaries)

        params = self.ElevationStationsOutput(
                self.AdcircMesh, runtype, self.forcing_start_date,
                self.start_date, self.end_date, self.DTDP)
        (NOUTE, NOUTSTE, TOUTFE, NSPOOLE, stations) = params
        f += "{:<2d} ".format(NOUTE)
        f += "{:<4.2f} ".format(NOUTSTE)
        f += "{:<4.2f} ".format(TOUTFE)
        f += "{:<10d}\n".format(NSPOOLE)
        f += "{:<32d}! NSTAE\n".format(len(stations.keys()))
        for station_id, (x, y) in stations.items():
            f += "{:<10.4f} ".format(x)
            f += "{:<10.4f} ".format(y)
            f += "{:20}".format('')
            f += "! {}\n".format(station_id)
        params = self.VelocityStationsOutput(
                self.AdcircMesh, runtype, self.forcing_start_date,
                self.start_date, self.end_date, self.DTDP)
        (NOUTV, NOUTSTV, TOUTFV, NSPOOLV, stations) = params
        f += "{:<2d} ".format(NOUTV)
        f += "{:<4.2f} ".format(NOUTSTV)
        f += "{:<4.2f} ".format(TOUTFV)
        f += "{:<10d}\n".format(NSPOOLV)
        f += "{:<32d}! NSTAV\n".format(len(stations.keys()))
        for station_id, (x, y) in stations.items():
            f += "{:<10.4f} ".format(x)
            f += "{:<10.4f} ".format(y)
            f += "! {}\n".format(station_id)
        if self.IM == 10:
            raise NotImplementedError('Concentration outputs go here.')
        if runtype == 'hotstart':
            if self.WindForcing:
                params = self.MeteorologicalStationsOutput(
                        self.AdcircMesh, runtype, self.forcing_start_date,
                        self.start_date, self.end_date, self.DTDP)
                (NOUTW, NOUTSTW, TOUTFW, NSPOOLW, stations) = params
                f += "{:<2d} ".format(NOUTW)
                f += "{:<4.2f} ".format(NOUTSTW)
                f += "{:<4.2f} ".format(TOUTFW)
                f += "{:<10d}\n".format(NSPOOLW)
                f += "{:<32d}! NSTAM\n".format(len(stations.keys()))
                for station_id, (x, y) in stations.items():
                    f += "{:<10.4f} ".format(x)
                    f += "{:<10.4f} ".format(y)
                    f += "! {}\n".format(station_id)
        params = self.ElevationGlobalOutput(
                    runtype, self.forcing_start_date, self.start_date,
                    self.end_date, self.DTDP)
        (NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE) = params
        f += "{:<2d} ".format(NOUTGE)
        f += "{:<4.2f} ".format(TOUTSGE)
        f += "{:<4.2f} ".format(TOUTFGE)
        f += "{:<24d}".format(NSPOOLGE)
        f += "! NOUTGE TOUTSGE TOUTFGE NSPOOLGE\n"
        params = self.VelocityGlobalOutput(
                    runtype, self.forcing_start_date, self.start_date,
                    self.end_date, self.DTDP)
        (NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV) = params
        f += "{:<2d} ".format(NOUTGV)
        f += "{:<4.2f} ".format(TOUTSGV)
        f += "{:<4.2f} ".format(TOUTFGV)
        f += "{:<24d}".format(NSPOOLGV)
        f += "! NOUTGV TOUTSGV TOUTFGV NSPOOLGV\n"
        if self.IM == 10:
            raise NotImplementedError
        if runtype == 'hotstart':
            if NWS != 0:
                params = self.MeteorologicalGlobalOutput(
                        runtype, self.forcing_start_date, self.start_date,
                        self.end_date, self.DTDP)
                (NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW) = params
                f += "{:<2d} ".format(NOUTGW)
                f += "{:<4.2f} ".format(TOUTSGW)
                f += "{:<4.2f} ".format(TOUTFGW)
                f += "{:<24d}".format(NSPOOLGW)
                f += "! NOUTGW TOUTSGW TOUTFGW NSPOOLGW\n"
        harmonic_analysis = [self.ElevationStationsOutput,
                             self.VelocityStationsOutput,
                             self.ElevationGlobalOutput,
                             self.VelocityGlobalOutput]
        harmonic_analysis = list(filter(lambda x: x.harmonic_analysis is True,
                                        harmonic_analysis))
        if runtype == 'coldstart':
            spinup = list(filter(lambda x: x.spinup is True,
                                 harmonic_analysis))
            if len(spinup) > 0:
                if len(harmonic_analysis) > 0:
                    NFREQ = len(self.TidalForcing.constituents)
                    THAS = 0.
                    dt = self.start_date - self.forcing_start_date
                    THAF = dt.total_seconds()/(24.*60.*60.)
                    NHAINC = [x.sampling_frequency.total_seconds()
                              for x in spinup]
                    NHAINC = int(np.min(NHAINC))
                else:
                    NFREQ = 0
                    THAS = 0
                    THAF = 0
                    NHAINC = 0
            else:
                NFREQ = 0
                THAS = 0
                THAF = 0
                NHAINC = 0
        elif runtype == 'hotstart':
            if len(harmonic_analysis) > 0:
                NFREQ = len(self.TidalForcing.constituents)
                dt = self.start_date - self.forcing_start_date
                THAS = dt.total_seconds()/(24.*60.*60.)
                dt = self.end_date - self.forcing_start_date
                THAF = dt.total_seconds()/(24.*60.*60.)
                NHAINC = [x.sampling_frequency.total_seconds()
                          for x in harmonic_analysis]
                NHAINC = int(np.min(NHAINC))
            else:
                NFREQ = 0
                THAS = 0
                THAF = 0
                NHAINC = 0
        elif runtype == 'metonly':
            raise NotImplementedError(
                'Can we do harmonic analysis on velocity outputs for metonly?')
            NFREQ = 0
            THAS = 0
            THAF = 0
            NHAINC = 0
        f += "{:<32d}! NFREQ\n".format(NFREQ)
        if NFREQ > 0:
            for constituent in self.TidalForcing.constituents:
                if constituent in self.NTIF:
                    forcing = self.TidalForcing(constituent)
                    f += '{:<32}\n'.format(constituent)
                    f += '{:<.6E} '.format(forcing[1])
                    f += '{:<.5E} '.format(forcing[3])
                    f += '{:<.5E} '.format(forcing[4])
                    f += '\n'
            for constituent in self.TidalForcing.constituents:
                if constituent not in self.NTIF:
                    forcing = self.TidalForcing(constituent)
                    f += '{:<32}\n'.format(constituent)
                    f += '{:<.6E} '.format(forcing[1])
                    f += '{:<.5E} '.format(forcing[3])
                    f += '{:<.5E} '.format(forcing[4])
                    f += '\n'
        if self.NHAINC is not None:
            NHAINC = self.NHAINC
        f += "{:<6.2f} ".format(THAS)
        f += "{:<6.2f} ".format(THAF)
        f += "{:<6d} ".format(NHAINC)
        f += "{:<3.1f}".format(self.FMV)
        f += 7 * " "
        f += "! THAS THAF NHAINC FMV\n"
        NHASE = self.ElevationStationsOutput._get_NHA(runtype)
        NHASV = self.VelocityStationsOutput._get_NHA(runtype)
        NHAGE = self.ElevationGlobalOutput._get_NHA(runtype)
        NHAGV = self.VelocityGlobalOutput._get_NHA(runtype)
        f += "{:<2d} ".format(NHASE)
        f += "{:<2d} ".format(NHASV)
        f += "{:<2d} ".format(NHAGE)
        f += "{:<2d}".format(NHAGV)
        f += 7 * " "
        f += "! NHASE NHASV NHAGE NHAGV\n"
        if runtype == 'coldstart':
            if self.netcdf:
                NHSTAR = 5
            else:
                NHSTAR = 1
        elif runtype == 'hotstart':
            NHSTAR = 0
        f += "{:<4d}".format(NHSTAR)
        f += "{:<20d}\n".format(self.NHSINC)
        f += "{:<4d}".format(self.ITITER)
        f += "{:<4d}".format(self.ISLDIA)
        f += "{:<6.2f}".format(self.CONVCR)
        f += "{:<4d}\n".format(self.ITMAX)
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
        f += "{} UTC".format(self.NCDATE)
        f += 6*" "
        f += "! model start_date {}\n".format(
                            self.start_date.strftime('%Y-%m-%d %H:%M UTC'))
        return f

    def __set_start_date(self):
        if self.TidalForcing and self.WindForcing:
            assert self.TidalForcing.start_date == self.WindForcing.start_date
            start_date = self.TidalForcing.start_date
        elif self.WindForcing and not self.TidalForcing:
            start_date = self.WindForcing.start_date
        else:
            start_date = self.TidalForcing.start_date
        self.__start_date = start_date

    def __set_end_date(self):
        if self.TidalForcing and self.WindForcing:
            assert self.TidalForcing.end_date == self.WindForcing.end_date
            self.__end_date = self.__TidalForcing.end_date
            end_date = self.TidalForcing.end_date
        elif self.WindForcing and not self.TidalForcing:
            self.__end_date = self.__WindForcing.end_date
        else:
            end_date = self.TidalForcing.end_date
        self.__end_date = end_date

    def __set_forcing_start_date(self):
        if self.TidalForcing and self.WindForcing:
            forcing_start_date = self.TidalForcing.forcing_start_date
        elif self.WindForcing and not self.TidalForcing:
            forcing_start_date = self.WindForcing.start_date
        else:
            forcing_start_date = self.TidalForcing.forcing_start_date
        self.__forcing_start_date = forcing_start_date

    @property
    def AdcircMesh(self):
        return self.__AdcircMesh

    @property
    def TidalForcing(self):
        return self.__TidalForcing

    @property
    def WindForcing(self):
        return self.__WindForcing

    @property
    def start_date(self):
        return self.__start_date

    @property
    def end_date(self):
        return self.__end_date

    @property
    def forcing_start_date(self):
        return self.__forcing_start_date

    @property
    def netcdf(self):
        return self.__netcdf

    @property
    def ElevationGlobalOutput(self):
        return self.__ElevationGlobalOutput

    @property
    def VelocityGlobalOutput(self):
        return self.__VelocityGlobalOutput

    @property
    def MeteorologicalGlobalOutput(self):
        return self.__MeteorologicalGlobalOutput

    @property
    def ElevationStationsOutput(self):
        return self.__ElevationStationsOutput

    @property
    def VelocityStationsOutput(self):
        return self.__VelocityStationsOutput

    @property
    def MeteorologicalStationsOutput(self):
        return self.__MeteorologicalStationsOutput

    @AdcircMesh.setter
    def AdcircMesh(self, AdcircMesh):
        assert isinstance(AdcircMesh, Model.AdcircMesh)
        if AdcircMesh.SpatialReference is None:
            raise RuntimeError('Must set mesh SpatialReference before '
                               + 'instantiating a run.')
        self.__AdcircMesh = AdcircMesh

    @TidalForcing.setter
    def TidalForcing(self, TidalForcing):
        if TidalForcing is not None:
            assert isinstance(TidalForcing, Model._TidalForcing)
        else:
            TidalForcing = False
        self.__TidalForcing = TidalForcing
        if self.TidalForcing is not None:
            self.TidalForcing._set_TPXO_interp(self.AdcircMesh)

    @WindForcing.setter
    def WindForcing(self, WindForcing):
        if WindForcing is not None:
            assert isinstance(WindForcing, Model._WindForcing)
            raise NotImplementedError
        else:
            WindForcing = False
        self.__WindForcing = WindForcing

    @netcdf.setter
    def netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf

    @ElevationGlobalOutput.setter
    def ElevationGlobalOutput(self, ElevationGlobalOutput):
        if ElevationGlobalOutput is not None:
            assert isinstance(ElevationGlobalOutput,
                              Model.ElevationGlobalOutput)
        else:
            ElevationGlobalOutput = Model.ElevationGlobalOutput()
        self.__ElevationGlobalOutput = ElevationGlobalOutput

    @VelocityGlobalOutput.setter
    def VelocityGlobalOutput(self, VelocityGlobalOutput):
        if VelocityGlobalOutput is not None:
            assert isinstance(VelocityGlobalOutput, Model.VelocityGlobalOutput)
        else:
            VelocityGlobalOutput = Model.VelocityGlobalOutput()
        self.__VelocityGlobalOutput = VelocityGlobalOutput

    @MeteorologicalGlobalOutput.setter
    def MeteorologicalGlobalOutput(self, MeteorologicalGlobalOutput):
        if MeteorologicalGlobalOutput is not None:
            assert isinstance(MeteorologicalGlobalOutput,
                              Model.MeteorologicalGlobalOutput)
        else:
            MeteorologicalGlobalOutput = Model.MeteorologicalGlobalOutput()
        self.__MeteorologicalGlobalOutput = MeteorologicalGlobalOutput

    @ElevationStationsOutput.setter
    def ElevationStationsOutput(self, ElevationStationsOutput):
        if ElevationStationsOutput is not None:
            assert isinstance(ElevationStationsOutput,
                              Model.ElevationStationsOutput)
        else:
            ElevationStationsOutput = Model.ElevationStationsOutput()
        self.__ElevationStationsOutput = ElevationStationsOutput

    @VelocityStationsOutput.setter
    def VelocityStationsOutput(self, VelocityStationsOutput):
        if VelocityStationsOutput is not None:
            assert isinstance(VelocityStationsOutput,
                              Model.VelocityStationsOutput)
        else:
            VelocityStationsOutput = Model.VelocityStationsOutput()
        self.__VelocityStationsOutput = VelocityStationsOutput

    @MeteorologicalStationsOutput.setter
    def MeteorologicalStationsOutput(self, MeteorologicalStationsOutput):
        if MeteorologicalStationsOutput is not None:
            assert isinstance(MeteorologicalStationsOutput,
                              Model.MeteorologicalStationsOutput)
        else:
            MeteorologicalStationsOutput = Model.MeteorologicalStationsOutput()
        self.__MeteorologicalStationsOutput = MeteorologicalStationsOutput
