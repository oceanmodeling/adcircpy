# global imports
import argparse
from datetime import datetime, timedelta
import uuid


# local imports
import AdcircPy
from AdcircPy.Model import ElevationGlobalOutput as EGO
from AdcircPy.Model import VelocityGlobalOutput as VGO
from AdcircPy.Model import MeteorologicalGlobalOutput as MGO
from AdcircPy.Model import ElevationStationsOutput as ESO
from AdcircPy.Model import VelocityStationsOutput as VSO
from AdcircPy.Model import MeteorologicalStationsOutput as MSO


class GenerateBestTrackRun(GenerateTidalRun):

    def __init__(self):
        self.__set_args()
        self.__set_AdcircMesh()
        # self.__set_start_datetime()
        # self.__set_end_datetime()
        # self.__set_fort13()
        # self.__set_coldstart_attributes()
        # self.__set_hotstart_attributes()
        # self.__set_TAU0()
        # self.__set_FFACTOR()
        # self.__set_EGO()
        # self.__set_VGO()
        # self.__set_MGO()
        # self.__set_stations()
        # self.__set_ESO()
        # self.__set_VSO()
        # self.__set_MSO()
        # self.__set_constituents()
        # self.__set_TidalRun()
        # self.__write_to_disk()

    def __set_args(self):
        parser = argparse.ArgumentParser(
            description="Program to generate ADCIRC tidal only run.")
        parser.add_argument(
            "fort14", type=str,
            help="Full path to mesh file. Any filename can be used. (Filename "
            + "doesn't have to be fort.14 nor end with .14 extension.)")
        parser.add_argument("ATCF_name")
        parser.add_argument(
            "--no-tidal-forcing", type=bool, default=False,
            help="Full path to mesh file. Any filename can be used. (Filename "
            + "doesn't have to be fort.14 nor end with .14 extension.)")

    def __set_AdcircMesh(self):
        self.AdcircMesh = AdcircPy.read_mesh(self.args.fort14,
                                             self.args.mesh_epsg,
                                             self.args.mesh_vertical_datum)

    # def __set_start_datetime(self):
    #     self.start_datetime = datetime.strptime(
    #                         self.args.start_datetime, "%Y-%m-%dT%H:%M")

    # def __set_end_datetime(self):
    #     self.end_datetime = self.start_datetime + timedelta(
    #                             days=self.args.hotstart_duration)

    # def __set_fort13(self):
    #     if self.args.fort13 is not None:
    #         self.AdcircMesh.set_fort13(self.args.fort13)

    # def __set_TAU0(self):
    #     if self.args.TAU0 != 'None':
    #         if 'primitive_weighting_in_continuity_equation' \
    #           not in self.AdcircMesh.NodalAttributes._storage.keys():
    #             self.AdcircMesh.set_TAU0(3)
    #         else:
    #             self.AdcircMesh.set_nodal_attribute_state(
    #                 'primitive_weighting_in_continuity_equation',
    #                 coldstart=True, hotstart=True)
    #     if self.AdcircMesh.NodalAttributes.TAU0 in [3, -3]:
    #         self.AdcircMesh.set_nodal_attribute_state(
    #             'primitive_weighting_in_continuity_equation',
    #             coldstart=True, hotstart=True)

    # def __set_coldstart_attributes(self):
    #     if len(self.args.coldstart_attributes) == 0:
    #         self.args.coldstart_attributes \
    #             = [key for key, _ in self.AdcircMesh.NodalAttributes]
    #     else:
    #         for attribute in self.args.coldstart_attributes:
    #             self.AdcircMesh.set_nodal_attribute_state(attribute,
    #                                                       coldstart=True,
    #                                                       hotstart=False)

    # def __set_hotstart_attributes(self):
    #     if len(self.args.hotstart_attributes) == 0:
    #         self.args.coldstart_attributes \
    #             = [key for key, _ in self.AdcircMesh.NodalAttributes]
    #     for attribute in self.args.hotstart_attributes:
    #         self.AdcircMesh.set_nodal_attribute(attribute,  coldstart=True,
    #                                             hotstart=False)

    # def __set_FFACTOR(self):
    #     if self.args.TAU0 == 'None':
    #         if self.args.FFACTOR is None:
    #             raise Exception('Must pass --FFACTOR if --TAU0=None')

    # def __set_EGO(self):
    #     if self.args.ego_sampling_frequency > 0:
    #         self.EGO = EGO(sampling_frequency=timedelta(
    #                                 minutes=self.args.ego_sampling_frequency),
    #                        netcdf=self.args.netcdf,
    #                        spinup=self.args.ego_coldstart,
    #                        harmonic_analysis=self.args.ego_harmonic_analysis)
    #     else:
    #         self.EGO = None

    # def __set_VGO(self):
    #     if self.args.vgo_sampling_frequency > 0:
    #         self.VGO = VGO(sampling_frequency=timedelta(
    #                             minutes=self.args.vgo_sampling_frequency),
    #                        netcdf=self.args.netcdf,
    #                        spinup=self.args.vgo_coldstart,
    #                        harmonic_analysis=self.args.vgo_harmonic_analysis)
    #     else:
    #         self.VGO = None

    # def __set_MGO(self):
    #     if self.args.mgo_sampling_frequency > 0:
    #         self.MGO = MGO(sampling_frequency=timedelta(
    #                                 minutes=self.args.mgo_sampling_frequency),
    #                        netcdf=self.args.netcdf)
    #     else:
    #         self.MGO = None

    # def __set_stations(self):
    #     self.stations = dict()
    #     for station in self.args.station:
    #         _station = station.split(',')
    #         if len(_station) < 2 or len(_station) > 3:
    #             raise RuntimeError('Cannot parse station '
    #                                + '"{}". '.format(station)
    #                                + 'Argument must be of format "x,y" or '
    #                                + ' "x,y,id"')
    #         if len(_station) == 2:
    #             key = str(uuid.uuid4())[:8]
    #         elif len(_station) == 3:
    #             key = _station[2]
    #         if key in self.stations.keys():
    #             raise Exception('Non-unique key: {}'.format(key))
    #         self.stations[key] = dict()
    #         self.stations[key]['x'] = float(_station[0])
    #         self.stations[key]['y'] = float(_station[1])

    # def __set_ESO(self):
    #     sampling_frequency = timedelta(
    #                                 minutes=self.args.eso_sampling_frequency)
    #     df = sampling_frequency.total_seconds()
    #     if df > 0:
    #         netcdf = self.args.netcdf
    #         spinup = self.args.eso_coldstart
    #         harmonic_analysis = self.args.eso_harmonic_analysis
    #         self.ESO = ESO(sampling_frequency, netcdf, spinup,
    #                        harmonic_analysis)
    #         if self.args.stations_file is not None:
    #             self.ESO.add_stations_from_fort15(self.args.stations_file)

    #     else:
    #         self.ESO = None

    # def __set_VSO(self):
    #     sampling_frequency = timedelta(
    #                                 minutes=self.args.vso_sampling_frequency)
    #     df = sampling_frequency.total_seconds()
    #     if df > 0:
    #         netcdf = self.args.netcdf
    #         spinup = self.args.vso_coldstart
    #         harmonic_analysis = self.args.vso_harmonic_analysis
    #         self.VSO = VSO(sampling_frequency, netcdf, spinup,
    #                        harmonic_analysis)
    #         if self.args.stations_file is not None:
    #             self.VSO.add_stations_from_fort15(self.args.stations_file)

    #     else:
    #         self.VSO = None

    # def __set_MSO(self):
    #     sampling_frequency = timedelta(
    #                                 minutes=self.args.mso_sampling_frequency)
    #     df = sampling_frequency.total_seconds()
    #     if df > 0:
    #         netcdf = self.args.netcdf
    #         self.MSO = MSO(sampling_frequency, netcdf)
    #         if self.args.stations_file is not None:
    #             self.MSO.add_stations_from_fort15(self.args.stations_file)

    #     else:
    #         self.MSO = None

    # def __set_constituents(self):
    #     self.constituents = set()
    #     constituents = set(self.args.constituents)
    #     if len(constituents) == 0 or 'all' in constituents:
    #         self.constituents = 'all'
    #     elif 'major' in constituents:
    #         self.constituents = 'major'
    #     else:
    #         self.constituents = constituents

    # def __set_TidalRun(self):
    #     self.TidalRun = self.AdcircMesh.TidalRun(
    #                     self.start_datetime,
    #                     self.end_datetime,
    #                     spinup_days=self.args.coldstart_duration,
    #                     constituents=self.constituents,
    #                     ElevationGlobalOutput=self.EGO,
    #                     VelocityGlobalOutput=self.VGO,
    #                     MeteorologicalGlobalOutput=self.MGO,
    #                     ElevationStationsOutput=self.ESO,
    #                     VelocityStationsOutput=self.VSO,
    #                     MeteorologicalStationsOutput=self.MSO,
    #                     FFACTOR=self.args.FFACTOR,
    #                     **self.fort15_kwargs)

    # def __write_to_disk(self):
    #     self.TidalRun.dump(self.args.output_directory)


def main():
    GenerateBestTrackRun()























# from AdcircPy.Model._AdcircRun import _AdcircRun
# from AdcircPy.Tides.TidalForcing import TidalForcing
# from AdcircPy.Winds.BestTrack import BestTrack


# class BestTrackRun(_AdcircRun):
#     def __init__(self, AdcircMesh, storm_id, start_date=None, end_date=None,
#                  spinup_days=None, TidalForcing=None, netcdf=True):










































#     self.storm_id = storm_id
#     self.start_date = start_date
#     self.end_date = end_date
#     self.spinup_date = spinup_date
#     self.TidalForcing = tides
#     self._init_fort22()
#     self._init_TidalForcing(kwargs.pop('constituents', None))
#     super(BestTrackRun, self).__init__(AdcircMesh, **kwargs)

# def _init_fort22(self):
#     self.fort22 = BestTrack(self.storm_id, self.start_date, self.end_date)

# def _init_TidalForcing(self, constituents):
#     if self.TidalForcing == True:
#         self.TidalForcing = TidalForcing(self.fort22.start_date, self.fort22.end_date, self.spinup_date, constituents=constituents)
#     else:
#         self.TidalForcing = None

# def _init_NWS(self):
#     self.NWS=20

# def _init_DRAMP(self):
#     if self.DRAMP is None:
#         if self.TidalForcing is not None:
#             self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)
#         else:
#             self.DRAMP = 0
# def _write_NRAMP(self):
#     if self.IHOT>0:
#         self.NRAMP=8
#         self.f.write('{:<32d}'.format(self.NRAMP))
#     else:
#         self.NRAMP=0
#         self.f.write('{:<32d}'.format(self.NRAMP))
# def _write_RNDAY(self):
#     if self.TidalForcing is not None:
#         if self.IHOT==0:
#             RNDAY = (self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
#         elif self.IHOT==567:
#             RNDAY = (self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
#     else:
#         RNDAY = (self.fort22.end_date - self.fort22.start_date).total_seconds()/(60*60*24)
#     self.f.write('{:<32.2f}'.format(RNDAY))
# def _write_DRAMP(self):
#     if self.NRAMP>0 and self.TidalForcing is not None:
#         if self.DUnRampMete is None:
#             self.DUnRampMete = (self.TidalForcing.start_date - self.TidalForcing.spinup_date).days
#         else:
#             self.DUnRampMete = 0
#         if self.DRAMPElev is None:
#             self.DRAMPElev = self.DRAMP
#         if self.DRAMPTip is None:
#             self.DRAMPTip = self.DRAMP
#         self.f.write('{:<4.1f} '.format(self.DRAMP))
#         self.f.write('{:<2}'.format(int(self.DRAMPExtFlux)))
#         self.f.write('{:<2}'.format(int(self.FluxSettlingTime)))
#         self.f.write('{:<2}'.format(int(self.DRAMPIntFlux)))
#         self.f.write('{:<3}'.format(int(self.DRAMPElev)))
#         self.f.write('{:<3}'.format(int(self.DRAMPTip)))
#         self.f.write('{:<4.1f}'.format(self.DRAMPMete))
#         self.f.write('{:<2}'.format(int(self.DRAMPWRad)))
#         self.f.write('{:<2}'.format(int(self.DUnRampMete)))
#         self.f.write('{:<7}'.format(''))
#     elif self.NRAMP==0 and self.TidalForcing is not None:
#         self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)
#         self.f.write('{:<32.1f}'.format(self.DRAMP))
#     else:
#         raise NotImplementedError('I suppose this is a met-only run. We need to code more.')
# def dump(self, path):
#     """
#     Overloads the parent's dump method such that if there is tidal forcing present
#     it can call the parent, otherwise the met only run requires a slightly different
#     fort.15 generation as met-only has no coldstart phase.
#     """
#     if self.TidalForcing is not None:
#         _AdcircRun.dump(self, path)
#     else:
#         self.IHOT=0
#         self.NRAMP=8
#         self._write_fort15()
#     self.fort22.dump(path)
