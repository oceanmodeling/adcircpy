# global imports
from pathlib import Path
from netCDF4 import Dataset

# local imports
from AdcircPy import Outputs


class _OutputFactory:

    def __init__(self, path, fort14=None, SpatialReference=4326,
                 vertical_datum=None, datum_mesh=None):
        self._path = path
        self._fort14 = fort14
        self._SpatialReference = SpatialReference
        self._vertical_datum = vertical_datum
        self._datum_mesh = datum_mesh

    def __netcdf_factory(self):
        nc = Dataset(self.path)
        if 'adcirc_mesh' in nc.variables.keys():
            if 'zeta_max' in nc.variables.keys():
                return Outputs.Maxele.from_netcdf(
                            self.path, self.fort14, self.vertical_datum,
                            self.SpatialReference)
            elif 'zeta' in nc.variables.keys():
                return Outputs.ElevationSurfaceTimeseries.from_netcdf(
                                self.path, self.fort14,
                                self.vertical_datum, self.SpatialReference,
                                self.datum_mesh)
            else:
                raise NotImplementedError(
                    'The OutputFactory class has not implemented this output '
                    + 'type yet, or this is not an Adcirc output file.')
        else:
            if 'zeta' in nc.variables.keys():
                return Outputs.ElevationStationsTimeseries.from_netcdf(
                                                                    self.path)

            else:
                raise NotImplementedError(
                 'The OutputFactory class has not implemented this output '
                 + 'type yet, or this is not an Adcirc output file.')

    def __ascii_factory(self):
        raise NotImplementedError('Ascii files not yet supported.')

    @property
    def output(self):
        try:
            return self.__netcdf_factory()
        except FileNotFoundError:
            raise
        except OSError:
            return self.__ascii_factory()

    @property
    def path(self):
        return str(self._path)

    @property
    def fort14(self):
        return self._fort14

    @property
    def SpatialReference(self):
        return self._SpatialReference

    @property
    def vertical_datum(self):
        return self._vertical_datum

    @property
    def datum_mesh(self):
        return self._datum_mesh

    @property
    def _path(self):
        return self.__path

    @property
    def _fort14(self):
        return self.__fort14

    @property
    def _SpatialReference(self):
        return self.__SpatialReference

    @property
    def _vertical_datum(self):
        return self.__vertical_datum

    @property
    def _datum_mesh(self):
        return self.__datum_mesh

    @_path.setter
    def _path(self, path):
        self.__path = Path(path)

    @_fort14.setter
    def _fort14(self, fort14):
        self.__fort14 = fort14

    @_SpatialReference.setter
    def _SpatialReference(self, SpatialReference):
        self.__SpatialReference = SpatialReference

    @_vertical_datum.setter
    def _vertical_datum(self, vertical_datum):
        self.__vertical_datum = vertical_datum

    @_datum_mesh.setter
    def _datum_mesh(self, datum_mesh):
        self.__datum_mesh = datum_mesh





    # def get_output_instance(self):
    #     if self._is_ncfile() is True:
    #         return self._netcdf_factory()
    #     else:
    #         return self._ascii_factory()

    # def _init_fort14(self):
    #     if isinstance(self.fort14, str):
    #         self.fort14 = AdcircMesh.from_fort14(
    #                         fort14=self.fort14,
    #                         vertical_datum=self.vertical_datum,
    #                         epsg=self.SpatialReference,
    #                         datum_grid=self.datum_mesh)

    # def _is_ncfile(self):
    #     try:
    #         Dataset(self.path)
    #         return True
    #     except:
    #         return False

    # def _netcdf_factory(self):
    #     nc = Dataset(self.path)
    #     if 'adcirc_mesh' in nc.variables.keys():
    #         if 'zeta_max' in nc.variables.keys():
    #             return Maxele.from_netcdf(self.path, self.fort14,
    #                                       self.vertical_datum, self.SpatialReference,
    #                                       self.datum_mesh)
    #         elif 'zeta' in nc.variables.keys():
    #             return ElevationSurfaceTimeseries.from_netcdf(
    #                             self.path, self.fort14,
    #                             self.vertical_datum, self.SpatialReference,
    #                             self.datum_mesh)
    #         else:
    #             raise NotImplementedError(
    #                                 'The OutputFactory class has not'
    #                                 + 'implemented this output type yet,'
    #                                 + 'or this is not an Adcirc output file.')

    #     else:
    #         if 'zeta' in nc.variables.keys():
    #             return ElevationStations.from_netcdf(self.path)
    #         elif 'phs' in nc.variables.keys():
    #             return HarmonicConstituentsElevationStations.from_netcdf(
    #                                                                 self.path)
    #         else:
    #             raise NotImplementedError(
    #              'The OutputFactory class has not implemented this output '
    #              + 'type yet, or this is not an Adcirc output file.')

    # def _ascii_factory(self):
    #     self.f = open(self.path, 'r')
    #     # Might be a harmonic constituents file...
    #     if self._check_is_HarmonicConstituentsFile()==True:
    #         number_of_points = self._check_is_SurfaceOrStations()
    #         if self.fort14 is not None and self.fort14.x.size==number_of_points:
    #             raise NotImplementedError('Guessed an ASCII HarmonicConstituentsOutputSurface but instantiantion has not yet been implemented.')
    #         elif self.fort15 is not None:
    #             # should probably parse fort.15 and compare the number of stations first.
    #             return HarmonicConstituentsElevationStations.from_ascii(self.path, self.fort15)
    #         else:
    #             raise Exception('When loading an ASCII Harmonic constituents file, the fort.14 needs to be provided for fort.53 files or the fort.15 needs to be provided for fort.51 files.')    

    #     # else, could be a general surface file...
    #     elif self._check_is_SurfaceFile()==True:
    #         if self.number_of_points == self.fort14.x.size:
    #             # Could be a surface time series
    #             if self.number_of_datasets > 2:
    #                 raise NotImplementedError("Guessed Scalar Surface Timeseries output, but instantiantion is not yet implemented")
    #             # or a surface maxima file
    #             elif self.number_of_datasets in [1,2]:
    #                 return self._init_ScalarSurfaceMaxima()

    #     # else has to be a stations output file.
    #     else:
    #         raise NotImplementedError(
    #               'Guessed an Output Stations file but instantiantion has not yet been implemented.')

    # def _check_is_HarmonicConstituentsFile(self):
    #     self.line = self.f.readline().strip()
    #     try:
    #         int(self.line)
    #         return True
    #     except:
    #         return False

    # def _check_is_SurfaceOrStations(self):
    #     for i in range(int(self.line)):
    #         self.f.readline()
    #     return int(self.f.readline())

    # def _check_is_SurfaceFile(self):
    #     # we keep pushing this check further and further below...
    #     if self.fort14 is None:
    #         raise Exception(
    #                'A fort.14 is required for reading some output file types.')
    #     self.line = self.f.readline().split()
    #     self.number_of_datasets = int(self.line[0])
    #     self.number_of_points = int(self.line[1])
    #     if self.number_of_points == self.fort14.x.size:
    #         return True

    # def _init_ScalarSurfaceMaxima(self):
    #     """ """
    #     self.f.readline().split()
    #     nodeID = list()
    #     values = list()
    #     for i in range(self.number_of_points):
    #         self.line = self.f.readline().split()
    #         nodeID.append(int(self.line[0].strip(' \n')))
    #         values.append(float(self.line[1].strip(' \n')))
    #     extrema_time_vector = list()
    #     if self.number_of_datasets == 2:
    #         for i in range(self.number_of_points):
    #             self.line = self.f.readline().split()
    #             extrema_time_vector.append(float(self.line[1].strip(' \n')))
    #     nodeID = np.asarray(nodeID)
    #     values = np.ma.masked_equal(values, -99999.)
    #     return ScalarSurfaceExtrema(x=self.fort14.x,
    #                                 y=self.fort14.y,
    #                                 elements=self.fort14.elements,
    #                                 values=values,
    #                                 epsg=self.SpatialReference,
    #                                 vertical_datum=self.vertical_datum,
    #                                 times=extrema_time_vector,
    #                                 nodeID=nodeID)

    # def __del__(self):
    #     if hasattr(self, 'f'):
    #         self.f.close()
