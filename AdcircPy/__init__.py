name = "AdcircPy"
from AdcircPy.Mesh   import Mesh
from AdcircPy.Outputs import Outputs

def read_mesh(fort14, datum=None, epsg=4326):
    """
    Reads ADCIRC input files.
    
    -----------
    Parameters:
    -----------
        fort14 (string) :  Path to fort.14 file.
    
    -------------------
    Optional arguments:
    -------------------
        fort13 (string) <None> : Path to optional fort.13 file.
        fort15 (string) <None> : Path to optional fort.15 file.
        datum  (string) <'MSL' (default)|'NAVD88'> : Mesh vertical datum.
    
    -------
    return:
    -------
        AdcirPy.Mesh instance.
    """
    return Mesh.init_from_fort14(fort14, datum, epsg)
    

def read_output(path, **kwargs):
    """
    Reads ADCIRC output files. Supports both ASCII and NetCDF.
    
    -----------
    Parameters:
    -----------
        path (str): Path to output file.
    
    -------------------
    Optional arguments:
    -------------------
        fort14 (str) <None> : Path to fort.14 file.
                              Required for ASCII gridded field outputs (e.g. maxele).
                              Optional in case of NetCDF files (used for inclusion of boundary data).
        fort13 (str) <None> : Path to optional fort.13 file.
        fort15 (str) <None> : Path to optional fort.15 file.
        datum (str) <'MSL'(default)|'NAVD88'> : Mesh vertical datum.
    
    -------
    return:
    -------
        AdcirPy.<output>  where <output> is the output type.
    """
    return Outputs.read_outputs(path, **kwargs)





















































# class Model():
#     def __init__(self, Mesh, Fort15):
#         self.Mesh   = Mesh
#         self.Fort15 = Fort15

#     def write_files()




    # @staticmethod
    # def _read_config_file(path)
    #     return _Fort15._read_config_file(path)
# class topobathy(Mesh):

#     def __init__(self, **kwargs):
#         super(topobathy, self).__init__(**kwargs)


#     def read_ascii_output(self, path, output_type=None):
#         """
#         Reads an ascii output file associated with the mesh instance.
        
#         args:
#             path (string) : Path to ascii output file
#         """
#         return outputs.read_ascii_output(path, self, output_type=output_type)

#     def write_fort14(self, path):
#         """ Writes the mesh to path. """
#         writters.write_fort14(self, path)

#     def interpolate_DEM(self, DEM, **kwargs):
#         """
#         Averages the DEM data contained on the Finite Volume surrounding the node and that
#         passes through the midpoints and centroids of the elements surrounding the node.
        
#         args:
#             DEM (object) : Must be a initialized instance of adcpy.demtools.DEM

#         Note: 
#             Function will modify values internally. It is recommended that a copy of the grid
#             instance is made before interpolation, so posterior comparisons can be made.
#         """

#         utilities.interpolate_DEM(self, DEM, **kwargs)

#     def make_plot(self, **kwargs):
#         """
#         Generates a ready-made plot with matplotlib.
        
#         kwargs:
#             extent  (list) <default None> : [minlon, minlat, maxlon, maxlat]
#             axes (matplotlib.Axes object) <default None> : )Pass user generated axes to plot.
#             title=None
#             shoreline=False
#             mesh=False
#             vmin=None
#             vmax=None
#         """
#         return plotters.plot_bathy(self, **kwargs)

#     def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
#         return raster.get_raster_from_geoTransform(self, geoTransform, shape, **kwargs)
    
#     def rasterize_to_extent(self, extent, dx, dy, epsg, padding=None):
#         return raster.get_raster_from_extent(self, extent, dx, dy, epsg, padding)

#     def plot_shoreline(self, **kwargs):
#         return plotters.plot_shoreline(self, **kwargs)

#     def msl_to_navd88(self, **kwargs):
#         datum.msl_to_navd88(self, **kwargs)

#     def navd88_to_msl(self, **kwargs):
#         datum.navd88_to_msl(self, **kwargs)

#     def deploy_to_PostGIS(self, dbname, **kwargs):
#         """
#         Writes the mesh file to a PostGIS database.
#         Database must be initialized on the server prior to deployment
#         and PostGIS support must be enabled on the database.
        
#         args:
#             dbname (string) : Name of the PostgreSQL database on which to deploy the mesh.

#         kwargs:
#             user (string)     <'postgres'> : Username to connect to database.
#                                               Make sure the user has write privileges on database.
#             password (string|bool)  <True> : If True password will be asked interactively. 
#                                              If False a passwordless login is assumed.
#                                              If string, literal string corresponds to password.
#                                              string literal passwords and passwordless users are not recommended.
#             host (string)    <'localhost'> : Domain name or IP address of Postgresql database.
#             port         (int)      <5432> : Postgresql port number.
#             progressbar (bool)     <False> : Shows a progressbar during deployment of elements to database.
#             overwrite   (bool)     <False> : Allow to dump existing tables and recreate them.
#             boundaries_only (bool) <False> : This argument skips the deployment of mesh elements.
        
#         Notes:
#             Deployment of mesh elements can take considerable time, depending on the size of the mesh.
#             You can turn off deployment of mesh elements by setting boundaries_only=True
#             Boundaries will be deployed and committed to the database before the elements, which will give
#             the user a chance to start exploring the data before the element deployment is completed.
#             The elements are deployed as "LinestringZ", they do not share sides or nodes, therefore do not
#             constitute a formal topology. A topological version has been tested and exists as a comment on
#             the source code. The main disadvantage of using the topological version is that deployment is 
#             extremely slow. Parallelization of element deployment is currently being considered, but so
#             far it looks like Postgresql does not support parallelization.
#         """
#         utilities.deploy_to_PostGIS(self, dbname, **kwargs)

# class outputSurface(Mesh):
#     def __init__(self, **kwargs):
#         self.time = kwargs.pop("time", None)
#         self.timestep = kwargs.pop("timestep", None)
#         super(outputSurface, self).__init__(**kwargs)

#     def make_plot(self, **kwargs):
#         return plotters.plot_surface(self, **kwargs)

#     def make_animation(self, **kwargs):
#         return plotters.surface_animation(self, **kwargs)
    
#     def __sub__(self, other):
#         return grid.get_difference(self, other)

#     def msl_to_navd88(self, datum_grid_path):
#         datum.msl_to_navd88(self, datum_grid_path)

#     def navd88_to_msl(self):
#         datum.navd88_to_msl(self)
    
# class stationTimeseries(object):
#     def __init__(self, **kwargs):
#         pass

# class surfaceDifference(Mesh):
#     def __init__(self, **kwargs):
#         super(surfaceDifference, self).__init__(**kwargs)
 
#     def make_plot(self, extent=None, axes=None, vmin=None, vmax=None, title=None, **kwargs):
#         plotters.plot_diff(self, extent, axes, vmin, vmax, title, **kwargs)


# # class vectorSurface(Mesh):
# #     self.time = kwargs.pop("time", None)
# #     self.timesteps = kwargs.pop("timesteps", None)
# #     self.bathymetry = kwargs.pop("bathymetry", None)
# #     super(vectorSurface, self).__init__(**kwargs)


# class ADCIRC():
#     def read_mesh(fort14, fort13=None, fort15=None, datum='MSL', epsg=4326):
#         """
#         Reads ADCIRC input files.
#         -----------
#         Parameters:
#         -----------
#             fort14 (string) :  Path to fort.14 file. 
#         -------------------
#         Optional arguments:
#         -------------------
#             fort13 (string) <None> : Path to optional fort.13 file.
#             fort15 (string) <None> : Path to optional fort.15 file.
#             datum  (string) <'MSL' (default)|'NAVD88'> : Mesh vertical datum.
#             epsg   (int)    <4326> : EPSG code for grid. 
#         -------
#         return:
#         -------
#             AdcirPy.adcirc.grid object instance.            
#         """
#         return _Mesh.init_from_files(fort14, fort13, fort15, datum, epsg)

#     def read_output(path, fort14=None, fort13=None, fort15=None, datum='MSL', epsg=4326):
#         """
#         Reads ADCIRC output files. Supports both ASCII and NetCDF.
#         -----------
#         Parameters:
#         -----------
#             path (str): Path to output file.
#         -------------------
#         Optional arguments:
#         -------------------
#             fort14 (str) <None> : Path to fort.14 file. Required for ASCII gridded field outputs (e.g. maxele).
#                                   Optional in case of NetCDF files (used for inclusion of boundary data).
#             fort13 (str) <None> : Path to optional fort.13 file.
#             fort15 (str) <None> : Path to optional fort.15 file.
#             datum (str) <'MSL'(default)|'NAVD88'> : Mesh vertical datum.
#             epsg (int) <4326> : EPSG code for grid. 
#         -------
#         return:
#         -------
#             AdcirPy.adcirc.<output> object where <output> depends on the type of output provided.
#         """
#         return adcpy.adcirc.outputs.read_outputs(path, fort14, fort13, fort15, datum, epsg)
