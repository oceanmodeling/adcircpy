# globa imports
import matplotlib as mpl
import numpy as np
from pyproj import Proj, transform
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

# local imports
from AdcircPy.Mesh import Boundaries
from AdcircPy.Mesh.Boundaries._ModelDomain import _ModelDomain
from AdcircPy._FixPointNormalize import _FixPointNormalize


# unittest imports
import unittest

mpl.rcParams['agg.path.chunksize'] = 10000


class UnstructuredMesh(object):

    def __init__(self, xy=None, z=None, elements=None, epsg=None, nodeId=None,
                 elementId=None, vertical_datum=None, OceanBoundaries=None,
                 LandBoundaries=None, InnerBoundaries=None,
                 InflowBoundaries=None, OutflowBoundaries=None,
                 WeirBoundaries=None, CulvertBoundaries=None):
        super(UnstructuredMesh, self).__init__()
        self._xy = xy
        self._z = z
        self._elements = elements
        self._epsg = epsg
        self._OceanBoundaries = OceanBoundaries
        self._LandBoundaries = LandBoundaries
        self._InnerBoundaries = InnerBoundaries
        self._InflowBoundaries = InflowBoundaries
        self._OutflowBoundaries = OutflowBoundaries
        self._WeirBoundaries = WeirBoundaries
        self._CulvertBoundaries = CulvertBoundaries
        self.__init_ModelDomain()

    def get_x(self, epsg=None):
        if self.xy is not None:
            if epsg is not None:
                if epsg != self.epsg:
                    return self.__get_transformed_array(epsg, self.x)
            else:
                return self.x

    def get_y(self, epsg=None):
        if self.y is not None:
            if epsg is not None:
                if epsg != self.epsg:
                    return self.__get_transformed_array(epsg, self.y)
            else:
                return self.y

    def get_xy(self, epsg=None):
        if self.xy is not None:
            if epsg is not None:
                if epsg != self.epsg:
                    return self.__get_transformed_array(epsg, self.xy)
            else:
                return self.xy

    def get_xyz(self, epsg=None):
        if self.xyz is not None:
            if epsg is not None:
                if epsg != self.epsg:
                    xy = self.__get_transformed_array(epsg, self.xy)
                    return np.c_[xy, self.z]
            else:
                return self.xyz

    def get_boundaries(self, epsg=None, bbox=None, h0=None):
        return self.ModelDomain.get_boundaries(epsg, bbox, h0)

    def get_outer_boundary(self, epsg=None, bbox=None, h0=None):
        return self.ModelDomain.get_outer_boundary(epsg, bbox, h0)

    def get_inner_boundaries(self, epsg=None, bbox=None, h0=None):
        return self.ModelDomain.get_inner_boundaries(epsg, bbox, h0)

    def transform(self, epsg):
        if self.epsg is None:
            raise Exception("No projection defined on mesh. "
                            + "Set Projection before running transform()")
        if epsg != self.epsg:
            self._xy = self.__get_transformed_array(epsg, self.xy)
            self._epsg = epsg

    def make_plot(self, elements=True, axes=None, vmin=None, vmax=None,
                  cmap=None, levels=None, show=False, title=None, figsize=None,
                  colors=None, extent=None, cbar_label=None, norm=None,
                  tricontourf_kwargs=dict(), triplot_kwargs=dict()):
        if axes is None:
            axes = plt.figure(figsize=figsize).add_subplot(111)
        if vmin is None:
            vmin = np.ceil(np.min(self.z))
        if vmax is None:
            vmax = np.ceil(np.max(self.z))
        if colors is None:
            colors = 256
        else:
            assert isinstance(colors, int)
        if cmap is None:
            cmap = plt.cm.get_cmap('viridis')
            if levels is None:
                levels = np.linspace(vmin, vmax, colors)
            _col_val = 0.
        elif cmap == 'topobathy':
            if vmax < 0.:
                cmap = plt.cm.seismic
                _col_val = 0.
            else:
                wet_count = int(np.floor(colors*(float((self.z < 0.).sum())
                                                 / float(self.z.size))))
                _col_val = float(wet_count)/colors
                dry_count = colors - wet_count
                colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
                colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
                colors = np.vstack((colors_undersea, colors_land))
                cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
                wlevels = np.linspace(vmin, 0.0, wet_count, endpoint=False)
                dlevels = np.linspace(0.0, vmax, dry_count)
                levels = np.hstack((wlevels, dlevels))
        else:
            cmap = plt.cm.get_cmap(cmap)
            levels = np.linspace(vmin, vmax, colors)
            _col_val = 0.
        if norm is None:
            norm = _FixPointNormalize(sealevel=0.0, vmax=vmax, vmin=vmin,
                                      col_val=_col_val)
        if tricontourf_kwargs is None:
            tricontourf_kwargs = dict()
        axes.tricontourf(self.x, self.y, self.elements, self.z, levels=levels,
                         cmap=cmap, norm=norm)
        if elements is True:
            if triplot_kwargs is None:
                triplot_kwargs = dict()
            color = tricontourf_kwargs.pop('color', 'k')
            alpha = tricontourf_kwargs.pop('alpha', 0.25)
            linewidth = tricontourf_kwargs.pop('linewidth', 0.25)
            axes.triplot(self.x, self.y, self.elements, color=color,
                         linewidth=linewidth, alpha=alpha, **triplot_kwargs)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        mappable = ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(mappable, cax=cax,  # extend=cmap_extend,
                            orientation='horizontal')
        if _col_val != 0:
            cbar.set_ticks([vmin, vmin + _col_val * (vmax-vmin), vmax])
            cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
        else:
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        if show is True:
            plt.show()
        return axes

    def __init_x(self):
        if self.xy is not None:
            self.__x = self.xy[:, 0]
        else:
            self.__x = None

    def __init_y(self):
        if self.xy is not None:
            self.__y = self.xy[:, 1]
        else:
            self.__y = None

    def __init_xyz(self):
        if self.z is not None and self.xy is not None:
            self.__xyz = np.vstack([self.x, self.y, self.z]).T
        else:
            self.__xyz = None

    def __init_ModelDomain(self):
        ExternalBoundaries = list(filter(lambda x: x is not None,
                                         [self.OceanBoundaries,
                                          self.LandBoundaries,
                                          self.InflowBoundaries,
                                          self.OutflowBoundaries]))
        InternalBoundaries = list(filter(lambda x: x is not None,
                                         [self.InnerBoundaries,
                                          self.WeirBoundaries,
                                          self.CulvertBoundaries]))
        if len(ExternalBoundaries) > 0:
            self.__ModelDomain = _ModelDomain(
                                        self.epsg,
                                        *ExternalBoundaries,
                                        InternalBoundaries=InternalBoundaries)
        else:
            self.__ModelDomain = None

    def __get_transformed_array(self, target_epsg, xy):
        target_proj = Proj(init="epsg:{}".format(target_epsg))
        source_Proj = Proj(init="epsg:{}".format(self.epsg))
        x, y = transform(source_Proj, target_proj, xy[:, 0], xy[:, 1])
        return np.vstack([x, y]).T

    @property
    def xy(self):
        return self._xy

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def xyz(self):
        return self.__xyz

    @property
    def elements(self):
        return self._elements

    @property
    def epsg(self):
        return self._epsg

    @property
    def BoundaryCollection(self):
        return self.__BoundaryCollection

    @property
    def OceanBoundaries(self):
        return self._OceanBoundaries

    @property
    def LandBoundaries(self):
        return self._LandBoundaries

    @property
    def InnerBoundaries(self):
        return self._InnerBoundaries

    @property
    def InflowBoundaries(self):
        return self._InflowBoundaries

    @property
    def OutflowBoundaries(self):
        return self._OutflowBoundaries

    @property
    def WeirBoundaries(self):
        return self._WeirBoundaries

    @property
    def CulvertBoundaries(self):
        return self._CulvertBoundaries

    @property
    def ModelDomain(self):
        return self._ModelDomain

    @property
    def plt(self):
        return self._plt

    @property
    def _xy(self):
        return self.__xy

    @property
    def _x(self):
        return self.__x

    @property
    def _y(self):
        return self.__y

    @property
    def _z(self):
        return self.__z

    @property
    def _xyz(self):
        return self.__xyz

    @property
    def _elements(self):
        return self.__elements

    @property
    def _epsg(self):
        return self.__epsg

    @property
    def _OceanBoundaries(self):
        return self.__OceanBoundaries

    @property
    def _LandBoundaries(self):
        return self.__LandBoundaries

    @property
    def _InnerBoundaries(self):
        return self.__InnerBoundaries

    @property
    def _InflowBoundaries(self):
        return self.__InflowBoundaries

    @property
    def _OutflowBoundaries(self):
        return self.__OutflowBoundaries

    @property
    def _WeirBoundaries(self):
        return self.__WeirBoundaries

    @property
    def _CulvertBoundaries(self):
        return self.__CulvertBoundaries

    @property
    def _ModelDomain(self):
        return self.__ModelDomain

    @_xy.setter
    def _xy(self, xy):
        if xy is not None:
            xy = np.asarray(xy)
            assert xy.ndim == 2
            self.__xy = xy
        else:
            self.__xy = None
        self.__init_x()
        self.__init_y()

    @_z.setter
    def _z(self, z):
        if z is not None:
            z = np.asarray(z)
            assert z.shape != self.xy.shape[0]
        self.__z = z
        self.__init_xyz()



    @_elements.setter
    def _elements(self, elements):
        if elements is not None:
            elements = np.asarray(elements)
            assert elements.shape[1] == 3
            self.__elements = elements
        else:
            self.__elements = None

    @_epsg.setter
    def _epsg(self, _epsg):
        if _epsg is not None:
            assert isinstance(_epsg, int)
            self.__epsg = _epsg
        else:
            self.__epsg = None

    @_OceanBoundaries.setter
    def _OceanBoundaries(self, OceanBoundaries):
        if OceanBoundaries is not None:
            assert isinstance(OceanBoundaries, Boundaries.OceanBoundaries)
        else:
            OceanBoundaries = Boundaries.OceanBoundaries()
        self.__OceanBoundaries = OceanBoundaries

    @_LandBoundaries.setter
    def _LandBoundaries(self, LandBoundaries):
        if LandBoundaries is not None:
            assert isinstance(LandBoundaries, Boundaries.LandBoundaries)
        else:
            LandBoundaries = Boundaries.LandBoundaries()
        self.__LandBoundaries = LandBoundaries

    @_InnerBoundaries.setter
    def _InnerBoundaries(self, InnerBoundaries):
        if InnerBoundaries is not None:
            assert isinstance(InnerBoundaries, Boundaries.InnerBoundaries)
        else:
            InnerBoundaries = Boundaries.InnerBoundaries()
        self.__InnerBoundaries = InnerBoundaries

    @_InflowBoundaries.setter
    def _InflowBoundaries(self, InflowBoundaries):
        if InflowBoundaries is not None:
            assert isinstance(InflowBoundaries, Boundaries.InflowBoundaries)
        else:
            InflowBoundaries = Boundaries.InflowBoundaries()
        self.__InflowBoundaries = InflowBoundaries

    @_OutflowBoundaries.setter
    def _OutflowBoundaries(self, OutflowBoundaries):
        if OutflowBoundaries is not None:
            assert isinstance(OutflowBoundaries, Boundaries.OutflowBoundaries)
        else:
            OutflowBoundaries = Boundaries.OutflowBoundaries()
        self.__OutflowBoundaries = OutflowBoundaries

    @_WeirBoundaries.setter
    def _WeirBoundaries(self, WeirBoundaries):
        if WeirBoundaries is not None:
            assert isinstance(WeirBoundaries, Boundaries.WeirBoundaries)
        else:
            WeirBoundaries = Boundaries.WeirBoundaries()
        self.__WeirBoundaries = WeirBoundaries

    @_CulvertBoundaries.setter
    def _CulvertBoundaries(self, CulvertBoundaries):
        if CulvertBoundaries is not None:
            assert isinstance(CulvertBoundaries, Boundaries.CulvertBoundaries)
        else:
            CulvertBoundaries = Boundaries.CulvertBoundaries()
        self.__CulvertBoundaries = CulvertBoundaries


class UnstructuredMeshTestCase(unittest.TestCase):

    def setUp(self):
        self.UnstructuredMesh = UnstructuredMesh

    def test_empty(self):
        self.UnstructuredMesh()
