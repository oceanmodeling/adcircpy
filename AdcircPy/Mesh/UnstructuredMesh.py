# global imports
import numpy as np
from osgeo import osr
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
from scipy.spatial import cKDTree

# local imports
from AdcircPy.Mesh import Boundaries, NodalAttributes
from AdcircPy._FixPointNormalize import _FixPointNormalize


# unittest imports
import unittest
import os


class UnstructuredMesh(object):

    def __init__(self, xy=None, z=None, elements=None, SpatialReference=None,
                 node_id=None, element_id=None, vertical_datum=None,
                 OceanBoundaries=None, LandBoundaries=None,
                 InnerBoundaries=None, InflowBoundaries=None,
                 OutflowBoundaries=None, WeirBoundaries=None,
                 CulvertBoundaries=None, description='', mask_value=-99999.):
        # this needs to be improved.
        # these attributes must exist befor anything else can happen
        self.__OceanBoundaries = None
        self.__LandBoundaries = None
        self.__InnerBoundaries = None
        self.__InflowBoundaries = None
        self.__OutflowBoundaries = None
        self.__WeirBoundaries = None
        self.__CulvertBoundaries = None
        self._xy = xy
        self._mask_value = mask_value
        self._z = z
        self._elements = elements
        self._SpatialReference = SpatialReference
        self._OceanBoundaries = OceanBoundaries
        self._LandBoundaries = LandBoundaries
        self._InnerBoundaries = InnerBoundaries
        self._InflowBoundaries = InflowBoundaries
        self._OutflowBoundaries = OutflowBoundaries
        self._WeirBoundaries = WeirBoundaries
        self._CulvertBoundaries = CulvertBoundaries
        self._description = description
        self.__set_OuterRing()
        self.__set_InnerRings()
        self.__set_ConvexHull()
        self.__set_NodalAttributes()
        self.__set_KDTree()

    def get_x(self, SpatialReference=None):
        if self.xy is not None:
            if SpatialReference is not None:
                return self.__get_transformed_array(SpatialReference,
                                                    self.xy)[:, 0]
            else:
                return self.x

    def get_y(self, SpatialReference=None):
        if self.y is not None:
            if SpatialReference is not None:
                return self.__get_transformed_array(SpatialReference,
                                                    self.xy)[:, 1]
            else:
                return self.y

    def get_xy(self, SpatialReference=None):
        if self.xy is not None:
            if SpatialReference is not None:
                return self.__get_transformed_array(SpatialReference, self.xy)
            else:
                return self.xy

    def get_xyz(self, SpatialReference=None):
        if self.xyz is not None:
            xy = self.get_xy(SpatialReference)
            return np.c_[xy, self.z]
        else:
            return self.xyz

    def get_extent(self, SpatialReference=None):
        xy = self.get_xy(SpatialReference)
        return (np.min(xy[:, 0]), np.max(xy[:, 0]),
                np.min(xy[:, 1]), np.max(xy[:, 1]))

    # def get_boundaries(self, SpatialReference=None, bbox=None, h0=None):
    #     return self.ModelDomain.get_boundaries(epsg, bbox, h0)

    # def get_outer_boundary(self, SpatialReference=None, bbox=None, h0=None):
    #     return self.ModelDomain.get_outer_boundary(epsg, bbox, h0)

    # def get_inner_boundaries(self, SpatialReference=None, bbox=None, h0=None):
    #     return self.ModelDomain.get_inner_boundaries(epsg, bbox, h0)

    def transform(self, SpatialReference):
        if self.SpatialReference is None:
            raise Exception("No spatial reference defined on mesh. "
                            + "must instantiate with SpatialReference kwarg "
                            + "in order to be able to use transform()")
        if not self.SpatialReference.IsSame(SpatialReference):
            self._xy = self.__get_transformed_array(SpatialReference, self.xy)
            self._SpatialReference = SpatialReference

    def add_boundary(self, ibtype, **kwargs):
        """
        Boundary type is identified by ibtype argument. Pass literal None
        for ocean boundaries. All other boundaries follow ADCIRC
        specifications found here:
        https://adcirc.org/home/documentation/users-manual-v53/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14/
        """
        if self.SpatialReference is None:
            raise Exception('add_boundary() not supported for meshes with no '
                            + 'spatial reference.')
        assert isinstance(ibtype, type(None)) or isinstance(ibtype, int)
        if ibtype is None:
            name = 'ocean boundary'
            kwargs = self.__get_boundary_kwargs(['indexes'], kwargs, name,
                                                ibtype)
            self.OceanBoundaries._add_boundary(**kwargs)
            self.__set_OuterRing()
        elif ibtype in [0, 10, 20]:
            name = 'land boundary'
            kwargs = self.__get_boundary_kwargs(['indexes'], kwargs, name,
                                                ibtype)
            self.LandBoundaries._add_boundary(**kwargs)
            self.__set_OuterRing()
        elif ibtype in [1, 11, 21]:
            name = 'inner boundary'
            kwargs = self.__get_boundary_kwargs(['indexes'], kwargs, name,
                                                ibtype)
            self.InnerBoundaries._add_boundary(**kwargs)
            self.__set_InnerRings()
        elif ibtype in [2, 12, 22, 102, 122]:
            name = 'inflow boundary'
            kwargs = self.__get_boundary_kwargs(['indexes'], kwargs, name,
                                                ibtype)
            self.InflowBoundaries._add_boundary(**kwargs)
        elif ibtype in [3, 13, 23]:
            name = 'outflow boundary'
            kwargs = self.__get_boundary_kwargs([
                                        'indexes',
                                        'barrier_heights',
                                        'supercritical_flow_coefficients'],
                                        kwargs, name, ibtype)
            self.OutflowBoundaries._add_boundary(**kwargs)
            self.__set_InnerRings()
        elif ibtype in [4, 24]:
            name = 'weir boundary'
            kwargs = self.__get_boundary_kwargs([
                                        'front_face_indexes',
                                        'back_face_indexes',
                                        'barrier_heights',
                                        'supercritical_flow_coefficients',
                                        'subcritical_flow_coefficients'],
                                        kwargs, name, ibtype)
            self.WeirBoundaries._add_boundary(**kwargs)
            self.__set_InnerRings()
        elif ibtype in [5, 25]:
            name = 'culvert boundary'
            kwargs = self.__get_boundary_kwargs([
                                        'front_face_indexes',
                                        'back_face_indexes',
                                        'barrier_heights',
                                        'supercritical_flow_coefficients',
                                        'subcritical_flow_coefficients',
                                        'cross_barrier_pipe_heights',
                                        'friction_factors',
                                        'pipe_diameters'],
                                        kwargs, name, ibtype)
            self.CulvertBoundaries._add_boundary(**kwargs)
            self.__set_InnerRings()
        else:
            raise TypeError("Unrecognized ibtype={}.".format(ibtype)
                            + "Check https://adcirc.org/home/documentation/"
                            + "users-manual-v53/input-file-descriptions/"
                            + "adcirc-grid-and-boundary-information-"
                            + "file-fort-14/ for valid ibtypes")

    def make_plot(self, elements=False, axes=None, vmin=None, vmax=None,
                  cmap=None, levels=None, show=False, title=None, figsize=None,
                  colors=256, extent=None, cbar_label=None, norm=None,
                  tricontourf_kwargs=dict(), triplot_kwargs=dict()):
        if axes is None:
            axes = plt.figure(figsize=figsize).add_subplot(111)
        if vmin is None:
            vmin = np.ceil(np.min(self.z))
        if vmax is None:
            vmax = np.ceil(np.max(self.z))
        cmap, norm, levels, col_val = self.__get_cmap(vmin, vmax, cmap, levels,
                                                      colors, norm)
        if tricontourf_kwargs is None:
            tricontourf_kwargs = dict()
        axes.tricontourf(self.x, self.y, self.elements, self.z, levels=levels,
                         cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
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
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(mappable, cax=cax,  # extend=cmap_extend,
                            orientation='horizontal')
        if col_val != 0:
            cbar.set_ticks([vmin, vmin + col_val * (vmax-vmin), vmax])
            cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
        else:
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        if show is True:
            plt.show()
        return axes

    def set_nodal_attribute_state(self, attribute_name, coldstart=True,
                                  hotstart=True):
        self.NodalAttributes.set_attribute_state(attribute_name,
                                                 coldstart=coldstart,
                                                 hotstart=hotstart)

    def set_xy(self, xy, elements=None):
        self._xy = xy
        self.set_elements(elements)

    def set_elements(self, elements):
        if self.xy is None:
            raise Exception('Cannot set elements on empty mesh.')
        self._elements = elements

    def set_SpatialReference(self, SpatialReference):
        self._SpatialReference = SpatialReference

    def set_node_id(self, node_id):
        if self.xy is None:
            raise Exception('Cannot set node_id on empty mesh.')
        self._node_id = node_id

    def set_element_id(self, element_id):
        if self.elements is None:
            raise Exception('Cannot set element_id on mesh with no '
                            + 'elements.')
        self._element_id = element_id

    def set_vertical_datum(self, vertical_datum):
        self._vertical_datum = vertical_datum

    def set_description(self, description):
        self._description = description

    def get_element_containing_coord(self, coord):
        x = coord[0]
        y = coord[1]
        distance, node_idx = self.KDTree.query([x, y])
        elements = self.get_finitive_volume_as_Path_list(node_idx)
        for i, element in enumerate(elements):
            if element.contains_point((x, y)):
                return self.get_finitive_volume_indexes(node_idx)[i]

    def get_finitive_volume_as_Path_list(self, index):
        return [self.get_Path_from_element(element) for element
                in self.get_finitive_volume_indexes(index)]

    def get_finitive_volume_indexes(self, index):
        return self.elements[np.where(np.any(np.isin(self.elements, index),
                                             axis=1))[0]]

    def get_Path_from_element(self, element):
        return Path([[self.x[element[0]], self.y[element[0]]],
                     [self.x[element[1]], self.y[element[1]]],
                     [self.x[element[2]], self.y[element[2]]],
                     [self.x[element[0]], self.y[element[0]]]], closed=True)

    def __set_xyz(self):
        if self.z is not None and self.xy is not None:
            self.__xyz = np.vstack([self.x, self.y, self.z]).T
        else:
            self.__xyz = None

    def __set_OuterRing(self):
        if self.SpatialReference is not None:
            ExternalBoundaries = list(filter(lambda x: x is not None,
                                             [self.OceanBoundaries,
                                              self.LandBoundaries,
                                              self.InflowBoundaries,
                                              self.OutflowBoundaries]))
            open_boundaries = list()
            for _BaseBoundary in ExternalBoundaries:
                for boundary in _BaseBoundary:
                    SpatialReference = boundary['SpatialReference']
                    vertices = boundary['vertices']
                    if not self.SpatialReference.IsSame(SpatialReference):
                        vertices = self.__transform_vertices(vertices,
                                                             SpatialReference)
                    open_boundaries.append(vertices)
            if len(open_boundaries) > 0:
                ordered_vertices = open_boundaries.pop()
                while open_boundaries:
                    # Take tail of ordered_vertices.
                    x0, y0 = ordered_vertices[-1]
                    # Find a matching head.
                    heads = [boundary[0, :] for boundary in open_boundaries]
                    heads_tree = cKDTree(heads)
                    head_dist, hidx = heads_tree.query([x0, y0])
                    # Find a matching tail
                    tails = [boundary[-1, :] for boundary in open_boundaries]
                    tails_tree = cKDTree(tails)
                    tail_dist, tidx = tails_tree.query([x0, y0])
                    if head_dist < tail_dist:
                        ordered_vertices = np.vstack(
                                                [ordered_vertices,
                                                 open_boundaries.pop(hidx)])
                    else:
                        ordered_vertices = np.vstack([
                                        ordered_vertices,
                                        np.flipud(open_boundaries.pop(tidx))])
                self.__OuterRing = Path(ordered_vertices, closed=True)
            else:
                self.__OuterRing = None
        else:
            self.__OuterRing = None

    def __set_InnerRings(self):
        if self.SpatialReference is not None:
            InternalBoundaries = list(filter(lambda x: x is not None,
                                             [self.InnerBoundaries,
                                              self.WeirBoundaries,
                                              self.CulvertBoundaries]))
            inner_boundaries = list()
            for _BaseBoundary in InternalBoundaries:
                if isinstance(_BaseBoundary, Boundaries.InnerBoundaries):
                    for boundary in _BaseBoundary:
                        SpatialReference = boundary['SpatialReference']
                        vertices = boundary['vertices']
                        if not self.SpatialReference.IsSame(SpatialReference):
                            vertices = self.__transform_vertices(
                                                            vertices,
                                                            SpatialReference)
                        inner_boundaries.append(vertices)
                elif isinstance(_BaseBoundary, (Boundaries.WeirBoundaries,
                                                Boundaries.CulvertBoundaries)):
                    for boundary in _BaseBoundary:
                        SpatialReference = boundary['SpatialReference']
                        front_face_vertices = boundary['front_face_vertices']
                        back_face_vertices = boundary['back_face_vertices']
                        if not self.SpatialReference.IsSame(SpatialReference):
                            front_face_vertices = self.__transform_vertices(
                                                        front_face_vertices,
                                                        SpatialReference)
                            back_face_vertices = self.__transform_vertices(
                                                        back_face_vertices,
                                                        SpatialReference)
                        vertices = np.vstack([front_face_vertices,
                                              np.flipud(back_face_vertices)])
                        inner_boundaries.append(vertices)
            if len(inner_boundaries) > 0:
                self.__InnerRings = [Path(vertices, closed=True) for vertices
                                     in inner_boundaries]
        else:
            self.__InnerRings = list()

        if len(inner_boundaries) == 0:
            self.__InnerRings = list()

    def __set_ConvexHull(self):
        self.__ConvexHull = (self.OuterRing, self.InnerRings)

    def __set_NodalAttributes(self):
        self.__NodalAttributes = NodalAttributes(self)

    def __set_KDTree(self):
        self.__KDTree = cKDTree(self.xy)

    def __get_transformed_array(self, SpatialReference, xy):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        if not self.SpatialReference.IsSame(SpatialReference):
            CoordinateTransform = osr.CoordinateTransformation(
                                                        self.SpatialReference,
                                                        SpatialReference)
            xy = [(x, y) for x, y in xy]
            xy = CoordinateTransform.TransformPoints(xy)
            xy = np.asarray([(x, y) for x, y, _ in xy])
        return xy

    def __get_boundary_kwargs(self, keys, kwargs, name, ibtype):
        kwargs['UnstructuredMesh'] = self
        kwargs['ibtype'] = ibtype
        for key in keys:
            try:
                kwargs[key] = kwargs.pop('{}'.format(key))
            except KeyError:
                raise Exception("add_boundary() missing "
                                + "'{}' ".format(key)
                                + " in kwargs for ibtype="
                                + "{} ".format(ibtype)
                                + "({})".format(name))
        return kwargs

    def __get_cmap(self, vmin, vmax, cmap=None, levels=None, colors=256,
                   norm=None):
        colors = int(colors)
        if cmap is None:
            cmap = plt.cm.get_cmap('jet')
            if levels is None:
                levels = np.linspace(vmin, vmax, colors)
            col_val = 0.
        elif cmap == 'topobathy':
            if vmax <= 0.:
                cmap = plt.cm.seismic
                col_val = 0.
                levels = np.linspace(vmin, vmax, colors)
            else:
                wet_count = int(np.floor(colors*(float((self.z < 0.).sum())
                                                 / float(self.z.size))))
                col_val = float(wet_count)/colors
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
            col_val = 0.
        if vmax > 0:
            if norm is None:
                norm = _FixPointNormalize(sealevel=0.0, vmax=vmax, vmin=vmin,
                                          col_val=col_val)
        return cmap, norm, levels, col_val

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
    def SpatialReference(self):
        return self._SpatialReference

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
    def OuterRing(self):
        return self._OuterRing

    @property
    def InnerRings(self):
        return self._InnerRings

    @property
    def ConvexHull(self):
        return self.__ConvexHull

    @property
    def NodalAttributes(self):
        return self._NodalAttributes

    @property
    def description(self):
        return self._description

    @property
    def KDTree(self):
        return self.__KDTree

    @property
    def mask_value(self):
        return self._mask_value

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
    def _SpatialReference(self):
        return self.__SpatialReference

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
    def _OuterRing(self):
        return self.__OuterRing

    @property
    def _InnerRings(self):
        return self.__InnerRings

    @property
    def _NodalAttributes(self):
        return self.__NodalAttributes

    @property
    def _description(self):
        return self.__description

    @property
    def _mask_value(self):
        return self.__mask_value

    @_xy.setter
    def _xy(self, xy):
        if xy is not None:
            xy = np.asarray(xy)
            assert xy.shape[1] == 2
            self.__xy = xy
            self._x = self.xy[:, 0]
            self._y = self.xy[:, 1]
        else:
            self.__xy = None
            self._x = None
            self._y = None

    @_x.setter
    def _x(self, x):
        self.__x = x

    @_y.setter
    def _y(self, y):
        self.__y = y

    @_z.setter
    def _z(self, z):
        if z is not None:
            z = np.asarray(z)
            assert z.shape != self.xy.shape[0]
        z = np.ma.masked_equal(z, self.mask_value)
        self.__z = z
        self.__set_xyz()

    @_elements.setter
    def _elements(self, elements):
        if elements is not None:
            elements = np.asarray(elements)
            assert elements.shape[1] == 3
            self.__elements = elements
        else:
            self.__elements = None

    @_SpatialReference.setter
    def _SpatialReference(self, SpatialReference):
        if SpatialReference is not None:
            if isinstance(SpatialReference, int):
                EPSG = SpatialReference
                SpatialReference = osr.SpatialReference()
                SpatialReference.ImportFromEPSG(EPSG)
            elif isinstance(SpatialReference, str):
                try:
                    EPSG = int(SpatialReference)
                    SpatialReference = osr.SpatialReference()
                    SpatialReference.ImportFromEPSG(EPSG)
                except ValueError:
                    raise
            else:
                assert isinstance(SpatialReference, osr.SpatialReference)
        self.__SpatialReference = SpatialReference

    @_OceanBoundaries.setter
    def _OceanBoundaries(self, OceanBoundaries):
        self.__OceanBoundaries = Boundaries.OceanBoundaries()
        if OceanBoundaries is not None:
            OceanBoundaries = list(OceanBoundaries)
            for boundary in OceanBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(None, **boundary)

    @_LandBoundaries.setter
    def _LandBoundaries(self, LandBoundaries):
        self.__LandBoundaries = Boundaries.LandBoundaries()
        if LandBoundaries is not None:
            LandBoundaries = list(LandBoundaries)
            for boundary in LandBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_InnerBoundaries.setter
    def _InnerBoundaries(self, InnerBoundaries):
        self.__InnerBoundaries = Boundaries.InnerBoundaries()
        if InnerBoundaries is not None:
            InnerBoundaries = list(InnerBoundaries)
            for boundary in InnerBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_InflowBoundaries.setter
    def _InflowBoundaries(self, InflowBoundaries):
        self.__InflowBoundaries = Boundaries.InflowBoundaries()
        if InflowBoundaries is not None:
            InflowBoundaries = list(InflowBoundaries)
            for boundary in InflowBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_OutflowBoundaries.setter
    def _OutflowBoundaries(self, OutflowBoundaries):
        self.__OutflowBoundaries = Boundaries.OutflowBoundaries()
        if OutflowBoundaries is not None:
            OutflowBoundaries = list(OutflowBoundaries)
            for boundary in OutflowBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_WeirBoundaries.setter
    def _WeirBoundaries(self, WeirBoundaries):
        self.__WeirBoundaries = Boundaries.WeirBoundaries()
        if WeirBoundaries is not None:
            WeirBoundaries = list(WeirBoundaries)
            for boundary in WeirBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_CulvertBoundaries.setter
    def _CulvertBoundaries(self, CulvertBoundaries):
        self.__CulvertBoundaries = Boundaries.CulvertBoundaries()
        if CulvertBoundaries is not None:
            CulvertBoundaries = list(CulvertBoundaries)
            for boundary in CulvertBoundaries:
                assert isinstance(boundary, dict)
                self.add_boundary(**boundary)

    @_description.setter
    def _description(self, description):
        self.__description = str(description).strip('\n')

    @_mask_value.setter
    def _mask_value(self, mask_value):
        self.__mask_value = float(mask_value)

class UnstructuredMeshTestCase(unittest.TestCase):

    def setUp(self):
        from AdcircPy.Model.AdcircMesh import AdcircMesh
        self.fort14 = AdcircMesh.parse_fort14(os.getenv('FORT14'))
        self.UnstructuredMesh = UnstructuredMesh()
        self.xy = np.vstack([self.fort14['x'],  self.fort14['y']]).T
        self.SpatialReference = os.getenv('FORT14_EPSG')

    def test_set_xy(self):
        self.UnstructuredMesh.set_xy(self.xy)

    def test_set_node_id(self):
        self.UnstructuredMesh.set_xy(self.xy)
        self.UnstructuredMesh.set_node_id(self.fort14['node_id'])

    def test_set_elements(self):
        self.UnstructuredMesh.set_xy(self.xy)
        self.UnstructuredMesh.set_elements(self.fort14['elements'])

    def test_set_element_id(self):
        self.UnstructuredMesh.set_xy(self.xy)
        self.UnstructuredMesh.set_elements(self.fort14['elements'])
        self.UnstructuredMesh.set_element_id(self.fort14['element_id'])

    def test_set_SpatialReference(self):
        self.UnstructuredMesh.set_SpatialReference(self.SpatialReference)

    def test_set_vertical_datum(self):
        self.UnstructuredMesh.set_vertical_datum('NAVD88')

    def test_add_OceanBoundary(self):
        self.UnstructuredMesh.set_xy(self.xy)
        self.UnstructuredMesh.set_SpatialReference(self.SpatialReference)
        for boundary in self.fort14['OceanBoundaries']:
            self.UnstructuredMesh.add_boundary(None, **boundary)

    def test_set_node_id_raise_empty(self):
        with self.assertRaises(Exception) as context:
            self.UnstructuredMesh.set_node_id(self.fort14['node_id'])
        self.assertTrue('Cannot set node_id on empty mesh.'
                        in str(context.exception))

    def test_set_elements_raise_empty(self):
        with self.assertRaises(Exception) as context:
            self.UnstructuredMesh.set_elements(self.fort14['elements'])
        self.assertTrue('Cannot set elements on empty mesh.'
                        in str(context.exception))

    def test_set_element_id_raise_empty(self):
        with self.assertRaises(Exception) as context:
            self.UnstructuredMesh.set_element_id(self.fort14['element_id'])
        self.assertTrue('Cannot set element_id on mesh with no elements.'
                        in str(context.exception))

    def test_set_SpatialReference_raise_empty(self):
        self.UnstructuredMesh.set_SpatialReference(4326)

    def test_add_OceanBoundary_raise_no_spatial_reference(self):
        self.UnstructuredMesh.set_xy(self.xy)
        for boundary in self.fort14['OceanBoundaries']:
            with self.assertRaises(Exception) as context:
                self.UnstructuredMesh.add_boundary(None, **boundary)
            self.assertTrue('add_boundary() not supported for meshes with no '
                            + 'spatial reference.' in str(context.exception))

    # TODO
    def _test_set_vertical_datum_raise_empty(self):
        self.UnstructuredMesh.set_vertical_datum('NAVD88')

    def _test_add_LandBoundary(self):
        if len(self.fort14['LandBoundaries']) > 0:
            self.UnstructuredMesh.add_boundary(
                                    1, self.fort14['LandBoundaries'].pop())

    def _test_add_InflowBoundary(self):
        if len(self.fort14['InflowBoundaries']) > 0:
            self.UnstructuredMesh.add_boundary(
                                    2, self.fort14['InflowBoundaries'].pop())

    def _test_add_OutflowBoundary(self):
        if len(self.fort14['OutflowBoundaries']) > 0:
            self.UnstructuredMesh.add_boundary(
                                3, self.fort14['OoutflowBoundaries'].pop())

    def _test_add_WeirBoundary(self):
        if len(self.fort14['WeirBoundaries']) > 0:
            self.UnstructuredMesh.add_boundary(
                                    4, self.fort14['WeirBoundaries'].pop())

    def _test_add_CulvertBoundary(self):
        if len(self.fort14['LandBoundaries']) > 0:
            self.UnstructuredMesh.add_boundary(
                                5, self.fort14['CulvertBoundaries'].pop())
