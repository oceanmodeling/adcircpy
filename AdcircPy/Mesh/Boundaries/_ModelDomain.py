# global imports
from matplotlib.path import Path, get_paths_extents
from osgeo import ogr, osr
from scipy.spatial import cKDTree
import numpy as np

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest
import os
import matplotlib.pyplot as plt


class _ModelDomain(object):

    def __init__(self, SpatialReference, *ExternalBoundaries,
                 InternalBoundaries=None):
        self._SpatialReference = SpatialReference
        self.__init_OuterBoundary(ExternalBoundaries)
        self.__init_InnerBoundaries(InternalBoundaries)
        self.__init_DomainPath()

    def get_boundaries(self, SpatialReference=None, bbox=None, h0=None):
        OuterBoundary = self.get_outer_boundary(SpatialReference, bbox, h0)
        InnerBoundaries = self.get_inner_boundaries(SpatialReference, bbox, h0)
        return OuterBoundary, InnerBoundaries

    def get_outer_boundary(self, SpatialReference=None, bbox=None, h0=None):
        OuterBoundary = self.OuterBoundary
        if h0 is not None:
            OuterBoundary = self.__resample_Path_to_h0(OuterBoundary, h0)
        if bbox is not None:
            OuterBoundary = OuterBoundary.clip_to_bbox(bbox)
        if SpatialReference is not None:
            OuterBoundary = self.__transform_boundary(SpatialReference,
                                                      OuterBoundary)
        return OuterBoundary

    def get_inner_boundaries(self, SpatialReference=None, bbox=None, h0=None):
        InnerBoundaries = self.InnerBoundaries
        for i, boundary in enumerate(InnerBoundaries):
            if h0 is not None:
                boundary = self.__resample_Path_to_h0(boundary, h0)
            if bbox is not None:
                InnerBoundaries[i] = boundary.clip_to_bbox(bbox)
            if SpatialReference is not None:
                InnerBoundaries[i] = self.__transform_boundary(
                                                    SpatialReference, boundary)
        return InnerBoundaries

    def __init_OuterBoundary(self, ExternalBoundaries):
        """
        This function takes the external boundaries and builds a closed outer
        boundary by stitching together the open external boundaries. The outer
        boundary object is necessary for defining the concave hull of the data.
        """
        open_boundaries = list()
        for boundary in ExternalBoundaries:
            if isinstance(boundary, str):
                DataSource = ogr.Open(boundary)
                for Layer in DataSource:
                    for Feature in Layer:
                        FeatureDefn = Feature.GetDefnRef()
                        if FeatureDefn.GetGeomFieldCount() != 1:
                            raise RuntimeError('Feature table must contain '
                                               + 'exactly one geometry column.'
                                               )
                        Geometry = Feature.GetGeometryRef()
                        GeomDefn = FeatureDefn.GetGeomFieldDefn(0)
                        GeometryName = Geometry.GetGeometryName()
                        if 'LINEARRING' in GeometryName:
                            if Layer.GetFeatureCount() > 1:
                                raise TypeError('Cannot use more than one '
                                                + 'Feature when passing a '
                                                + 'LinearRing as Outer '
                                                + 'Boundary.')
                        elif 'LINESTRING' in GeometryName:
                            pass
                        else:
                            raise TypeError('Can only pass LineString as '
                                            + 'Outer Boundary')
                        SpatialReference = GeomDefn.GetSpatialRef()
                        if not SpatialReference.IsSame(self.SpatialReference):
                            Geometry.TransformTo(self.SpatialReference)
                        open_boundaries.append(
                                            np.asarray(Geometry.GetPoints()))

            elif isinstance(boundary, ogr.Geometry):
                SpatialReference = boundary.GetSpatialReference()
                if not SpatialReference.IsSame(self.SpatialReference):
                    _boundary = boundary.Clone()
                    _boundary.TransformTo(self.SpatialReference)
                open_boundaries.append(np.asarray(_boundary.GetPoints()))

            elif isinstance(boundary, _BaseBoundary):
                for boundary_type in boundary.keys():
                    for _boundary in boundary[boundary_type]:
                        Geometry = _boundary['Geometry']
                        SpatialReference = Geometry.GetSpatialReference()
                        if not SpatialReference.IsSame(self.SpatialReference):
                            Geometry.TransformTo(self.SpatialReference)
                        open_boundaries.append(np.asarray(
                                                        Geometry.GetPoints()))
            else:
                raise TypeError('Unknown type for argument ExternalBoundaries')
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
                ordered_vertices = np.vstack([ordered_vertices,
                                             open_boundaries.pop(hidx)])
            else:
                ordered_vertices = np.vstack([
                                    ordered_vertices,
                                    np.flipud(open_boundaries.pop(tidx))])
        self.__OuterBoundary = Path(ordered_vertices, closed=True)

    def __init_InnerBoundaries(self, InternalBoundaries):
        """
        self.InnerBoundaries is a list of Path objects.
        """
        inner_boundaries = list()
        if InternalBoundaries is not None:
            for boundary in InternalBoundaries:
                if isinstance(boundary, str):
                    DataSource = ogr.Open(boundary)
                    for k, Layer in enumerate(DataSource):
                        for Feature in Layer:
                            Geometry = Feature.GetGeometryRef()
                            SpatialRef = Geometry.GetSpatialReference()
                            if not SpatialRef.IsSame(self.SpatialReference):
                                Geometry.TransformTo(self.SpatialReference)
                            GeometryCount = Geometry.GetGeometryCount()
                            if GeometryCount == 1:
                                inner_boundaries.append(
                                            np.asarray(
                                                Geometry.GetPoints())[:, :2])
                            if GeometryCount > 1:
                                _LinearRing = ogr.Geometry(ogr.wkbLinearRing)
                                for i, SubGeom in enumerate(Geometry):
                                    if i == 1:
                                        p = SubGeom.GetPoints()
                                        for x, y in p[::-1]:
                                            _LinearRing.AddPoint(x, y)
                                    else:
                                        for x, y in SubGeom.GetPoints():
                                            _LinearRing.AddPoint(x, y)
                                inner_boundaries.append(
                                        np.asarray(
                                            _LinearRing.GetPoints())[:, :2])

                elif isinstance(boundary, ogr.Geometry):
                    # Should probably restrict boundary types to LinearRing
                    # and LineString.
                    Geometry = boundary.Clone()
                    SpatialReference = boundary.GetSpatialReference()
                    if not SpatialReference.IsSame(self.SpatialReference):
                        Geometry = boundary.Clone()
                        Geometry.TransformTo(self.SpatialReference)
                    inner_boundaries.append(np.asarray(Geometry.GetPoints())
                                            [:, :2])
                    Geometry.Destroy()

                elif isinstance(boundary, _BaseBoundary):
                    for boundary_type in boundary.keys():
                        for _boundary in boundary[boundary_type]:
                            Geometry = _boundary['Geometry']
                            SpatialReference = Geometry.GetSpatialReference()
                            if not SpatialReference.IsSame(
                                                        self.SpatialReference):
                                Geometry.TransformTo(self.SpatialReference)
                            if Geometry.GetGeometryCount() > 1:
                                if Geometry.GetGeometryCount() > 2:
                                    raise TypeError('More than one geometry '
                                                    + 'was passed.')
                                _array = list()
                                for i, subgeom in enumerate(Geometry):
                                    if i == 1:
                                        _array.append(
                                            np.asarray(subgeom.GetPoints()))
                                    else:
                                        _array.append(
                                            np.flipud(
                                                np.asarray(
                                                    subgeom.GetPoints())))

                                _array = np.vstack(_array)
                            else:
                                _array = np.asarray(Geometry.GetPoints())
                            inner_boundaries.append(_array)
                else:
                    raise TypeError('Unknown argument type for '
                                    + 'InternalBoundaries')
        self.__InnerBoundaries = [Path(vertices, closed=True) for vertices in
                                  inner_boundaries]

    def __init_DomainPath(self):
        self.__DomainPath = Path.make_compound_path(*[self.OuterBoundary,
                                                    *self.InnerBoundaries])

    def __transform_boundary(self, SpatialReference, boundary):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        elif isinstance(SpatialReference, osr.SpatialReference):
            pass
        else:
            raise TypeError('Unknown argument type for SpatialReference.')
        if not SpatialReference.IsSame(self.SpatialReference):
            CoordinateTransformation = osr.CoordinateTransformation(
                                                        self.SpatialReference,
                                                        SpatialReference)
            p = CoordinateTransformation.TransformPoints(boundary.vertices)
            p = np.asarray(p)
            return Path(p[:, :2])
        else:
            return boundary

    @staticmethod
    def __resample_Path_to_h0(Path, h0):
        if Path.codes is not None and 79 in Path.codes:
            Geom = ogr.Geometry(ogr.wkbLinearRing)
        else:
            Geom = ogr.Geometry(ogr.wkbLineString)
        for x, y in Path.vertices:
            Geom.AddPoint_2D(x, y)
        Geom.Segmentize(h0)
        return Path(Geom.GetPoints())

    @property
    def bbox(self):
        return get_paths_extents([self.DomainPath])

    @property
    def DomainPath(self):
        return self.__DomainPath

    @property
    def InnerBoundaries(self):
        return self.__InnerBoundaries

    @property
    def OuterBoundary(self):
        return self.__OuterBoundary

    @property
    def SpatialReference(self):
        return self._SpatialReference

    @property
    def _SpatialReference(self):
        return self.__SpatialReference

    @_SpatialReference.setter
    def _SpatialReference(self, SpatialReference):
        if isinstance(SpatialReference, osr.SpatialReference):
            self.__SpatialReference = SpatialReference
        elif isinstance(SpatialReference, int):
            self.__SpatialReference = osr.SpatialReference()
            self.SpatialReference.ImportFromEPSG(SpatialReference)
        else:
            raise TypeError('Unknown type for SpatialReference argument.')


class ModelDomainTestCase(unittest.TestCase):

    def setUp(self):
        self.ExternalBoundaries = [
                        os.getenv('HSOFS_OCEAN_BOUNDARY_SHAPEFILE'),
                        os.getenv('HSOFS_LAND_BOUNDARIES_SHAPEFILE')]
        self.InternalBoundaries = [
                        os.getenv('HSOFS_INNER_BOUNDARIES_SHAPEFILE'),
                        os.getenv('HSOFS_WEIR_BOUNDARIES_SHAPEFILE')]
        self.epsg = 4326

    def test_external_boundaries_shapefile_instantiation(self):
        _ModelDomain(self.epsg, *self.ExternalBoundaries)

    def test_internal_boundaries_shapefile_instantiation(self):
        _ModelDomain(self.epsg, *self.ExternalBoundaries,
                     InternalBoundaries=self.InternalBoundaries)

    def test_get_outer_boundary(self):
        _ModelDomain(self.epsg, *self.ExternalBoundaries).get_outer_boundary()

    def test_get_outer_boundary_epsg(self):
        _ModelDomain(self.epsg, *self.ExternalBoundaries).get_outer_boundary(
                          SpatialReference=3395)

    def test_get_inner_boundaries(self):
        _ModelDomain(self.epsg, *self.ExternalBoundaries,
                     InternalBoundaries=self.InternalBoundaries) \
                     .get_inner_boundaries()

    def test_get_inner_boundaries_epsg(self):
        _ModelDomain(
            self.epsg,
            *self.ExternalBoundaries,
            InternalBoundaries=self.InternalBoundaries) \
           .get_inner_boundaries(SpatialReference=3395)

    def test_get_outer_boundary_from_shapefile_and_plot(self):
        OuterBoundary = _ModelDomain(self.epsg, *self.ExternalBoundaries) \
                       .get_outer_boundary()
        plt.plot(OuterBoundary.vertices[:, 0], OuterBoundary.vertices[:, 1])
        plt.title('shapefile test case')
        plt.gca().axis('scaled')
        plt.show(block=False)
