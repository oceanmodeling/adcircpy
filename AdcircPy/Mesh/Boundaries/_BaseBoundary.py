from collections.abc import Mapping
import abc
import numpy as np
from osgeo import ogr, osr


class _BaseBoundary(Mapping, metaclass=abc.ABCMeta):

    def __init__(self):
        super(_BaseBoundary, self).__init__()

    def __get_LineStringTypeGeometry(self, i, SpatialReference=None):
        vertices = self[i]['vertices']
        _SpatialReference = self[i]['SpatialReference']
        if SpatialReference is not None:
            if _SpatialReference is None:
                raise RuntimeError('Cannot reproject, boundary has no '
                                   + 'projection information.')
            elif not _SpatialReference.IsSame(SpatialReference):
                CoordinateTransform = osr.CoordinateTransformation(
                                                            _SpatialReference,
                                                            SpatialReference)
                vertices = CoordinateTransform.TransformPoints(vertices)
                vertices = np.asarray([(x, y) for x, y, _ in vertices])
        else:
            SpatialReference = _SpatialReference
        Geometry = ogr.Geometry(ogr.wkbLineString)
        if SpatialReference is not None:
            Geometry.AssignSpatialReference(SpatialReference)
        for x, y in vertices:
            Geometry.AddPoint_2D(x, y)
        return Geometry

    def __get_MultiLineStringTypeGeometry(self, i, SpatialReference=None):
        front_face_vertices = self[i]['front_face_vertices']
        back_face_vertices = self[i]['back_face_vertices']
        _SpatialReference = self[i]['SpatialReference']
        if SpatialReference is not None:
            if _SpatialReference is None:
                raise RuntimeError('Cannot reproject, boundary has no '
                                   + 'projection information.')
            elif not _SpatialReference.IsSame(SpatialReference):
                CoordinateTransform = osr.CoordinateTransformation(
                                                            _SpatialReference,
                                                            SpatialReference)
                front_face_vertices = CoordinateTransform.TransformPoints(
                                                        front_face_vertices)
                back_face_vertices = CoordinateTransform.TransformPoints(
                                                            back_face_vertices)
                front_face_vertices = np.asarray(
                                   [(x, y) for x, y, _ in front_face_vertices])
                back_face_vertices = np.asarray(
                                    [(x, y) for x, y, _ in back_face_vertices])
        else:
            SpatialReference = _SpatialReference
        front_face_geom = ogr.Geometry(ogr.wkbLineString)
        back_face_geom = ogr.Geometry(ogr.wkbLineString)
        for x, y in front_face_vertices:
            front_face_geom.AddPoint_2D(x, y)
        for x, y in back_face_vertices:
            back_face_geom.AddPoint_2D(x, y)
        Geometry = ogr.Geometry(ogr.wkbMultiLineString)
        if SpatialReference is not None:
            Geometry.AssignSpatialReference(SpatialReference)
        Geometry.AddGeometry(front_face_geom)
        Geometry.AddGeometry(back_face_geom)
        return Geometry

    def __get_LinearRingTypeGeometry(self, i, SpatialReference=None):
        vertices = self[i]['vertices']
        _SpatialReference = self[i]['SpatialReference']
        if SpatialReference is not None:
            if _SpatialReference is None:
                raise RuntimeError('Cannot reproject, boundary has no '
                                   + 'projection information.')
            elif not _SpatialReference.IsSame(SpatialReference):
                CoordinateTransform = osr.CoordinateTransformation(
                                                            _SpatialReference,
                                                            SpatialReference)
                vertices = CoordinateTransform.TransformPoints(vertices)
                vertices = np.asarray([(x, y) for x, y, _ in vertices])
        else:
            SpatialReference = _SpatialReference
        Geometry = ogr.Geometry(ogr.wkbLinearRing)
        if SpatialReference is not None:
            Geometry.AssignSpatialReference(SpatialReference)
        for x, y in vertices:
            Geometry.AddPoint_2D(x, y)
        return Geometry

    @abc.abstractmethod
    def _add_boundary(self, *args, **kwargs):
        """ """

    @property
    @abc.abstractmethod
    def storage(self):
        """
        Create a class variable __storage and use this property to access it
        within the child object.
        """
        return self.__storage
