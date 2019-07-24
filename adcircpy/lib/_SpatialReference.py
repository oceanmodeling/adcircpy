import numpy as np
from osgeo import osr


class _SpatialReference(object):

    def __init__(self):
        self.__SpatialReference = None

    @staticmethod
    def transform_vertices(
        vertices,
        input_SpatialReference,
        output_SpatialReference
    ):

        if isinstance(input_SpatialReference, int):
            EPSG = input_SpatialReference
            input_SpatialReference = osr.SpatialReference()
            input_SpatialReference.ImportFromEPSG(EPSG)
        assert isinstance(input_SpatialReference, osr.SpatialReference)

        if isinstance(output_SpatialReference, int):
            EPSG = output_SpatialReference
            output_SpatialReference = osr.SpatialReference()
            output_SpatialReference.ImportFromEPSG(EPSG)
        assert isinstance(output_SpatialReference, osr.SpatialReference)

        if not output_SpatialReference.IsSame(input_SpatialReference):
            CoordinateTransform = osr.CoordinateTransformation(
                                                    input_SpatialReference,
                                                    output_SpatialReference)
            if input_SpatialReference.IsGeographic():
                vertices = [(y, x) for x, y in vertices]
            else:
                vertices = [(x, y) for x, y in vertices]
            vertices = CoordinateTransform.TransformPoints(vertices)
            vertices = np.asarray([(x, y) for x, y, _ in vertices])
        return vertices

    def _get_spatial_reference(self):
        return self.__SpatialReference

    def _set_spatial_reference(self, SpatialReference):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        assert isinstance(SpatialReference, osr.SpatialReference)
        self.__SpatialReference = SpatialReference

    def _clear_spatial_reference(self):
        self.__SpatialReference = None

    # @property
    # def SpatialReference(self):
    #     return self._get_spatial_reference()

    # @property
    # def _SpatialReference(self):
    #     return self.__SpatialReference

    # @SpatialReference.setter
    # def SpatialReference(self, SpatialReference):
    #     if self.__SpatialReference is None:
    #         raise Exception(
    #             'Cannot transform SpatialReference of object: initial '
    #             + 'SpatialReference is not set.')
    #     if SpatialReference is None:
    #         SpatialReference = osr.SpatialReference()
    #     elif isinstance(SpatialReference, int):
    #         EPSG = SpatialReference
    #         SpatialReference = osr.SpatialReference()
    #         SpatialReference.ImportFromEPSG(EPSG)
    #     assert isinstance(SpatialReference, osr.SpatialReference)
    #     self._SpatialReference = SpatialReference

    # @_SpatialReference.setter
    # def _SpatialReference(self, SpatialReference):
    #     raise NotImplementedError(
    #         'Inheriting class must implement _SpatialReference.setter')
        # if SpatialReference is not None:
        #     if isinstance(SpatialReference, int):
        #         EPSG = SpatialReference
        #         SpatialReference = osr.SpatialReference()
        #         SpatialReference.ImportFromEPSG(EPSG)
        #     assert isinstance(SpatialReference, osr.SpatialReference)
        # self.__SpatialReference = SpatialReference
