from typing import Union

import geopandas as gpd
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from shapely.geometry import LineString, MultiLineString

from adcircpy.figures import figure, get_topobathy_kwargs
from adcircpy.mesh.base import Grd  # , sort_edges, signed_polygon_area
from adcircpy.mesh.parsers import grd


class BaseBoundaries:
    def __init__(self, mesh, data):
        self._mesh = mesh
        self._data = data

    @property
    def ids(self):
        if not hasattr(self, '_ids'):
            self._ids = range(len(self._data))
        return self._ids

    @property
    def indexes(self):
        if not hasattr(self, '_indexes'):
            self._indexes = [self._mesh.nodes.index[data['node_id']] for data in self._data]
        return np.array(self._indexes)

    @property
    def node_id(self):
        if not hasattr(self, '_node_id'):
            self._node_id = list()
            for data in self._data:
                self._node_id.append(data['node_id'])
        return self._node_id

    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            data = []
            for i, boundary in enumerate(self._data):
                data.append(
                    {
                        'geometry': LineString(
                            self._mesh.coords.loc[self.indexes[i], :].values
                        ),
                        'key': f'{boundary.get("ibtype")}:{id}',
                        'indexes': self.indexes[i],
                        **boundary,
                    }
                )
            self._gdf = gpd.GeoDataFrame(data, crs=self._mesh.crs)
        return self._gdf

    def __eq__(self, other: 'BaseBoundaries') -> bool:
        return self._data == other._data


class OceanBoundaries(BaseBoundaries):
    pass


class LandBoundaries(BaseBoundaries):
    pass


class InteriorBoundaries(BaseBoundaries):
    pass


class InflowBoundaries(BaseBoundaries):
    pass


class BarrierBaseBoundaries(BaseBoundaries):
    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            data = []
            for i, (id, boundary) in enumerate(self._data.items()):
                front_face, back_face = list(zip(*self.indexes[i]))
                data.append(
                    {
                        'geometry': MultiLineString(
                            [
                                LineString(self._mesh.coords.iloc[front_face, :].values),
                                LineString(self._mesh.coords.iloc[back_face, :].values),
                            ]
                        ),
                        'key': f'{boundary.get("ibtype")}:{id}',
                        **boundary,
                    }
                )
            self._gdf = gpd.GeoDataFrame(data, crs=self._mesh.crs)
        return self._gdf


class OutflowBoundaries(BarrierBaseBoundaries):
    pass


class WeirBoundaries(BarrierBaseBoundaries):
    pass


class CulvertBoundaries(BarrierBaseBoundaries):
    pass


class Fort14Boundaries:
    def __init__(self, fort14: 'Fort14', boundaries: Union[dict, None]):
        self._data = {} if boundaries is None else boundaries
        self._mesh = fort14

    def to_dict(self):
        return self._data

    @figure
    def plot(
        self,
        *args,
        ocean=True,
        land=True,
        interior=True,
        inflow=True,
        outflow=True,
        weir=True,
        culvert=True,
        **kwargs,
    ):
        ax = kwargs['axes']
        if ocean is True:
            self.ocean.gdf.plot(ax=ax)

        if land is True:
            self.land.gdf.plot(ax=ax)

    @property
    def ocean(self):
        if not hasattr(self, '_ocean'):
            self._ocean = OceanBoundaries(self._mesh, self._data.get(None, {}))
        return self._ocean

    @property
    def land(self):
        if not hasattr(self, '_land'):
            self._land = LandBoundaries(self._mesh, self._aggregate_boundaries('0'))
        return self._land

    @property
    def interior(self):
        if not hasattr(self, '_interior'):
            self._interior = InteriorBoundaries(self._mesh, self._aggregate_boundaries('1'))
        return self._interior

    @property
    def inflow(self):
        if not hasattr(self, '_inflow'):
            self._inflow = InflowBoundaries(self._mesh, self._aggregate_boundaries('2'))
        return self._inflow

    @property
    def outflow(self):
        if not hasattr(self, '_outflow'):
            self._outflow = OutflowBoundaries(self._mesh, self._aggregate_boundaries('3'))
        return self._outflow

    @property
    def weir(self):
        if not hasattr(self, '_weir'):
            self._weir = WeirBoundaries(self._mesh, self._aggregate_boundaries('4'),)
        return self._weir

    @property
    def culvert(self):
        if not hasattr(self, '_culvert'):
            self._culvert = CulvertBoundaries(self.mesh, self._aggregate_boundaries('5'))
        return self._culvert

    def _aggregate_boundaries(self, endswith):
        boundaries = {}
        for ibtype, _boundaries in self._data.items():
            if ibtype is None:
                continue
            if ibtype.endswith(endswith):
                for bdata in list(_boundaries.values()):
                    boundaries.update({len(boundaries) + 1: {'ibtype': ibtype, **bdata,}})
        return boundaries

    def __eq__(self, other: 'Fort14Boundaries') -> bool:
        return self._data == other._data


class Fort14(Grd):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """

    # _boundaries = BoundariesDescriptor()

    def __init__(self, *args, boundaries=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._boundaries = Fort14Boundaries(self, boundaries)

    @classmethod
    def open(cls, path, crs=None):
        _grd = grd.read(path, crs=crs)
        _grd['nodes'].iloc[:, 2:] *= -1
        return cls(**_grd)

    def write(self, path, overwrite=False, format='fort.14'):
        if format in ['fort.14']:
            _grd = self.to_dict()
            nodes = self.nodes
            nodes.iloc[:, 2:] *= -1
            _grd['nodes'] = nodes

            grd.write(
                grd=_grd, path=path, overwrite=overwrite,
            )
        else:
            super().write(path=path, overwrite=overwrite, format=format)

    def to_dict(self, boundaries=True):
        _grd = super().to_dict()
        if boundaries is True:
            _grd.update({'nodes': self.nodes, 'boundaries': self.boundaries.to_dict()})
        return _grd

    @figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        # figsize=rcParams["figure.figsize"],
        extent=None,
        cbar_label=None,
        **kwargs,
    ):
        if vmin is None:
            vmin = np.min(self.values.values)
        if vmax is None:
            vmax = np.max(self.values.values)
        kwargs.update(**get_topobathy_kwargs(self.values, vmin, vmax))
        kwargs.pop('col_val')
        levels = kwargs.pop('levels')
        if vmin != vmax:
            self.tricontourf(axes=axes, levels=levels, vmin=vmin, vmax=vmax, **kwargs)
        else:
            self.tripcolor(axes=axes, **kwargs)
        self.quadface(axes=axes, **kwargs)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=kwargs['cmap'])
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes('bottom', size='2%', pad=0.5)
        cbar = plt.colorbar(mappable, cax=cax, orientation='horizontal')
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        return axes

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def ocean_boundaries(self):
        return self.boundaries.ocean

    @property
    def land_boundaries(self):
        return self.boundaries.land

    @property
    def interior_boundaries(self):
        return self.boundaries.interior

    @property
    def inflow_boundaries(self):
        return self.boundaries.inflow

    @property
    def outflow_boundaries(self):
        return self.boundaries.outflow

    @property
    def weir_boundaries(self):
        return self.boundaries.weir

    @property
    def culvert_boundaries(self):
        return self.boundaries.culvert

    def __eq__(self, other: 'Fort14') -> bool:
        return super().__eq__(other) and self.boundaries == other.boundaries
