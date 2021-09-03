from typing import Union

from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.cm import get_cmap
import numpy
from shapely.geometry import MultiPoint, MultiPolygon, Polygon
from shapely.geometry import shape as shapely_shape


def plot_polygon(
    geometry: Union[Polygon, MultiPolygon],
    fill: bool = False,
    axis: Axes = None,
    show: bool = False,
    **kwargs,
) -> Axes:
    """
    Plot the given polygon.

    :param geometry: Shapely polygon (or multipolygon)
    :param axis: `pyplot` axis to plot to
    :param show:  whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    if 'c' not in kwargs:
        try:
            color = next(axis._get_lines.color_cycle)
        except AttributeError:
            color = 'r'
        kwargs['c'] = color

    if isinstance(geometry, dict):
        geometry = shapely_shape(geometry)

    if type(geometry) is Polygon:
        if fill:
            axis.fill(*geometry.exterior.xy, **kwargs)
            kwargs['c'] = 'w'
            for interior in geometry.interiors:
                axis.fill(*interior.xy, **kwargs)
        else:
            axis.plot(*geometry.exterior.xy, **kwargs)
            for interior in geometry.interiors:
                axis.plot(*interior.xy, **kwargs)
    elif type(geometry) is MultiPolygon:
        for polygon in geometry:
            plot_polygon(geometry=polygon, axis=axis, fill=fill, show=False, **kwargs)
    else:
        if fill:
            axis.fill(*geometry.xy, **kwargs)
        else:
            axis.plot(*geometry.xy, **kwargs)

    if show:
        pyplot.show()

    return axis


def plot_polygons(
    geometries: [Polygon],
    colors: [str] = None,
    fill: bool = False,
    axis: Axes = None,
    show: bool = False,
    **kwargs,
) -> Axes:
    """
    Plot the given polygons using the given colors.

    :param geometries: list of shapely polygons or multipolygons
    :param colors: colors to plot each region
    :param axis: `pyplot` axis to plot to
    :param show: whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    if 'c' in kwargs:
        colors = [kwargs['c'] for _ in range(len(geometries))]
    elif colors is None:
        colors = [
            get_cmap('gist_rainbow')(color_index / len(geometries))
            for color_index in range(len(geometries))
        ]

    for geometry_index, geometry in enumerate(geometries):
        kwargs['c'] = colors[geometry_index]
        plot_polygon(geometry=geometry, fill=fill, axis=axis, **kwargs)

    if show:
        pyplot.show()

    return axis


def plot_bounding_box(
    sw: (float, float), ne: (float, float), axis: Axes = None, show: bool = False, **kwargs,
) -> Axes:
    """
    Plot the bounding box of the given extent.

    :param sw: XY coordinates of southwest corner
    :param ne: XY coordinates of northeast corner
    :param axis: `pyplot` axis to plot to
    :param show: whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    corner_points = numpy.array([sw, (ne[0], sw[1]), ne, (sw[0], ne[1]), sw])

    axis.plot(corner_points[:, 0], corner_points[:, 1], **kwargs)

    if show:
        pyplot.show()

    return axis


def plot_points(
    points: Union[numpy.array, MultiPoint],
    index: int = 0,
    axis: Axes = None,
    show: bool = False,
    **kwargs,
) -> Axes:
    """
    Create a scatter plot of the given points.

    :param points: N x M array of points
    :param index: zero-based index of vector layer to read
    :param axis: `pyplot` axis to plot to
    :param show: whether to show the plot
    """

    if type(points) is MultiPoint:
        points = numpy.squeeze(numpy.stack((point._get_coords() for point in points), axis=0))

    if axis is None:
        axis = pyplot.gca()

    if 'c' not in kwargs and points.shape[1] > 2:
        kwargs['c'] = points[:, index + 2]

    if 's' not in kwargs:
        kwargs['s'] = 2

    axis.scatter(points[:, 0], points[:, 1], **kwargs)

    if show:
        pyplot.show()

    return axis
