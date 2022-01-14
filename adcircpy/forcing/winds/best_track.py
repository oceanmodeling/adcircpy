from datetime import datetime
import io
import os
from os import PathLike
import pathlib
from typing import Union

from matplotlib import pyplot
from matplotlib.axis import Axis
from matplotlib.transforms import Bbox
import numpy as numpy
from pandas import DataFrame
import pandas.util
from pyproj import CRS, Transformer
from shapely import ops
from shapely.geometry import Point, Polygon
from stormevents import VortexTrack
from stormevents.nhc.atcf import ATCF_FileDeck, ATCF_Mode
from stormevents.nhc.track import read_atcf
import utm

from adcircpy.forcing.winds.base import WindForcing
from adcircpy.plotting import plot_coastline, plot_polygon, plot_polygons
from adcircpy.utilities import get_logger

LOGGER = get_logger(__name__)

# suppress `SettingWithCopyWarning` when instantiating `VortexForcing`
# TODO: figure out why that is happening
pandas.options.mode.chained_assignment = None  # default='warn'


class BestTrackForcing(VortexTrack, WindForcing):
    def __init__(
        self,
        storm: Union[str, PathLike, DataFrame, io.BytesIO],
        nws: int = None,
        interval_seconds: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
        mode: ATCF_Mode = None,
        filename: PathLike = None,
    ):
        if nws is None:
            nws = 20

        valid_nws_values = [8, 19, 20]
        assert (
            nws in valid_nws_values
        ), f'ATCF BestTrack can only use `nws` values in {valid_nws_values}'

        if interval_seconds is None:
            interval_seconds = 3600

        VortexTrack.__init__(
            self,
            storm=storm,
            start_date=start_date,
            end_date=end_date,
            file_deck=ATCF_FileDeck.b,
            mode=mode,
            record_type='BEST',
            filename=filename,
        )
        WindForcing.__init__(self, nws=nws, interval_seconds=interval_seconds)

    def summary(
        self, output: Union[str, os.PathLike] = None, overwrite: bool = False,
    ):
        min_storm_speed = numpy.min(self.speed)
        max_storm_speed = numpy.max(self.speed)
        track_length = self.track_length
        duration = self.duration
        min_central_pressure = numpy.min(self.central_pressure)
        max_wind_speed = numpy.max(self.data['max_sustained_wind_speed'])
        start_loc = (self.data['longitude'][0], self.data['latitude'][0])
        end_loc = (self.data['longitude'].iloc[-1], self.data['latitude'].iloc[-1])
        f = [
            f'Summary of storm: {self.storm_id}',
            f'min./max. track speed: {min_storm_speed} m/s, {max_storm_speed} m/s',
            f'min. central pressure: {min_central_pressure} hPa',
            f'max. wind speed: {max_wind_speed} kts',
            f'Starting at: {start_loc} and ended at: {end_loc}',
            f'Total track length: {track_length:.2f} km',
            f'Total track duration: {duration:.2f} days',
        ]
        summary = '\n'.join(f)
        if output is not None:
            if not isinstance(output, pathlib.Path):
                output = pathlib.Path(output)
            if overwrite or not output.exists():
                with open(output, 'w+') as fh:
                    fh.write(summary)
            else:
                LOGGER.warning(f'skipping existing file "{output}"')
        return summary

    @property
    def NWS(self) -> int:
        try:
            return self.__NWS
        except AttributeError:
            return 20

    @NWS.setter
    def NWS(self, NWS: int):
        assert NWS in [8, 19, 20]
        self.__NWS = int(NWS)

    @property
    def BLADj(self) -> float:
        try:
            return self.__BLADj
        except AttributeError:
            return 0.9

    @BLADj.setter
    def BLADj(self, BLADj: float):
        BLADj = float(BLADj)
        assert BLADj >= 0 and BLADj <= 1
        self.__BLADj = BLADj

    @property
    def geofactor(self) -> float:
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @geofactor.setter
    def geofactor(self, geofactor: float):
        geofactor = float(geofactor)
        assert geofactor >= 0 and geofactor <= 1
        self.__geofactor = geofactor

    @classmethod
    def from_fort22(
        cls,
        fort22: PathLike,
        nws: int = None,
        interval_seconds: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
    ) -> 'BestTrackForcing':
        filename = None
        if pathlib.Path(fort22).exists():
            filename = fort22
        return cls(
            storm=read_atcf(fort22),
            start_date=start_date,
            end_date=end_date,
            nws=nws,
            interval_seconds=interval_seconds,
            mode=None,
            filename=filename,
        )

    @classmethod
    def from_atcf_file(
        cls,
        atcf: PathLike,
        nws: int = None,
        interval_seconds: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
    ) -> 'BestTrackForcing':
        filename = None
        if pathlib.Path(atcf).exists():
            filename = atcf
        return cls(
            storm=atcf,
            start_date=start_date,
            end_date=end_date,
            nws=nws,
            interval_seconds=interval_seconds,
            mode=None,
            filename=filename,
        )

    def clip_to_bbox(self, bbox, bbox_crs):
        msg = f'bbox must be a {Bbox} instance.'
        assert isinstance(bbox, Bbox), msg
        bbox_pol = Polygon(
            [
                [bbox.xmin, bbox.ymin],
                [bbox.xmax, bbox.ymin],
                [bbox.xmax, bbox.ymax],
                [bbox.xmin, bbox.ymax],
                [bbox.xmin, bbox.ymin],
            ]
        )
        _switch = True
        unique_dates = numpy.unique(self.dataframe['datetime'])
        _found_start_date = False
        for _datetime in unique_dates:
            records = self.dataframe[self.dataframe['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852.0 * radii  # convert to meters
            lon = records['longitude'].iloc[0]
            lat = records['latitude'].iloc[0]
            _, _, number, letter = utm.from_latlon(lat, lon)
            df_crs = CRS.from_epsg(4326)
            utm_crs = CRS.from_epsg(f'326{number}')
            transformer = Transformer.from_crs(df_crs, utm_crs, always_xy=True)
            p = Point(*transformer.transform(lon, lat))
            pol = p.buffer(radii)
            transformer = Transformer.from_crs(utm_crs, bbox_crs, always_xy=True)
            pol = ops.transform(transformer.transform, pol)
            if _switch is True:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _found_start_date = True
                    _switch = False
                    continue

            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

        if _found_start_date is False:
            raise Exception(f'No data within mesh bounding box for storm {self.storm_id}.')

    def plot_track(
        self,
        axis: Axis = None,
        show: bool = False,
        color: str = 'k',
        coastline: bool = True,
        **kwargs,
    ):
        kwargs.update({'color': color})
        if axis is None:
            fig = pyplot.figure()
            axis = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * numpy.sin(numpy.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * numpy.cos(numpy.deg2rad(self.direction.iloc[i]))
            axis.quiver(self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            if i % 6 == 0:
                axis.annotate(
                    self.data['datetime'].iloc[i],
                    (self.longitude.iloc[i], self.latitude.iloc[i]),
                )
        if show:
            axis.axis('scaled')
        if bool(coastline) is True:
            plot_coastline(axis, show)

    def plot_wind_swath(self, isotach: int, segments: int = 91):
        isotachs = self.isotachs(wind_speed=isotach, segments=segments)
        swath = self.wind_swath(isotach=isotach, segments=segments)

        plot_polygons(isotachs)
        plot_polygon(swath)
        pyplot.suptitle(
            f'{self.storm_id} - isotach {isotach} kt ({self.start_date} - {self.end_date})'
        )
        pyplot.show()
