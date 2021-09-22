from collections.abc import Collection
from datetime import datetime, timedelta
from enum import Enum
import ftplib
from functools import wraps
import gzip
import io
import logging
import os
from os import PathLike
import pathlib
import socket
import time
from typing import Any, TextIO, Union
from urllib.error import URLError
import zipfile

import appdirs
from dateutil.parser import parse as parse_date
import geopandas
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.transforms import Bbox
import numpy as numpy
from pandas import DataFrame, read_csv, Series
import pandas.util
from pyproj import CRS, Geod, Transformer
import requests
from shapely import ops
from shapely.geometry import Point, Polygon
import utm

from adcircpy.forcing.winds.base import WindForcing
from adcircpy.plotting import plot_polygon, plot_polygons

logger = logging.getLogger(__name__)

# suppress `SettingWithCopyWarning` when instantiating `VortexForcing`
# TODO: figure out why that is happening
pandas.options.mode.chained_assignment = None  # default='warn'


def download_coastline(overwrite: bool = False) -> pathlib.Path:
    data_directory = pathlib.Path(appdirs.user_data_dir('ne_coastline'))
    if not data_directory.exists():
        data_directory.mkdir(exist_ok=True, parents=True)

    coastline_filename = data_directory / 'ne_110m_coastline.shp'

    if not coastline_filename.exists() or overwrite:
        # download and save if not present
        url = 'http://naciscdn.org/naturalearth/110m/physical/ne_110m_coastline.zip'
        response = requests.get(url, stream=True)
        with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
            for member_filename in zip_file.namelist():
                file_data = zip_file.read(member_filename)
                with open(data_directory / member_filename, 'wb') as output_file:
                    output_file.write(file_data)
        assert coastline_filename.exists(), 'coastline file not downloaded'

    return coastline_filename


def plot_coastline(axis: Axes = None, show: bool = False, save_filename: PathLike = None):
    if axis is None:
        figure = pyplot.figure()
        axis = figure.add_subplot(1, 1, 1)

    coastline_filename = download_coastline()
    dataframe = geopandas.read_file(coastline_filename)
    dataframe.plot(ax=axis)

    if save_filename is not None:
        pyplot.savefig(save_filename)

    if show:
        pyplot.show()


class FileDeck(Enum):
    """http://hurricanes.ral.ucar.edu/realtime/"""

    a = 'a'
    b = 'b'
    c = 'c'
    d = 'd'
    e = 'e'
    f = 'f'


class Mode(Enum):
    historical = 'ARCHIVE'
    realtime = 'aid_public'


class VortexForcing:
    def __init__(
        self,
        storm: Union[str, PathLike, DataFrame, io.BytesIO],
        start_date: datetime = None,
        end_date: datetime = None,
        file_deck: FileDeck = None,
        mode: Mode = None,
        record_type: str = None,
        filename: PathLike = None,
    ):
        self.__dataframe = None
        self.__filename = None

        self.__atcf = None
        self.__storm_id = None
        self.__start_date = start_date  # initially used to filter A-deck here
        self.__end_date = None
        self.__file_deck = None
        self.__mode = None
        self.__record_type = None

        self.__invalid_storm_name = False
        self.__location_hash = None

        self.file_deck = file_deck
        self.mode = mode
        self.record_type = record_type
        self.filename = filename

        if isinstance(storm, DataFrame):
            self.dataframe = storm
        elif isinstance(storm, io.BytesIO):
            self.__atcf = storm
        elif isinstance(storm, (str, PathLike, pathlib.Path)):
            if os.path.exists(storm):
                self.__atcf = io.open(storm, 'rb')
            else:
                try:
                    self.storm_id = storm
                except ValueError:
                    if pathlib.Path(storm).exists():
                        self.filename = storm
                    else:
                        raise

        self.__previous_configuration = {
            'storm_id': self.storm_id,
            'file_deck': self.file_deck,
            'mode': self.mode,
            'filename': self.filename,
        }

        # use start and end dates to mask dataframe here
        self.start_date = start_date
        self.end_date = end_date

    @property
    def filename(self) -> pathlib.Path:
        return self.__filename

    @filename.setter
    def filename(self, filename: PathLike):
        if filename is not None and not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)
        self.__filename = filename

    def __str__(self):
        record_number = self.__generate_record_numbers()
        lines = []

        dataframe = self.data

        for column in dataframe.select_dtypes(include=['float']):
            if column not in ['latitude', 'longitude']:
                dataframe.loc[:, column] = (
                    dataframe[column].round(0).astype('Int64', copy=False)
                )

        for i, (_, row) in enumerate(dataframe.iterrows()):
            line = []

            line.extend(
                [
                    f'{row["basin"]:<2}',
                    f'{row["storm_number"]:>3}',
                    f'{row["datetime"]:%Y%m%d%H}'.rjust(11),
                    f'{"":3}',
                    f'{row["record_type"]:>5}',
                    f'{convert_value((row["datetime"] - self.start_date) / timedelta(hours=1), to_type=int):>4}',
                ]
            )

            latitude = convert_value(row['latitude'] / 0.1, to_type=int, round_digits=1)
            if latitude >= 0:
                line.append(f'{latitude:>4}N')
            else:
                line.append(f'{latitude * -.1:>4}S')

            longitude = convert_value(row['longitude'] / 0.1, to_type=int, round_digits=1)
            if longitude >= 0:
                line.append(f'{longitude:>5}E')
            else:
                line.append(f'{longitude * -1:>5}W')

            line.extend(
                [
                    f'{convert_value(row["max_sustained_wind_speed"], to_type=int, round_digits=0):>4}',
                    f'{convert_value(row["central_pressure"], to_type=int, round_digits=0):>5}',
                    f'{row["development_level"]:>3}',
                    f'{convert_value(row["isotach"], to_type=int, round_digits=0):>4}',
                    f'{row["quadrant"]:>4}',
                    f'{convert_value(row["radius_for_NEQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_SEQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_SWQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_NWQ"], to_type=int, round_digits=0):>5}',
                ]
            )

            if row['background_pressure'] is None:
                row['background_pressure'] = self.data['background_pressure'].iloc[i - 1]
            if (
                row['background_pressure'] <= row['central_pressure']
                and 1013 > row['central_pressure']
            ):
                background_pressure = 1013
            elif (
                row['background_pressure'] <= row['central_pressure']
                and 1013 <= row['central_pressure']
            ):
                background_pressure = convert_value(
                    row['central_pressure'] + 1, to_type=int, round_digits=0,
                )
            else:
                background_pressure = convert_value(
                    row['background_pressure'], to_type=int, round_digits=0,
                )
            line.append(f'{background_pressure:>5}')

            line.extend(
                [
                    f'{convert_value(row["radius_of_last_closed_isobar"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_of_maximum_winds"], to_type=int, round_digits=0):>4}',
                    f'{"":>5}',  # gust
                    f'{"":>4}',  # eye
                    f'{"":>4}',  # subregion
                    f'{"":>4}',  # maxseas
                    f'{"":>4}',  # initials
                    f'{row["direction"]:>3}',
                    f'{row["speed"]:>4}',
                    f'{row["name"]:^12}',
                ]
            )

            # from this point forwards it's all aswip
            line.append(f'{record_number[i]:>4}')

            lines.append(','.join(line))

        return '\n'.join(lines)

    def write(self, path: PathLike, overwrite: bool = False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        if overwrite or not path.exists():
            with open(path, 'w') as f:
                f.write(str(self))
        else:
            logger.warning(f'skipping existing file "{path}"')

    @property
    def storm_id(self) -> str:
        if self.__storm_id is None and not self.__invalid_storm_name:
            if self.__dataframe is not None:
                storm_id = (
                    f'{self.__dataframe["basin"].iloc[-1]}'
                    f'{self.__dataframe["storm_number"].iloc[-1]}'
                    f'{self.__dataframe["datetime"].iloc[-1].year}'
                )
                try:
                    self.storm_id = storm_id
                except ValueError:
                    try:
                        storm_id = get_atcf_id(
                            storm_name=self.__dataframe['name'].tolist()[-1],
                            year=self.__dataframe['datetime'].tolist()[-1].year,
                        )
                        self.storm_id = storm_id
                    except ValueError:
                        self.__invalid_storm_name = True
        return self.__storm_id

    @storm_id.setter
    def storm_id(self, storm_id: str):
        if storm_id is not None:
            digits = sum([1 for character in storm_id if character.isdigit()])

            if digits == 4:
                atcf_id = get_atcf_id(storm_name=storm_id[:-4], year=int(storm_id[-4:]))
                if atcf_id is None:
                    raise ValueError(f'No storm with id: {storm_id}')
                storm_id = atcf_id
        self.__storm_id = storm_id

    @property
    def track_length(self) -> float:
        geodetic = Geod(ellps='WGS84')

        forward_azimuths, inverse_azimuths, distances = geodetic.inv(
            self.data['longitude'].iloc[:-1],
            self.data['latitude'].iloc[:-1],
            self.data['longitude'].iloc[1:],
            self.data['latitude'].iloc[1:],
        )

        return numpy.sum(distances)

    @property
    def start_date(self) -> datetime:
        return self.__start_date

    @start_date.setter
    def start_date(self, start_date: datetime):
        if start_date is None:
            start_date = self.dataframe['datetime'][0]
        else:
            if not isinstance(start_date, datetime):
                start_date = parse_date(start_date)
            if (
                start_date < self.dataframe['datetime'].iloc[0]
                or start_date > self.dataframe['datetime'].iloc[-1]
            ):
                raise ValueError(
                    f'given start date is outside of data bounds ({self.dataframe["datetime"].iloc[0]} - {self.dataframe["datetime"].iloc[-1]})'
                )
        self.__start_date = start_date

    @property
    def end_date(self) -> datetime:
        return self.__end_date

    @end_date.setter
    def end_date(self, end_date: datetime):
        if end_date is None:
            end_date = self.dataframe['datetime'].iloc[-1]
        else:
            if not isinstance(end_date, datetime):
                end_date = parse_date(end_date)
            if (
                end_date < self.dataframe['datetime'].iloc[0]
                or end_date > self.dataframe['datetime'].iloc[-1]
            ):
                raise ValueError(
                    f'given end date is outside of data bounds ({self.dataframe["datetime"].iloc[0]} - {self.dataframe["datetime"].iloc[-1]})'
                )
            if end_date <= self.start_date:
                raise ValueError(f'end date must be after start date ({self.start_date})')
        self.__end_date = end_date

    @property
    def duration(self) -> float:
        d = (self.end_date - self.start_date) / timedelta(days=1)
        return d

    @property
    def file_deck(self) -> FileDeck:
        return self.__file_deck

    @property
    def name(self) -> str:
        return self.data['name'].value_counts()[:].index.tolist()[0]

    @property
    def basin(self) -> str:
        return self.data['basin'].iloc[0]

    @property
    def storm_number(self) -> str:
        return self.data['storm_number'].iloc[0]

    @property
    def year(self) -> int:
        return self.data['datetime'].iloc[0].year

    @property
    def datetime(self):
        return self.data['datetime']

    @property
    def central_pressure(self):
        return self.data['central_pressure']

    @property
    def speed(self):
        return self.data['speed']

    @property
    def direction(self):
        return self.data['direction']

    @property
    def longitude(self):
        return self.data['longitude']

    @property
    def latitude(self):
        return self.data['latitude']

    @file_deck.setter
    def file_deck(self, file_deck: FileDeck):
        if file_deck is None:
            file_deck = FileDeck.a
        elif not isinstance(file_deck, FileDeck):
            file_deck = convert_value(file_deck, FileDeck)
        self.__file_deck = file_deck

    @property
    def mode(self) -> Mode:
        return self.__mode

    @mode.setter
    def mode(self, mode: Mode):
        if mode is None:
            mode = Mode.historical
        elif not isinstance(mode, Mode):
            mode = convert_value(mode, Mode)
        self.__mode = mode

    @property
    def record_type(self) -> str:
        return self.__record_type

    @record_type.setter
    def record_type(self, record_type: str):
        # e.g. BEST, OFCL, HWRF, etc.
        if record_type is not None:
            if self.file_deck == FileDeck.a:
                # see ftp://ftp.nhc.noaa.gov/atcf/docs/nhc_techlist.dat
                # there are more but they may not have enough columns
                record_types_list = ['OFCL', 'OFCP', 'HWRF', 'HMON', 'CARQ']
            elif self.file_deck == FileDeck.b:
                record_types_list = ['BEST']
            else:
                raise NotImplementedError(f'file deck {self.file_deck.value} not implemented')
            if record_type not in record_types_list:
                raise ValueError(
                    f'request_record_type = {record_type} not allowed, select from {record_types_list}'
                )
        self.__record_type = record_type

    @property
    def data(self):
        start_date_mask = self.dataframe['datetime'] >= self.start_date
        if self.end_date is None:
            return self.dataframe[start_date_mask]
        else:
            return self.dataframe[
                start_date_mask & (self.dataframe['datetime'] <= self.__file_end_date)
            ]

    @property
    def atcf(self) -> open:
        configuration = {
            'storm_id': self.storm_id,
            'file_deck': self.file_deck,
            'mode': self.mode,
            'filename': self.filename,
        }

        if (
            self.storm_id is not None
            and self.__atcf is None
            or configuration != self.__previous_configuration
        ):
            self.__atcf = get_atcf_file(self.storm_id, self.file_deck, self.mode)

        return self.__atcf

    @property
    def dataframe(self):
        configuration = {
            'storm_id': self.storm_id,
            'file_deck': self.file_deck,
            'mode': self.mode,
            'filename': self.filename,
        }

        # only download new file if the configuration has changed since the last download
        if (
            self.__dataframe is None
            or len(self.__dataframe) == 0
            or configuration != self.__previous_configuration
        ):
            if configuration['filename'] is not None:
                dataframe = read_atcf(configuration['filename'])
            else:
                # https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt

                columns = [
                    'basin',
                    'storm_number',
                    'datetime',
                    'record_type',
                    'latitude',
                    'longitude',
                    'max_sustained_wind_speed',
                    'central_pressure',
                    'development_level',
                    'isotach',
                    'quadrant',
                    'radius_for_NEQ',
                    'radius_for_SEQ',
                    'radius_for_SWQ',
                    'radius_for_NWQ',
                    'background_pressure',
                    'radius_of_last_closed_isobar',
                    'radius_of_maximum_winds',
                    'name',
                    'direction',
                    'speed',
                ]

                atcf = self.atcf
                if isinstance(atcf, io.BytesIO):
                    # test if Gzip file
                    atcf.seek(0)  # rewind
                    if atcf.read(2) == b'\x1f\x8b':
                        atcf.seek(0)  # rewind
                        atcf = gzip.GzipFile(fileobj=atcf)
                    else:
                        atcf.seek(0)  # rewind

                start_date = self.start_date
                # Only accept request record type or
                # BEST track or OFCL (official) advisory by default
                allowed_record_types = self.record_type
                if allowed_record_types is None:
                    allowed_record_types = ['BEST', 'OFCL']
                records = []

                for line_index, line in enumerate(atcf):
                    line = line.decode('UTF-8').split(',')

                    record = {
                        'basin': line[0],
                        'storm_number': line[1].strip(' '),
                    }

                    record['record_type'] = line[4].strip(' ')

                    if record['record_type'] not in allowed_record_types:
                        continue

                    # computing the actual datetime based on record_type
                    if record['record_type'] == 'BEST':
                        # Add minutes line to base datetime
                        minutes = line[3].strip(' ')
                        if minutes == "":
                            minutes = '00'
                        record['datetime'] = parse_date(line[2].strip(' ') + minutes)
                    else:
                        # Add validation time to base datetime
                        minutes = '00'
                        record['datetime'] = parse_date(line[2].strip(' ') + minutes)
                        if start_date is not None:
                            # Only keep records where base date == start time for advisories
                            if start_date != record['datetime']:
                                continue
                        validation_time = int(line[5].strip(' '))
                        record['datetime'] = record['datetime'] + timedelta(
                            hours=validation_time
                        )

                    latitude = line[6]
                    if 'N' in latitude:
                        latitude = float(latitude.strip('N '))
                    elif 'S' in latitude:
                        latitude = float(latitude.strip('S ')) * -1
                    latitude *= 0.1
                    record['latitude'] = latitude

                    longitude = line[7]
                    if 'E' in longitude:
                        longitude = float(longitude.strip('E ')) * 0.1
                    elif 'W' in longitude:
                        longitude = float(longitude.strip('W ')) * -0.1
                    record['longitude'] = longitude

                    record.update(
                        {
                            'max_sustained_wind_speed': float(line[8].strip(' ')),
                            'central_pressure': float(line[9].strip(' ')),
                            'development_level': line[10].strip(' '),
                        }
                    )

                    try:
                        record['isotach'] = int(line[11].strip(' '))
                    except ValueError:
                        raise Exception(
                            'Error: No radial wind information for this storm; '
                            'parametric wind model cannot be built.'
                        )

                    record.update(
                        {
                            'quadrant': line[12].strip(' '),
                            'radius_for_NEQ': int(line[13].strip(' ')),
                            'radius_for_SEQ': int(line[14].strip(' ')),
                            'radius_for_SWQ': int(line[15].strip(' ')),
                            'radius_for_NWQ': int(line[16].strip(' ')),
                        }
                    )

                    if len(line) > 18:
                        record.update(
                            {
                                'background_pressure': int(line[17].strip(' ')),
                                'radius_of_last_closed_isobar': int(line[18].strip(' ')),
                                'radius_of_maximum_winds': int(line[19].strip(' ')),
                            }
                        )

                        if len(line) > 23:
                            record['name'] = line[27].strip(' ')
                        else:
                            record['name'] = ""
                    else:
                        record.update(
                            {
                                'background_pressure': records[-1]['background_pressure'],
                                'radius_of_last_closed_isobar': records[-1][
                                    'radius_of_last_closed_isobar'
                                ],
                                'radius_of_maximum_winds': records[-1][
                                    'radius_of_maximum_winds'
                                ],
                                'name': records[-1]['name'],
                            }
                        )

                    records.append(record)

                if len(records) == 0:
                    raise ValueError(f'no records found with type(s) "{allowed_record_types}"')

                dataframe = DataFrame.from_records(data=records, columns=columns)

            self.__dataframe = dataframe
            self.__previous_configuration = configuration

        # if location values have changed, recompute velocity
        location_hash = pandas.util.hash_pandas_object(
            self.__dataframe[['longitude', 'latitude']]
        ).sum()
        if self.__location_hash is None or location_hash != self.__location_hash:
            self.__dataframe = self.__compute_velocity(self.__dataframe)
            self.__location_hash = location_hash

        return self.__dataframe

    @dataframe.setter
    def dataframe(self, dataframe: DataFrame):
        self.__dataframe = dataframe

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
        axis: Axes = None,
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

    def isotachs(self, wind_speed: float, segments: int = 91) -> [Polygon]:
        """
        Get the isotach at the given speed at every time in the dataset.

        :param wind_speed: wind speed to extract (in knots)
        :param segments: number of discretization points per quadrant
        :return: list of isotachs as polygons
        """

        ## Collect the attributes needed from the forcing to generate swath
        data = self.data[self.data['isotach'] == wind_speed]

        # convert quadrant radii from nautical miles to meters
        quadrants = ['radius_for_NEQ', 'radius_for_NWQ', 'radius_for_SWQ', 'radius_for_SEQ']
        data[quadrants] *= 1852.0  # nautical miles to meters

        geodetic = Geod(ellps='WGS84')

        ## Generate overall swath based on the desired isotach
        polygons = []
        for index, row in data.iterrows():
            # get the starting angle range for NEQ based on storm direction
            rot_angle = 360 - row['direction']
            start_angle = 0 + rot_angle
            end_angle = 90 + rot_angle

            # append quadrants in counter-clockwise direction from NEQ
            arcs = []
            for quadrant in quadrants:
                # enter the angle range for this quadrant
                theta = numpy.linspace(start_angle, end_angle, segments)
                # move angle to next quadrant
                start_angle = start_angle + 90
                end_angle = end_angle + 90
                # skip if quadrant radius is zero
                if row[quadrant] <= 1.0:
                    continue
                # make the coordinate list for this quadrant
                ## entering origin
                coords = [row[['longitude', 'latitude']].tolist()]
                # using forward geodetic (origin,angle,dist)
                for az12 in theta:
                    lont, latt, backaz = geodetic.fwd(
                        lons=row['longitude'],
                        lats=row['latitude'],
                        az=az12,
                        dist=row[quadrant],
                    )
                    coords.append((lont, latt))
                ## start point equals last point
                coords.append(coords[0])
                # enter quadrant as new polygon
                arcs.append(Polygon(coords))
            polygons.append(ops.unary_union(arcs))
        return polygons

    def wind_swath(
        self, isotach: int, segments: int = 91, plot_swath: bool = False,
    ) -> Polygon:
        """
        extracts the wind swath of the BestTrackForcing class object as a Polygon object

        :param isotach: the wind swath to extract (34-kt, 50-kt, or 64-kt)
        :param segments: number of discretization points per quadrant (default = 91)
        :param plot_swath: plot the swath before returning? (default = False)
        """

        # parameter
        # isotach should be one of 34, 50, 64
        valid_isotach_values = [34, 50, 64]
        assert (
            isotach in valid_isotach_values
        ), f'`isotach` value in `get_wind_swath` must be one of {valid_isotach_values}'

        isotachs = self.isotachs(wind_speed=isotach, segments=segments)

        convex_hulls = []
        for index in range(len(isotachs) - 1):
            convex_hulls.append(
                ops.unary_union([isotachs[index], isotachs[index + 1]]).convex_hull
            )

        # get the union of polygons
        swath = ops.unary_union(convex_hulls)

        if plot_swath:
            plot_polygons(isotachs)
            plot_polygon(swath)
            pyplot.suptitle(
                f'{self.storm_id} - isotach {isotach} kt ({self.start_date} - {self.end_date})'
            )
            pyplot.show()

        return swath

    def __generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i - 1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number

    @property
    def __file_end_date(self):
        unique_dates = numpy.unique(self.dataframe['datetime'])
        for date in unique_dates:
            if date >= numpy.datetime64(self.end_date):
                return date

    def __copy__(self) -> 'VortexForcing':
        return self.__class__(
            storm=self.dataframe.copy(),
            start_date=self.start_date,
            end_date=self.end_date,
            file_deck=self.file_deck,
            record_type=self.record_type,
        )

    def __eq__(self, other: 'VortexForcing') -> bool:
        return numpy.all(self.dataframe == other.dataframe)

    @staticmethod
    def __compute_velocity(data: DataFrame) -> DataFrame:
        """Output has units of meters per second."""

        geodetic = Geod(ellps='WGS84')

        unique_datetimes = numpy.unique(data['datetime'])
        for datetime_index, unique_datetime in enumerate(unique_datetimes):
            unique_datetime_indices = numpy.where(
                numpy.asarray(data['datetime']) == unique_datetime
            )[0]
            for unique_datetime_index in unique_datetime_indices:
                if unique_datetime_indices[-1] + 1 < len(data['datetime']):
                    dt = (
                        data['datetime'].iloc[unique_datetime_indices[-1] + 1]
                        - data['datetime'].iloc[unique_datetime_index]
                    )
                    forward_azimuth, inverse_azimuth, distance = geodetic.inv(
                        data['longitude'].iloc[unique_datetime_indices[-1] + 1],
                        data['latitude'].iloc[unique_datetime_indices[-1] + 1],
                        data['longitude'].iloc[unique_datetime_index],
                        data['latitude'].iloc[unique_datetime_index],
                    )
                else:
                    dt = (
                        data['datetime'].iloc[unique_datetime_index]
                        - data['datetime'].iloc[unique_datetime_indices[0] - 1]
                    )
                    forward_azimuth, inverse_azimuth, distance = geodetic.inv(
                        data['longitude'].iloc[unique_datetime_indices[0] - 1],
                        data['latitude'].iloc[unique_datetime_indices[0] - 1],
                        data['longitude'].iloc[unique_datetime_index],
                        data['latitude'].iloc[unique_datetime_index],
                    )

                speed = distance / (dt / timedelta(seconds=1))
                bearing = inverse_azimuth % 360

                data['speed'].iloc[unique_datetime_index] = speed
                data['direction'].iloc[unique_datetime_index] = bearing

                data['speed'] = data['speed'].astype('float', copy=False)
                data['direction'] = data['direction'].astype('float', copy=False)

        return data

    @classmethod
    def from_fort22(
        cls, fort22: PathLike, start_date: datetime = None, end_date: datetime = None,
    ) -> 'VortexForcing':
        filename = None
        if pathlib.Path(fort22).exists():
            filename = fort22
        return cls(
            storm=read_atcf(fort22),
            start_date=start_date,
            end_date=end_date,
            file_deck=None,
            mode=None,
            record_type=None,
            filename=filename,
        )

    @classmethod
    def from_atcf_file(
        cls, atcf: PathLike, start_date: datetime = None, end_date: datetime = None,
    ) -> 'VortexForcing':
        filename = None
        if pathlib.Path(atcf).exists():
            filename = atcf
        return cls(
            storm=atcf,
            start_date=start_date,
            end_date=end_date,
            file_deck=None,
            mode=None,
            record_type=None,
            filename=filename,
        )


class BestTrackForcing(VortexForcing, WindForcing):
    def __init__(
        self,
        storm: Union[str, PathLike, DataFrame, io.BytesIO],
        nws: int = None,
        interval_seconds: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
        mode: Mode = None,
        filename: PathLike = None,
    ) -> 'BestTrackForcing':
        if nws is None:
            nws = 20

        valid_nws_values = [8, 19, 20]
        assert (
            nws in valid_nws_values
        ), f'ATCF BestTrack can only use `nws` values in {valid_nws_values}'

        if interval_seconds is None:
            interval_seconds = 3600

        VortexForcing.__init__(
            self,
            storm=storm,
            start_date=start_date,
            end_date=end_date,
            file_deck=FileDeck.b,
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
                logger.warning(f'skipping existing file "{output}"')
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


def convert_value(value: Any, to_type: type, round_digits: int = None) -> Any:
    if type(value).__name__ == 'Quantity':
        value = value.magnitude
    if issubclass(to_type, Enum):
        try:
            value = to_type[value]
        except (KeyError, ValueError):
            try:
                value = to_type(value)
            except (KeyError, ValueError):
                raise ValueError(
                    f'unrecognized entry "{value}"; must be one of {list(to_type)}'
                )
    elif value is not None and value != "":
        if round_digits is not None and issubclass(to_type, (int, float)):
            if isinstance(value, str):
                value = float(value)
            value = round(value, round_digits)
        value = to_type(value)
    return value


def retry(ExceptionToCheck, tries=4, delay=3, backoff=2, logger=None):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param ExceptionToCheck: the exception to check. may be a tuple of
        exceptions to check
    :type ExceptionToCheck: Exception or tuple
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """

    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    msg = '%s, Retrying in %d seconds...' % (str(e), mdelay)
                    if logger:
                        logger.warning(msg)
                    # else:
                    #     print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry


def get_atcf_entry(
    year: int, basin: str = None, storm_number: int = None, storm_name: str = None,
) -> Series:
    url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/storm.table'

    try:
        storm_table = read_csv(url, header=None)
    except URLError:
        raise ConnectionError(f'cannot connect to "{url}"')

    if basin is not None and storm_number is not None:
        row = storm_table[
            (storm_table[1] == f'{basin.upper():>3}')
            & (storm_table[7] == storm_number)
            & (storm_table[8] == int(year))
        ]
    elif storm_name is not None:
        row = storm_table[
            (storm_table[0] == f'{storm_name.upper():>10}') & (storm_table[8] == int(year))
        ]
    else:
        raise ValueError('need either storm name or basin + storm number')

    if len(row) > 0:
        return list(row.iterrows())[0][1]


def get_atcf_id(storm_name: str, year: int) -> str:
    return get_atcf_entry(storm_name=storm_name, year=year)[20].strip()


def get_atcf_file(storm_id: str, file_deck: FileDeck = None, mode: Mode = None) -> io.BytesIO:
    url = atcf_url(file_deck=file_deck, storm_id=storm_id, mode=mode).replace('ftp://', "")
    logger.info(f'Downloading storm data from {url}')

    hostname, filename = url.split('/', 1)

    handle = io.BytesIO()

    try:
        ftp = ftplib.FTP(hostname, 'anonymous', "")
        ftp.encoding = 'utf-8'
        ftp.retrbinary(f'RETR {filename}', handle.write)
    except socket.gaierror:
        raise ConnectionError(f'cannot connect to {hostname}')

    return handle


def read_atcf(track: PathLike) -> DataFrame:
    try:
        if not isinstance(track, TextIO):
            track = open(track)
        track = track.readlines()
    except FileNotFoundError:
        # check if the entire track file was passed as a string
        track = str(track).splitlines()
        if len(track) == 1:
            raise

    data = {
        'basin': [],
        'storm_number': [],
        'datetime': [],
        'record_type': [],
        'latitude': [],
        'longitude': [],
        'max_sustained_wind_speed': [],
        'central_pressure': [],
        'development_level': [],
        'isotach': [],
        'quadrant': [],
        'radius_for_NEQ': [],
        'radius_for_SEQ': [],
        'radius_for_SWQ': [],
        'radius_for_NWQ': [],
        'background_pressure': [],
        'radius_of_last_closed_isobar': [],
        'radius_of_maximum_winds': [],
        'name': [],
        'direction': [],
        'speed': [],
    }

    for index, row in enumerate(track):
        row = [value.strip() for value in row.split(',')]

        row_data = {key: None for key in data}

        row_data['basin'] = row[0]
        row_data['storm_number'] = row[1]
        row_data['datetime'] = datetime.strptime(row[2], '%Y%m%d%H')
        row_data['record_type'] = row[4]

        latitude = row[6]
        if 'N' in latitude:
            latitude = float(latitude[:-1]) * 0.1
        elif 'S' in latitude:
            latitude = float(latitude[:-1]) * -0.1
        row_data['latitude'] = latitude

        longitude = row[7]
        if 'E' in longitude:
            longitude = float(longitude[:-1]) * 0.1
        elif 'W' in longitude:
            longitude = float(longitude[:-1]) * -0.1
        row_data['longitude'] = longitude

        row_data['max_sustained_wind_speed'] = convert_value(
            row[8], to_type=int, round_digits=0,
        )
        row_data['central_pressure'] = convert_value(row[9], to_type=int, round_digits=0)
        row_data['development_level'] = row[10]
        row_data['isotach'] = convert_value(row[11], to_type=int, round_digits=0)
        row_data['quadrant'] = row[12]
        row_data['radius_for_NEQ'] = convert_value(row[13], to_type=int, round_digits=0)
        row_data['radius_for_SEQ'] = convert_value(row[14], to_type=int, round_digits=0)
        row_data['radius_for_SWQ'] = convert_value(row[15], to_type=int, round_digits=0)
        row_data['radius_for_NWQ'] = convert_value(row[16], to_type=int, round_digits=0)
        row_data['background_pressure'] = convert_value(row[17], to_type=int, round_digits=0)
        row_data['radius_of_last_closed_isobar'] = convert_value(
            row[18], to_type=int, round_digits=0,
        )
        row_data['radius_of_maximum_winds'] = convert_value(
            row[19], to_type=int, round_digits=0,
        )
        row_data['direction'] = convert_value(row[25], to_type=int)
        row_data['speed'] = convert_value(row[26], to_type=int)
        row_data['name'] = row[27]

        for key, value in row_data.items():
            if isinstance(data[key], Collection):
                data[key].append(value)
            elif data[key] is None:
                data[key] = value

    return DataFrame(data=data)


def atcf_url(file_deck: FileDeck = None, storm_id: str = None, mode: Mode = None,) -> str:
    if storm_id is not None:
        year = int(storm_id[4:])
    else:
        year = None

    if mode is None:
        entry = get_atcf_entry(basin=storm_id[:2], storm_number=int(storm_id[2:4]), year=year)
        mode = entry[18].strip()

    if file_deck is not None and not isinstance(file_deck, FileDeck):
        try:
            file_deck = convert_value(file_deck, FileDeck)
        except ValueError:
            file_deck = None
    if file_deck is None:
        file_deck = FileDeck.a

    if mode is not None and not isinstance(mode, Mode):
        try:
            mode = convert_value(mode, Mode)
        except ValueError:
            mode = None
    if mode is None:
        mode = Mode.realtime

    if mode == Mode.historical:
        nhc_dir = f'archive/{year}'
        suffix = '.dat.gz'
    elif mode == Mode.realtime:
        if file_deck == FileDeck.a:
            nhc_dir = 'aid_public'
            suffix = '.dat.gz'
        elif file_deck == FileDeck.b:
            nhc_dir = 'btk'
            suffix = '.dat'

    url = f'ftp://ftp.nhc.noaa.gov/atcf/{nhc_dir}/'

    if storm_id is not None:
        url += f'{file_deck.value}{storm_id.lower()}{suffix}'

    return url


def atcf_storm_ids(file_deck: FileDeck = None, mode: Mode = None) -> [str]:
    if file_deck is None:
        file_deck = FileDeck.a
    elif not isinstance(file_deck, FileDeck):
        file_deck = convert_value(file_deck, FileDeck)

    url = atcf_url(file_deck=file_deck, mode=mode).replace('ftp://', "")
    hostname, directory = url.split('/', 1)
    ftp = ftplib.FTP(hostname, 'anonymous', "")

    filenames = [
        filename for filename, metadata in ftp.mlsd(directory) if metadata['type'] == 'file'
    ]
    if file_deck is not None:
        filenames = [filename for filename in filenames if filename[0] == file_deck.value]

    return sorted((filename.split('.')[0] for filename in filenames), reverse=True)
