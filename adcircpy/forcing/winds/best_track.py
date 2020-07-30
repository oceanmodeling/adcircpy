from datetime import datetime, timedelta
import gzip
from io import BytesIO
import pathlib
import urllib.error
import urllib.request

# import utm
from matplotlib import pyplot
from matplotlib.transforms import Bbox
import numpy
from pandas import DataFrame
from pyproj import CRS, Proj
from shapely.geometry import Point, Polygon

from adcircpy.forcing.winds import atcf_id
from adcircpy.forcing.winds.base import WindForcing


# import os
# from pathlib import Path
# import zipfile
# from adcircpy.lib._get_cache_directory import _get_cache_directory


class BestTrackForcing(WindForcing):
    def __init__(self, storm_id: str, start_date: datetime = None, end_date: datetime = None,
                 crs: CRS = None):
        self._storm_id = storm_id
        super().__init__(start_date, end_date, crs)

    def clip_to_bbox(self, bbox: Bbox):
        """
        Important: bbox must be expressed in Mercator projection (EPSG:3395)
        """
        assert isinstance(bbox, Bbox), f'bbox must be a {Bbox} instance.'
        bbox_pol = Polygon([
            [bbox.xmin, bbox.ymin],
            [bbox.xmax, bbox.ymin],
            [bbox.xmax, bbox.ymax],
            [bbox.xmin, bbox.ymax],
            [bbox.xmin, bbox.ymin]
        ])
        _switch = True
        unique_dates = numpy.unique(self._df['datetime'])
        for _datetime in unique_dates:
            records = self._df[self._df['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852. * radii  # convert to meters
            merc = Proj('EPSG:3395')
            x, y = merc(
                records['longitude'].iloc[0],
                records['latitude'].iloc[0])
            p = Point(x, y)
            pol = p.buffer(radii)
            if _switch:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _switch = False
                    continue
                    # self.start_date =
            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

    def plot_trajectory(self, ax: pyplot.Axes = None, show: bool = False, color='k', **kwargs):
        kwargs.update({'color': color})
        if ax is None:
            fig = pyplot.figure()
            ax = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * numpy.sin(numpy.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * numpy.cos(numpy.deg2rad(self.direction.iloc[i]))
            ax.quiver(
                self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            ax.annotate(
                self.df['datetime'].iloc[i],
                (self.longitude.iloc[i], self.latitude.iloc[i])
            )
        if show:
            ax.axis('scaled')
            pyplot.show()

    def write(self, path: str, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception('Files exist, set overwrite=True to allow overwrite.')
        with open(path, 'w') as f:
            f.write(self.fort22)

    @property
    def storm_id(self):
        return self._storm_id

    @property
    def _storm_id(self):
        return f"{self.basin}{self.storm_number}{self.year}"

    @_storm_id.setter
    def _storm_id(self, storm_id):
        chars = 0
        for char in storm_id:
            if char.isdigit():
                chars += 1

        if chars == 4:
            _atcf_id = atcf_id(storm_id)
            if _atcf_id is None:
                raise Exception(f'No storm with id: {storm_id}')
            storm_id = _atcf_id

        url = f'ftp://ftp.nhc.noaa.gov/atcf/archive/{storm_id[4:]}/b{storm_id[0:2].lower()}{storm_id[2:]}.dat.gz'
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.URLError:
            raise NameError(f'Did not find storm with id {storm_id} at url "{url}"')
        self.__atcf = BytesIO(response.read())

    @property
    def _start_date(self):
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is not None:
            assert isinstance(start_date, datetime)
        else:
            start_date = self._df['datetime'].iloc[0]
        assert self._df['datetime'].iloc[0] <= start_date < self._df['datetime'].iloc[-1], \
            f"start_date must be {self._df['datetime'].iloc[0]} <= start_date ({start_date}) < " \
            f"{self._df['datetime'].iloc[-1]}"
        self.__start_date = start_date

    @property
    def _end_date(self):
        return self.__end_date

    @_end_date.setter
    def _end_date(self, end_date):
        if end_date is not None:
            assert isinstance(end_date, datetime)
        else:
            end_date = self._df['datetime'].iloc[-1]
        assert self._df['datetime'].iloc[0] < end_date <= self._df['datetime'].iloc[-1], \
            f"end_date must be {self._df['datetime'].iloc[0]} <= end_date ({end_date}) <= " \
            f"{self._df['datetime'].iloc[-1]}"
        assert end_date > self.start_date, \
            f"end_date ({end_date}) must be after start_date ({self.start_date})"
        self.__end_date = end_date

    @property
    def _df(self):
        # https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt
        try:
            return self.__df
        except AttributeError:
            data = {
                'basin'                       : [],
                'storm_number'                : [],
                'datetime'                    : [],
                'record_type'                 : [],
                'latitude'                    : [],
                'longitude'                   : [],
                'max_sustained_wind_speed'    : [],
                'central_pressure'            : [],
                'development_level'           : [],
                'isotach'                     : [],
                'quadrant'                    : [],
                'radius_for_NEQ'              : [],
                'radius_for_SEQ'              : [],
                'radius_for_SWQ'              : [],
                'radius_for_NWQ'              : [],
                'background_pressure'         : [],
                'radius_of_last_closed_isobar': [],
                'radius_of_maximum_winds'     : [],
                'name'                        : [],
                'direction'                   : [],
                'speed'                       : []
            }
            for i, line in enumerate(gzip.GzipFile(fileobj=self.__atcf)):
                line = line.decode('UTF-8').split(',')
                data['basin'].append(line[0])
                data['storm_number'].append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == '':
                    _minutes = '00'
                _datetime = _datetime + _minutes
                data['datetime'].append(datetime.strptime(_datetime, '%Y%m%d%H%M'))
                data['record_type'].append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N ')) * 0.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S ')) * -0.1
                data['latitude'].append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E ')) * 0.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W ')) * -0.1
                data['longitude'].append(_lon)
                data['max_sustained_wind_speed'].append(float(line[8].strip(' ')))
                data['central_pressure'].append(float(line[9].strip(' ')))
                data['development_level'].append(line[10].strip(' '))
                data['isotach'].append(int(line[11].strip(' ')))
                data['quadrant'].append(line[12].strip(' '))
                data['radius_for_NEQ'].append(int(line[13].strip(' ')))
                data['radius_for_SEQ'].append(int(line[14].strip(' ')))
                data['radius_for_SWQ'].append(int(line[15].strip(' ')))
                data['radius_for_NWQ'].append(int(line[16].strip(' ')))
                if len(line) > 18:
                    data['background_pressure'].append(int(line[17].strip(' ')))
                    data['radius_of_last_closed_isobar'].append(int(line[18].strip(' ')))
                    data['radius_of_maximum_winds'].append(int(line[19].strip(' ')))
                    if len(line) > 23:
                        data['name'].append(line[27].strip(' '))
                    else:
                        data['name'].append('')
                else:
                    data['background_pressure'].append(data['background_pressure'][-1])
                    data['radius_of_last_closed_isobar'].append(
                        data['radius_of_last_closed_isobar'][-1])
                    data['radius_of_maximum_winds'].append(data['radius_of_maximum_winds'][-1])
                    data['name'].append('')
            data = self._compute_velocity(data)
            # data = self._transform_coordinates(data)
            self.__df = DataFrame(data=data)
            return self.__df

    @property
    def fort22(self):
        record_number = self._generate_record_numbers()
        fort22 = ''
        for i, (_, row) in enumerate(self.df.iterrows()):
            longitude = row['longitude']
            latitude = row['latitude']
            if longitude >= 0:
                longitude = f'{int(longitude / 0.1):>5}E'
            else:
                longitude = f'{int(longitude / -0.1):>5}W'
            if latitude >= 0:
                latitude = f'{int(latitude / 0.1):>4}N'
            else:
                latitude = f'{int(latitude / -0.1):>4}S'
            row['longitude'] = longitude
            row['latitude'] = latitude

            background_pressure = row['background_pressure']
            if background_pressure is None:
                background_pressure = self.df['background_pressure'].iloc[i - 1]
            if background_pressure is not None:
                if background_pressure <= row['central_pressure'] < 1013:
                    background_pressure = 1013
                elif background_pressure <= row['central_pressure'] >= 1013:
                    background_pressure = int(row['central_pressure'] + 1)
                else:
                    background_pressure = int(row['background_pressure'])
            else:
                background_pressure = ''
            row['background_pressure'] = background_pressure

            output_row = [
                # BASIN - basin, e.g. WP, IO, SH, CP, EP, AL, SL
                f'{row["basin"]:<2}',
                # CY - annual cyclone number: 1 through 99
                f'{row["storm_number"]:>3}',
                # YYYYMMDDHH - Warning Date-Time-Group: 0000010100 through 9999123123. (note, 4 digit year)
                f'{format(row["datetime"], "%Y%m%d%H"):>11}',
                # TECHNUM/MIN - objective technique sorting number, minutes for best track: 00 - 99
                f'{"":3}',
                # TECH - acronym for each objective technique or CARQ or WRNG, BEST for best track.
                f'{row["record_type"]:>5}',
                # TAU - forecast period: -24 through 240 hours, 0 for best-track, negative taus used for CARQ and WRNG records.
                f'{int((row["datetime"] - self.start_date) / timedelta(hours=1)):>4}',
                # LatN/S - Latitude (tenths of degrees) for the DTG: 0 through 900, N/S is the hemispheric index.
                f'{latitude:>5}',
                # LonE/W - Longitude (tenths of degrees) for the DTG: 0 through 1800, E/W is the hemispheric index.
                f'{longitude:>5}',
                # VMAX - Maximum sustained wind speed in knots: 0 through 300.
                f'{int(row["max_sustained_wind_speed"]):>4}',
                # MSLP - Minimum sea level pressure, 1 through 1100 MB.
                f'{int(row["central_pressure"]):>5}',
                # TY - Level of tc development:
                f'{row["development_level"] if row["development_level"] is not None else "":>3}',
                # RAD - Wind intensity (kts) for the radii defined in this record: 34, 50, 64.
                f'{int(row["isotach"]) if row["isotach"] is not None else "":>4}',
                # WINDCODE - Radius code: AAA - full circle, QQQ - quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
                f'{row["quadrant"] if row["quadrant"] is not None else "":>4}',
                # RAD1 - If full circle, radius of specified wind intensity, If semicircle or quadrant, radius of specified wind intensity of circle portion specified in radius code. 0 - 1200 nm.
                f'{int(row["radius_for_NEQ"]) if row["radius_for_NEQ"] is not None else "":>5}',
                # RAD2 - If full circle this field not used, If semicicle, radius (nm) of specified wind intensity for semicircle not specified in radius code, If quadrant, radius (nm) of specified wind intensity for 2nd quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.
                f'{int(row["radius_for_SEQ"]) if row["radius_for_SEQ"] is not None else "":>5}',
                # RAD3 - If full circle or semicircle this field not used, If quadrant, radius (nm) of specified wind intensity for 3rd quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.
                f'{int(row["radius_for_SWQ"]) if row["radius_for_SWQ"] is not None else "":>5}',
                # RAD4 - If full circle or semicircle this field not used, If quadrant, radius (nm) of specified wind intensity for 4th quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.
                f'{int(row["radius_for_NWQ"]) if row["radius_for_NWQ"] is not None else "":>5}',
                # RADP - pressure in millibars of the last closed isobar, 900 - 1050 mb.
                f'{row["background_pressure"]:>5}',
                # RRP - radius of the last closed isobar in nm, 0 - 9999 nm.
                f'{int(row["radius_of_last_closed_isobar"]) if row["radius_of_last_closed_isobar"] is not None else "":>5}',
                # MRD - radius of max winds, 0 - 999 nm.
                f'{int(row["radius_of_maximum_winds"]) if row["radius_of_maximum_winds"] is not None else "":>4}'
                # GUSTS - gusts, 0 through 995 kts.
                f'{"":>5}',
                # EYE - eye diameter, 0 through 999 nm.
                f'{"":>4}',
                # SUBREGION - subregion code: W, A, B, S, P, C, E, L, Q.
                f'{"":>4}',
                # MAXSEAS - max seas: 0 through 999 ft.
                f'{"":>4}',
                # INITIALS - Forecaster's initials, used for tau 0 WRNG, up to 3 chars.
                f'{"":>4}',
                # DIR - storm direction in compass coordinates, 0 - 359 degrees.
                f'{row["direction"] if row["direction"] is not None else "":>3}',
                # SPEED - storm speed, 0 - 999 kts.
                f'{row["speed"]:>4}',
                # STORMNAME - literal storm name, NONAME or INVEST. TCcyx used pre-1999, where:
                f'{row["name"]:^12}',
                # from this point forwards it's all aswip
                f'{record_number[i]:>4}'
            ]
            fort22 += f'{",".join(output_row)}\n'
        return fort22

    @property
    def WTIMINC(self):
        return f'{self.start_date:%Y %m %d %H} {self.df["storm_number"].iloc[0]} ' \
               f'{self.BLADj} {self.geofactor}'

    @property
    def BLADj(self):
        try:
            return self.__BLADj
        except AttributeError:
            return 0.9

    @BLADj.setter
    def BLADj(self, BLADj: float):
        BLADj = float(BLADj)
        assert 0 <= BLADj <= 1
        self.__BLADj = BLADj

    @property
    def geofactor(self):
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @geofactor.setter
    def geofactor(self, geofactor: float):
        geofactor = float(geofactor)
        assert 0 <= geofactor <= 1
        self.__geofactor = geofactor

    def _generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i - 1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number
