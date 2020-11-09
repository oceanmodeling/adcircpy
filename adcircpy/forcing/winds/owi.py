from datetime import datetime

import numpy as np

from adcircpy.forcing.winds.base import WindForcing


class OwiForcing(WindForcing):
    def __init__(self):
        self.__basin_scale_pressure = None
        self.__basin_scale_winds = None
        self.__regional_scale_pressure = None
        self.__regional_scale_winds = None

    def make_plot(self):
        pass

    def __set_datetime(self, datetime):
        msg = 'Dates of input files provided do not match.'
        if self.datetime is not None:
            assert np.array_equal(self.datetime, datetime), msg
        else:
            self.__datetime = np.asarray(datetime)

    @staticmethod
    def __parse_fort22_p(file):
        with open(file, 'r') as f:
            OWI = dict()
            OWI['datetime'] = list()
            OWI['values'] = list()
            OWI['header'] = f.readline().strip('\n')
            line = f.readline()
            OWI['iLat'] = int(line[6:9])
            OWI['iLon'] = int(line[15:19])
            OWI['DX'] = float(line[22:28])
            OWI['DY'] = float(line[31:37])
            OWI['SWLat'] = float(line[45:51])
            OWI['SWLon'] = float(line[57:65])
            OWI['datetime'].append(
                datetime.strptime(line[68:80], '%Y%m%d%H%M'))
            values = list()
            for line in f:
                if 'iLat' in line:
                    OWI['datetime'].append(
                        datetime.strptime(line[68:80], '%Y%m%d%H%M'))
                    OWI['values'].append(
                        np.asarray(values).reshape(
                            (OWI['iLon'], OWI['iLat'])))
                    values = list()
                else:
                    for n in [line[i:i + 10] for i in range(0, 80, 10)]:
                        values.append(float(n))
        return OWI

    @staticmethod
    def __parse_fort22_w(file):
        with open(file, 'r') as f:
            OWI = dict()
            OWI['datetime'] = list()
            OWI['values'] = dict()
            OWI['values']['u'] = list()
            OWI['values']['v'] = list()
            OWI['header'] = f.readline().strip('\n')
            line = f.readline()
            OWI['iLat'] = int(line[6:9])
            OWI['iLon'] = int(line[15:19])
            OWI['DX'] = float(line[22:28])
            OWI['DY'] = float(line[31:37])
            OWI['SWLat'] = float(line[45:51])
            OWI['SWLon'] = float(line[57:65])
            OWI['datetime'].append(
                datetime.strptime(line[68:80], '%Y%m%d%H%M'))
            values_u = list()
            values_v = list()
            shape = (OWI['iLon'], OWI['iLat'])
            size = OWI['iLon'] * OWI['iLat']
            for line in f:
                if 'iLat' in line:
                    OWI['datetime'].append(
                        datetime.strptime(line[68:80], '%Y%m%d%H%M'))
                    OWI['values']['u'].append(
                        np.asarray(values_u).reshape(shape))
                    OWI['values']['v'].append(
                        np.asarray(values_v).reshape(shape))
                    values_u = list()
                    values_v = list()
                else:
                    for n in [line[i:i + 10] for i in range(0, 80, 10)]:
                        if len(values_u) != size:
                            values_u.append(float(n))
                        else:
                            values_v.append(float(n))
        return OWI

    @property
    def start_date(self):
        try:
            return self.__start_date
        except AttributeError:
            return

    @property
    def end_date(self):
        try:
            return self.__end_date
        except AttributeError:
            return

    @property
    def datetime(self):
        try:
            return self.__datetime
        except AttributeError:
            return

    @property
    def basin_scale_pressure(self):
        if self.__basin_scale_pressure is None:
            raise AttributeError('Must set basin_scale_pressure attribute.')
        return self.__basin_scale_pressure

    @property
    def basin_scale_winds(self):
        if self.__basin_scale_winds is None:
            raise AttributeError('Must set basin_scale_winds attribute.')
        return self.__basin_scale_winds

    @property
    def regional_scale_pressure(self):
        if self.__regional_scale_pressure is None:
            raise AttributeError('Must set regional_scale_pressure attribute.')
        return self.__regional_scale_pressure

    @property
    def regional_scale_winds(self):
        if self.__regional_scale_winds is None:
            raise AttributeError('Must set regional_scale_winds attribute.')
        return self.__regional_scale_winds

    @property
    def NWS(self):
        return 12

    @property
    def fort22(self):
        raise NotImplementedError

    @property
    def fort221(self):
        raise NotImplementedError

    @property
    def fort222(self):
        raise NotImplementedError

    @property
    def fort223(self):
        raise NotImplementedError

    @property
    def fort224(self):
        raise NotImplementedError

    @basin_scale_pressure.setter
    def basin_scale_pressure(self, basin_scale_pressure):
        _ = self.__parse_fort22_p(basin_scale_pressure)
        self.__set_datetime(_['datetime'])

    @basin_scale_winds.setter
    def basin_scale_winds(self, basin_scale_winds):
        _ = self.__parse_fort22_w(basin_scale_winds)
        self.__set_datetime(_['datetime'])

    @regional_scale_pressure.setter
    def regional_scale_pressure(self, regional_scale_pressure):
        _ = self.__parse_fort22_p(regional_scale_pressure)
        self.__set_datetime(_['datetime'])

    @regional_scale_winds.setter
    def regional_scale_winds(self, regional_scale_winds):
        _ = self.__parse_fort22_w(regional_scale_winds)
        self.__set_datetime(_['datetime'])

    @fort221.setter
    def fort221(self, fort221):
        self.basin_scale_pressure = fort221

    @fort222.setter
    def fort222(self, fort222):
        self.basin_scale_winds = fort222

    @fort223.setter
    def fort223(self, fort223):
        self.regional_scale_pressure = fort223

    @fort224.setter
    def fort224(self, fort224):
        self.regional_scale_winds = fort224
