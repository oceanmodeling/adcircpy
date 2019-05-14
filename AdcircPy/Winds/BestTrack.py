import os
import io
import gzip
import urllib
import numpy as np
from datetime import datetime
from AdcircPy.Winds._WindVortex import _WindVortex


class BestTrack(_WindVortex):
    def __init__(self, storm_id, start_date=None, end_date=None):
        self.__init_url(storm_id)
        super(BestTrack, self).__init__(storm_id, start_date, end_date)

    @property
    def url(self):
        return self._url

    def __init_url(self, storm_id):
        self._url = 'http://ftp.nhc.noaa.gov/atcf/archive/'
        self._url = self._url+storm_id[4:]+'/b'+storm_id[0:2].lower()+storm_id[2:]+'.dat.gz'
