#! /usr/bin/env python
"""
Entry point for quickly looking at fort.61.nc files against COOPS data.
For a complete list of options use ./TidalTimeseriesStationsValidation.py -h

Example usage:
    TidalTimeseriesStationsValidation /path/to/fort.61.nc --save-path /path/to/directory/for/saving/plots
"""

import os
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from AdcircPy import AdcircPy
from AdcircPy.Validation import COOPS


class PlotTidalStationsOutput:
    def __init__(self):
        self.__set_args()
        self.__set_fort61()
        self.__set_COOPS()
        # self.__set_savedir()
        self.__generate_plots()

    def __set_args(self):
        parser = argparse.ArgumentParser(description="Program to see a quick plot of an ADCIRC mesh.")
        parser.add_argument("fort61", help="Path to ADCIRC fort.61 or fort.61.nc file.")
        show = parser.add_mutually_exclusive_group(required=False)
        show.add_argument('--show', dest='show', action='store_true', help='Shows plots to screen as they are generated (default).')
        show.add_argument('--no-show', dest='show', action='store_false', help='Prevents the plots from showing to screen. Useful for only saving the plots without showing them.')
        parser.set_defaults(show=True)
        parser.add_argument('--save-path', help="Directory where to save plots. Will be created if it doesn't exist.")
        self.args = parser.parse_args()

    def __set_fort61(self):
        self.fort61 = AdcircPy.read_output(self.args.fort61)

    def __set_COOPS(self):
        self.COOPS = COOPS.TidalStations()
        for station_id, station_data in self.fort61.items():
            start_date = station_data.time[0]
            end_date = station_data.time[-1]
            self.COOPS.fetch(station_id, start_date, end_date)

    def __set_save_dir(self):
        if self.args.save_path is not None:
            os.makedirs(self.args.save_path, exist_ok=True)

    def __generate_plots(self):

        for station, data in self.fort61.items():

            if station in self.COOPS.keys():
                plt.figure(figsize=(15, 11))
                plt.subplot(111)
                station_name = station
                coops_dates = [time.timestamp()
                               for time in self.COOPS[station].time]
                f = interp1d(coops_dates, self.COOPS[station].values)
                data_dates = [time.timestamp() for time in data.time]
                coops_values = f(data_dates)
                plt.plot(data.time, data.values, label='ADCIRC', color='r')
                plt.plot(data.time, coops_values, label='Observation',
                         color='b')
                plt.plot(data.time, data.values-coops_values, color='k',
                         label='difference')
                station_name = '{} '.format(self.COOPS[station].name)
                station_name += '{}'.format(station)
                plt.title('{}'.format(station_name))
                plt.legend()
                plt.show()
                # plt.gca().set_title(self.coops[station]['metadata']['name'])
                # plt.legend(loc='best')
                # if self.args.show is True:
                #     plt.show()
                # if self.args.save_path is not None:
                #     plt.savefig(self.args.save_path+'/{}_{}.png'.format(station, self.coops[station]['metadata']['name']))
                plt.close(plt.gcf())


def main():
    PlotTidalStationsOutput()


if __name__ == "__main__":
    main()
