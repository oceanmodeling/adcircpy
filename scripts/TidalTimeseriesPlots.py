#! /usr/bin/env python
"""
Standalone program for quickly looking at fort.61.nc files against COOPS data.
For a complete list of options use ./TidalTimeseriesPlots.py -h

Example usage:
    ./TidalTimeseriesPlots.py /path/to/fort.61.nc --save-path /path/to/directory/for/saving/plots
"""

import os
import argparse
import matplotlib.pyplot as plt
from AdcircPy import AdcircPy
from AdcircPy.Validation import COOPS

class TidalTimeseriesPlots(object):
  def __init__(self):
    self._parse_args()
    self._read_fort61()
    self._get_coops_timeseries()
    self._generate_plots()

  def _parse_args(self):
    parser = argparse.ArgumentParser(description="Program to see a quick plot of an ADCIRC mesh.")
    parser.add_argument("fort61", help="Path to ADCIRC fort.61 or fort.61.nc file.")
    show = parser.add_mutually_exclusive_group(required=False)
    show.add_argument('--show', dest='show', action='store_true', help='Shows plots to screen as they are generated (default).')
    show.add_argument('--no-show', dest='show', action='store_false', help='Prevents the plots from showing to screen. Useful for only saving the plots without showing them.')
    parser.set_defaults(show=True)
    parser.add_argument('--save-path', help="Directory where to save plots. Will be created if it doesn't exist.")
    self.args =  parser.parse_args()

  def _read_fort61(self):
    print('\n-> Reading ADCIRC tidal timeseries output file...')
    self.fort61 = AdcircPy.read_output(self.args.fort61)
    
  def _get_coops_timeseries(self):
    print('-> Fetching station data from COOPS...')
    self.coops = COOPS.TidalStations(self.fort61.keys(), self.fort61.time[0], self.fort61.time[-1])
    
  def _generate_plots(self):
    self.__init_save_dir()
    print('-> Generating plots...')
    for station in self.coops.keys():
      plt.figure(figsize=(15,11))
      plt.subplot(111)
      plt.plot(self.coops[station]['time'], self.coops[station]['zeta'], label='Observation')
      plt.plot(self.fort61.time, self.fort61[station]['zeta'], label='ADCIRC')
      plt.gca().set_title(self.coops[station]['metadata']['name'])
      plt.legend(loc='best')
      if self.args.show==True:
        plt.show()
      if self.args.save_path is not None:
        plt.savefig(self.args.save_path+'/{}_{}.png'.format(station, self.coops[station]['metadata']['name']))
      plt.close(plt.gcf())
  
  def __init_save_dir(self):
    if self.args.save_path is not None:
      os.makedirs(self.args.save_path, exist_ok=True)

if __name__ == "__main__":
  TidalTimeseriesPlots()