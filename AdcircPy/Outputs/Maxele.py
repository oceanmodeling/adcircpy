import argparse
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs.ScalarSurfaceExtrema import ScalarSurfaceExtrema


class Maxele(ScalarSurfaceExtrema):
    def __init__(self, x, y, elements, values, times, epsg=4326,
                 vertical_datum='LMSL', **kwargs):
        super(Maxele, self).__init__(x, y, elements, values, times,
                                     epsg, vertical_datum, **kwargs)

    @classmethod
    def from_netcdf(cls, path, fort14=None, vertical_datum='LMSL',
                    epsg=4326, datum_grid=None):
        nc = Dataset(path)
        if 'zeta_max' not in nc.variables.keys():
            raise Exception('Not a maxele file!')
        fort14 = cls.__init_fort14(fort14)
        return cls(nc['x'][:],
                   nc['y'][:],
                   nc['element'][:]-1,
                   nc['zeta_max'][:],
                   nc['time_of_zeta_max'],
                   epsg,
                   vertical_datum,
                   # **fort14.get_boundary_dictionary() # Need to add the additional infor provided by fort.14
                   datum_grid=datum_grid)

    @classmethod
    def from_ascii(cls, path, fort14, vertical_datum='LMSL',
                   epsg=4326, datum_grid=None):
        fort14 = cls.__init_fort14(fort14)
        with open(path, 'r') as f:
            line = f.readline()
            line = f.readline().split()
            NP = int(line[1])
            line = f.readline()
            values = list()
            for i in range(NP):
                values.append(float(f.readline().split()[1]))
            time_of_zeta_max = list()
            for i in range(NP):
                try:
                    time_of_zeta_max.append(float(f.readline().split()[1]))
                except IndexError:
                    time_of_zeta_max.append(-99999.)

        values = np.ma.masked_equal(values, -99999.)
        time_of_zeta_max = np.ma.masked_equal(time_of_zeta_max, -99999.)
        return cls(fort14.x,
                   fort14.y,
                   fort14.elements,
                   values,
                   time_of_zeta_max,
                   epsg=epsg,
                   vertical_datum=vertical_datum,
                   ocean_boundaries=fort14.ocean_boundaries,
                   land_boundaries=fort14.land_boundaries,
                   inner_boundaries=fort14.inner_boundaries,
                   weir_boundaries=fort14.weir_boundaries,
                   inflow_boundaries=fort14.inflow_boundaries,
                   outflow_boundaries=fort14.outflow_boundaries,
                   culvert_boundaries=fort14.culvert_boundaries,
                   datum_grid=datum_grid)

    @staticmethod
    def __init_fort14(fort14):
        if fort14 is not None:
            if isinstance(fort14, AdcircMesh) is False:
                return AdcircMesh.from_fort14(fort14)


def PlotMaxeleEntrypoint():
    class PlotMaxele(object):

        def __init__(self):
            self.parse_args()
            if self.check_is_netcdf():
                maxele = Maxele.from_netcdf(self.args.maxele)
            else:
                if self.args.fort14 is None:
                    raise Exception('A fort.14 file is required if maxele file'
                                    + ' is ascii.')
                maxele = Maxele.from_ascii(self.args.maxele, self.args.fort14)
            maxele.make_plot(title=self.args.title,
                             vmin=self.args.vmin,
                             vmax=self.args.vmax)
            plt.show()

        def parse_args(self):
            parser = argparse.ArgumentParser(
                        description="Program to see a quick plot of an ADCIRC "
                        + "maxele file.")
            parser.add_argument('maxele', help="Path to maxele file.")
            parser.add_argument('--fort14', help="Path to fort.14 file"
                                + "(required if maxele files is not netcdf).")
            parser.add_argument('--title', help="Plot title override.")
            parser.add_argument('--vmin', type=float)
            parser.add_argument('--vmax', type=float)
            self.args = parser.parse_args()

        def check_is_netcdf(self):
            try:
                Dataset(self.args.maxele)
                return True
            except FileNotFoundError:
                raise
            except OSError:
                return False
    PlotMaxele()
