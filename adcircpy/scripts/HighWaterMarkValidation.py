#! /usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import csv
import pandas as pd
from scipy import stats, spatial
import seaborn as sns
import numpy as np
from AdcircPy.Validation.USGS import HighWaterMarks
from AdcircPy.Outputs import Maxele
import warnings


class HighWaterMarkValidation(object):
    def __init__(self, args):
        self.args = args
        self._init_maxele()
        self._init_HighWaterMarks()
        self._init_dataframe()
        self._init_Validation()
        self._init_wl_threshold()
        self._init_CDH()
        self._init_POD_and_FAR()
        self._init_statistics()
        self._init_best_case_statistics()
        self._init_output_stats()
        self._print_output_stats()
        self._save_output_stats()
        self._dump_dataframe()

    def _init_maxele(self):
        if Maxele.is_netcdf(self.args.maxele):
            self.maxele = Maxele.from_netcdf(self.args.maxele)
        else:
            if self.args.fort14 is None:
                raise IOError("A mesh file is required if the maxele file is\
                              ASCII format.")
            self.maxele = Maxele.from_ascii(self.args.maxele, self.args.fort14,
                                            datum_grid=self.args.datum_grid)
        if self.args.units == 'feet':
            self.maxele._values *= 3.28084
        self.maxele.convert_datum(self.args.navd88_mesh, 'NAVD88')

    def _init_HighWaterMarks(self):
        """ """
        if self.args.csv is not None:
            self.HighWaterMarks = HighWaterMarks.from_csv(self.args.csv)
        else:
            self.HighWaterMarks = HighWaterMarks.from_event_name(
                                                  self.args.event_id)

    def _init_dataframe(self):
        if self.args.method == 'nearest':
            self.__init_nearest()
        elif self.args.method == 'nanmean':
            self.__init_nanmean()
        elif self.args.method == 'idw':
            self.__init_idw()

    def _init_Validation(self):
        self.__init_regression_plot()
        self.__init_regression_plot_best_case()

    def _init_wl_threshold(self):
        if self.args.units == 'feet':
            self.wl_thresholds = np.arange(3, 11)
        else:
            self.wl_thresholds = np.linspace(1, 4, num=6)

    def _init_CDH(self):
        sns.set_context("talk")
        CDH = list()
        for wl_threshold in self.wl_thresholds:
            cnt = 0
            for station in self.missed:
                if self.args.units == 'feet':
                    if self.HighWaterMarks[station]['elev_ft'] >= wl_threshold:
                        cnt += 1
                else:
                    if self.HighWaterMarks[station]['elev_ft']/3.28084 >= \
                                                                 wl_threshold:
                        cnt += 1
            CDH.append(cnt)
        plt.plot(self.wl_thresholds, CDH, color='gold', marker='.')
        title = 'Count of dry forecast points for {}, {}'.format(
                                                self.HighWaterMarks.event_name,
                                                self.HighWaterMarks.event_year)
        if self.args.subtitle is not None:
            title += '\n {}'.format(self.args.subtitle)
        plt.title(title)
        plt.ylim([0, np.max(CDH)+2])
        plt.xlabel('Observed Water Level Threshold [{} NAVD88]'.format(
                                                            self.args.units))
        plt.ylabel('Count Dry Forecasts')
        self.__save_plot('CDH.png')
        self.__show_plot()
        plt.close(plt.gcf())

    def _init_POD_and_FAR(self):
        sns.set_context("talk")
        self._POD = list()
        self._FAR = list()
        self._count = list()
        Obs = np.asarray(self.dataframe['Observation'])
        Adc = np.asarray(self.dataframe['ADCIRC'])
        # the total number of events
        for Thres in self.wl_thresholds:
            hit = np.where(np.logical_and(Adc >= Thres, Obs >= Thres))
            miss = np.where(np.logical_and(Adc < Thres, Obs >= Thres))
            fa = np.where(np.logical_and(Adc >= Thres, Obs < Thres))
            self._count.append(str(Obs[hit].size+Obs[miss].size))
            if (Obs[hit].size+np.float(Obs[miss].size)) > 0:
                self._POD.append(np.float(Obs[hit].size) /
                                 (np.float(Obs[hit].size) +
                                  np.float(Obs[miss].size)))
            else:
                self._POD.append(np.nan)
            if (np.float(Obs[fa].size)+np.float(Obs[hit].size)) > 0:
                self._FAR.append(np.float(Obs[fa].size) /
                                 (np.float(Obs[fa].size) +
                                  np.float(Obs[hit].size)))
            else:
                self._FAR.append(np.nan)
        plt.plot(self.wl_thresholds, self._POD, color='darkgreen', label='POD')
        plt.plot(self.wl_thresholds, self._FAR, color='darkblue', label='FAR')
        plt.plot(self.wl_thresholds, self.wl_thresholds*[0],
                 color='black', linestyle='--')
        for i, count in enumerate(self._count):
            plt.annotate(count, (self.wl_thresholds[i], 0),
                         verticalalignment='top', horizontalalignment='center')
        plt.ylim([-0.2, 1.2])
        plt.yticks(np.arange(0, 1.2, step=0.2))
        plt.legend(loc='best')
        title = 'Probability of detection for {}, {}'.format(
                                                self.HighWaterMarks.event_name,
                                                self.HighWaterMarks.event_year)
        if self.args.subtitle is not None:
            title += '\n {}'.format(self.args.subtitle)
        plt.title(title)
        plt.xlabel('Water Level Threshold [{} NAVD88]'.format(self.args.units))
        plt.ylabel('Skill')
        self.__save_plot('POD_FAR.png')
        self.__show_plot()
        plt.close(plt.gcf())

    def _init_statistics(self):
        self._statistics = dict()
        self._statistics['pearson_r'] = self.__get_pearson_r(
                                            self.dataframe["Observation"])
        self._statistics['bias'] = self.__get_bias(
                                            self.dataframe["Observation"])
        self._statistics['RMS_demeaned'] = self.__get_RMS_demeaned(
                                            self.dataframe["Observation"],
                                            self._statistics['bias'])
        self._statistics['scatter_index'] = self.__get_scatter_index(
                                            self.dataframe["Observation"],
                                            self._statistics['RMS_demeaned'])
        self._statistics['90percent_index'] = self.__get_90percent_index(
                                                self.dataframe["Observation"])

    def _init_best_case_statistics(self):
        self._best_case_statistics = dict()
        self._best_case_statistics['pearson_r'] = self.__get_pearson_r(
                                                self.dataframe["best_case"])
        self._best_case_statistics['bias'] = self.__get_bias(
                                                self.dataframe["best_case"])
        self._best_case_statistics['RMS_demeaned'] = self.__get_RMS_demeaned(
                                            self.dataframe["best_case"],
                                            self._best_case_statistics['bias'])
        self._best_case_statistics['scatter_index'] \
            = self.__get_scatter_index(
                            self.dataframe["best_case"],
                            self._best_case_statistics['RMS_demeaned'])
        self._best_case_statistics['90percent_index'] \
            = self.__get_90percent_index(self.dataframe["best_case"])

    def _init_output_stats(self):
        self.output_stats = list()
        self.output_stats.append("Pairing method: {}\n\n"
                                 .format(self.args.method))
        self.output_stats.append("Filtered points = {}\n".format(self.HighWaterMarks.filtered_count))
        self.output_stats.append("Missed points = {}\n".format(len(self.missed)))
        self.output_stats.append("Valid points = {}\n\n".format(len(self.dataframe)))
        self.output_stats.append("Pearson $r^2$ = {:.3f}\n".format(self._statistics['pearson_r']**2))
        self.output_stats.append("Bias = {:.3f}\n".format(self._statistics['bias']))
        self.output_stats.append("RMS demeaned = {:.3f}\n".format(self._statistics['RMS_demeaned']))
        self.output_stats.append("Scatter index = {:.3f}\n".format(self._statistics['scatter_index']))
        self.output_stats.append("90% index = {:.3f}\n\n".format(self._statistics['90percent_index']))
        self.output_stats.append("Pearson $r^2$ (best case) = {:.3f}\n".format(self._best_case_statistics['pearson_r']**2))
        self.output_stats.append("Bias (best case) = {:.3f}\n".format(self._best_case_statistics['bias']))
        self.output_stats.append("RMS demeaned (best case) = {:.3f}\n".format(self._best_case_statistics['RMS_demeaned']))
        self.output_stats.append("Scatter index (best case) = {:.3f}\n".format(self._best_case_statistics['scatter_index']))
        self.output_stats.append("90% index (best case) = {:.3f}".format(self._best_case_statistics['90percent_index']))

    def _save_output_stats(self):
        if self.args.save_path is not None:
            with open(self.args.save_path+'/statistics.txt', 'w') as f:
                for line in self.output_stats:
                    f.write(line)

    def _print_output_stats(self):
        if self.args.show is True:
            for line in self.output_stats:
                print(line, end='')

    def _dump_dataframe(self, filename='data.csv'):
        if self.args.save_path is not None:
            fname = self.args.save_path + '/' + filename
            with open(fname, 'w') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',')
                csvwriter.writerow(['longitude', 'latitude', 'adcirc', 'hwm',
                                    'uncertainty', 'units'])
                for i, station in enumerate(self.dataframe['station_id']):
                    row = list()
                    row.append(self.HighWaterMarks[station]['longitude'])
                    row.append(self.HighWaterMarks[station]['latitude'])
                    row.append(self.dataframe['ADCIRC'][i])
                    row.append(self.dataframe['Observation'][i])
                    row.append(self.dataframe['Uncertainty'][i])
                    row.append(self.args.units)
                    csvwriter.writerow(row)

    def __init_values_dict(self):
        self._values = {"Observation": list(),
                        "ADCIRC": list(),
                        "Quality_ID": list(),
                        "Uncertainty": list(),
                        "station_id": list(),
                        "best_case": list()}

    def __init_nearest(self):
        self.missed = list()
        self.__init_values_dict()
        for station in self.HighWaterMarks.keys():
            x = self.HighWaterMarks[station]['longitude']
            y = self.HighWaterMarks[station]['latitude']
            element = self.maxele.get_element_containing_coord((x,y))
            if element is None:
                self.missed.append(station)
            else:
                tree = spatial.KDTree(np.vstack([self.maxele.x[element], self.maxele.y[element]]).T)
                distances, indexes = tree.query([x,y], k=3)
                if self.maxele.values[element][indexes].mask.all():
                    self.missed.append(station)
                else:
                    for value in self.maxele.values[element][indexes]:
                        if np.ma.is_masked(value)==False:
                            break
                        else:
                            continue
                    self._values["ADCIRC"].append(value)
                    if self.args.units == 'meters':
                        self._values["Observation"].append(self.HighWaterMarks[station]["elev_ft"]/3.28084)
                    else:
                        self._values["Observation"].append(self.HighWaterMarks[station]["elev_ft"])
                    self._values["Quality_ID"].append(self.HighWaterMarks[station]["hwmQualityName"])
                    self._values["station_id"].append(station)
                    self.__init_uncertainties(station, value)
        self.dataframe = pd.DataFrame(self._values)

    def __init_nanmean(self):
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        self.missed = list()
        self.__init_values_dict()
        for station in self.HighWaterMarks.keys():
            x = self.HighWaterMarks[station]['longitude']
            y = self.HighWaterMarks[station]['latitude']
            element = self.maxele.get_element_containing_coord((x, y))
            if element is None:
                self.missed.append(station)
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    value = np.nanmean(np.ma.filled(
                                        self.maxele.values[element],
                                        fill_value=np.nan))
                if np.isnan(value) is True:
                    self.missed.append(station)
                else:
                    if self.args.units == 'meters':
                        self._values["Observation"].append(
                            self.HighWaterMarks[station]["elev_ft"]/3.28084)
                    else:
                        self._values["Observation"].append(
                            self.HighWaterMarks[station]["elev_ft"])
                    self._values["ADCIRC"].append(value)
                    self._values["Quality_ID"].append(
                                self.HighWaterMarks[station]["hwmQualityName"])
                    self._values["station_id"].append(station)
                    self.__init_uncertainties(station, value)
        self.dataframe = pd.DataFrame(self._values)

    def __init_idw(self):
        self.missed = list()
        self.__init_values_dict()
        for station in self.HighWaterMarks.keys():
            x = self.HighWaterMarks[station]['longitude']
            y = self.HighWaterMarks[station]['latitude']
            element = self.maxele.get_element_containing_coord((x, y))
            if element is None:
                self.missed.append(station)
            else:
                tree = spatial.KDTree(np.vstack([self.maxele.x[element],
                                                self.maxele.y[element]]).T)
                distances, indexes = tree.query([x, y], k=3)
                if self.maxele.values[element][indexes].mask.all():
                    self.missed.append(station)
                else:
                    _values = list()
                    _distances = list()
                    for i in range(3):
                        value = self.maxele.values[element][indexes[i]]
                        if np.ma.is_masked(value) is False:
                            _values.append(value)
                            _distances.append(distances[i])
                    value = np.divide(_values, _distances).sum() / \
                                     (1./np.asarray(_distances)).sum()
                    if self.args.units == 'meters':
                        self._values["Observation"].append(
                            self.HighWaterMarks[station]["elev_ft"]/3.28084)
                    else:
                        self._values["Observation"].append(
                            self.HighWaterMarks[station]["elev_ft"])
                    self._values["ADCIRC"].append(value)
                    self._values["Quality_ID"].append(
                            self.HighWaterMarks[station]["hwmQualityName"])
                    self._values["station_id"].append(station)
                    self.__init_uncertainties(station, value)
        self.dataframe = pd.DataFrame(self._values)

    def __init_uncertainties(self, station, value):
        if "Excellent" in self.HighWaterMarks[station]["hwmQualityName"]:
            uncertainty = 0.05
        elif "Good" in self.HighWaterMarks[station]["hwmQualityName"]:
            uncertainty = 0.1
        elif "Fair" in self.HighWaterMarks[station]["hwmQualityName"]:
            uncertainty = 0.2
        elif "Poor" in self.HighWaterMarks[station]["hwmQualityName"]:
            uncertainty = 0.4
        else:
            uncertainty = np.nan
        if self.args.units == 'meters':
            uncertainty /= 3.28084
        self._values["Uncertainty"].append(uncertainty)
        self.__init_best_case(station, uncertainty, value)

    def __init_best_case(self, station, uncertainty, value):
        if self.args.units == 'feet':
            cases = [(self.HighWaterMarks[station]["elev_ft"] + uncertainty),
                     (self.HighWaterMarks[station]["elev_ft"] - uncertainty)]
        else:
            cases = [((self.HighWaterMarks[station]["elev_ft"]/3.28084) +
                      uncertainty),
                     ((self.HighWaterMarks[station]["elev_ft"]/3.28084) -
                      uncertainty)]
        if (cases[0] < value < cases[1]) is False:
            cases_diff = [np.abs(value - cases[0]),
                          np.abs(value - cases[1])]
            case = cases[cases_diff.index(min(cases_diff))]
        else:
            if self.args.units == 'feet':
                case = self.HighWaterMarks[station]["elev_ft"]
            else:
                case = self.HighWaterMarks[station]["elev_ft"]/3.28084
        self._values["best_case"].append(case)

    def __init_regression_plot(self):
        sns.set(style="darkgrid", color_codes=True)
        sns.set_context("poster")
        colors = ['lightcoral', 'lightsalmon', 'lemonchiffon']
        x = self.dataframe["Observation"]
        y = self.dataframe["ADCIRC"]
        x_max = x.max()
        y_max = y.max()
        middle_line = np.linspace(0, 1.1*np.max([x_max, y_max]),
                                  num=1.1*np.max([x_max, y_max]))
        for i, percentage in enumerate([0.1, 0.2, 0.3]):
            lower = np.linspace(0,
                                1.1*(np.max([x_max, y_max])*(1.-percentage)),
                                num=1.1*np.max([x_max, y_max]))
            if percentage == 0.1:
                _lower = plt.fill_between(middle_line, lower, middle_line,
                                          facecolor=colors[i], alpha=0.5)
                _upper = plt.fill_betweenx(middle_line, lower, middle_line,
                                           facecolor=colors[i], alpha=0.5)
            plt.fill_between(middle_line, lower, middle_line,
                             facecolor=colors[i], alpha=0.5)
            plt.fill_betweenx(middle_line, lower, middle_line,
                              facecolor=colors[i], alpha=0.5)
            plt.plot(middle_line, middle_line, color=colors[i],
                     alpha=0.5)
        lower_path = _lower.get_paths().pop()
        upper_path = _upper.get_paths().pop()
        self._90_percent_index = 0
        for i, obs in enumerate(x):
            if lower_path.contains_point((obs, self.dataframe["ADCIRC"][i])) \
             or upper_path.contains_point((obs, self.dataframe["ADCIRC"][i])) \
             or obs == self.dataframe["ADCIRC"][i]:
                self._90_percent_index += 1
        plt.plot(middle_line, middle_line,
                 linestyle='--', color='grey', alpha=0.75)
        sns.regplot(x, y, fit_reg=False, ci=None)
        plt.errorbar(x, y, xerr=self.dataframe["Uncertainty"],
                     fmt='none', color='darkorange')
        plt.gca().set(xlabel='Observations [{}:NAVD88]'.format(
                                                            self.args.units),
                      ylabel='Hindcast [{}:NAVD88]'.format(self.args.units))
        title = 'High Water Mark Validation for {}, {}'.format(
                                            self.HighWaterMarks.event_name,
                                            self.HighWaterMarks.event_year)
        if self.args.subtitle is not None:
            title += '\n {}'.format(self.args.subtitle)
        plt.title(title)
        plt.axis('scaled')
        plt.axis([0, 1.1*np.max([x.max(), y.max()]),
                  0, 1.1*np.max([x.max(), y.max()])])
        plt.gcf().set_size_inches(10.5, 10.5)
        self.__save_plot('linear_regression.png')
        self.__show_plot()
        plt.close(plt.gcf())

    def __init_regression_plot_best_case(self):
        sns.set(style="darkgrid", color_codes=True)
        sns.set_context("poster")
        colors = ['lightcoral', 'lightsalmon', 'lemonchiffon']
        x = self.dataframe["best_case"]
        y = self.dataframe["ADCIRC"]
        x_max = x.max()
        y_max = y.max()
        middle_line = np.linspace(0, 1.1*np.max([x_max, y_max]), num=1.1*np.max([x_max, y_max]))
        for i, percentage in enumerate([0.1, 0.2, 0.3]):
            lower = np.linspace(0, 1.1*(np.max([x_max, y_max])*(1.-percentage)), num=1.1*np.max([x_max, y_max]))
            if percentage == 0.1:
                _lower = plt.fill_between(middle_line, lower, middle_line, facecolor=colors[i], alpha=0.5)
                _upper = plt.fill_betweenx(middle_line, lower, middle_line, facecolor=colors[i], alpha=0.5)
            plt.fill_between(middle_line, lower, middle_line, facecolor=colors[i], alpha=0.5)
            plt.fill_betweenx(middle_line, lower, middle_line, facecolor=colors[i], alpha=0.5)
            plt.plot(middle_line, middle_line, color=colors[i], alpha=0.5)
        lower_path = _lower.get_paths().pop()
        upper_path = _upper.get_paths().pop()
        self._90_percent_index = 0
        for i, obs in enumerate(x):
            if lower_path.contains_point((obs, self.dataframe["ADCIRC"][i])) or \
                 upper_path.contains_point((obs, self.dataframe["ADCIRC"][i])) or \
                 obs==self.dataframe["ADCIRC"][i]:
                self._90_percent_index+=1
        plt.plot(middle_line, middle_line, linestyle='--', color='grey', alpha=0.75)
        sns.regplot(x, y, fit_reg=False, ci=None)
        plt.gca().set(xlabel='Observations [{}:NAVD88]'.format(self.args.units), ylabel='Hindcast [{}:NAVD88]'.format(self.args.units))
        title = 'High Water Mark Validation for {}, {} (best case)'.format(self.HighWaterMarks.event_name, self.HighWaterMarks.event_year)
        if self.args.subtitle is not None:
            title += '\n {}'.format(self.args.subtitle)
        plt.title(title)
        plt.axis('scaled')
        plt.axis([0, 1.1*np.max([x.max(), y.max()]), 0, 1.1*np.max([x.max(), y.max()])])
        plt.gcf().set_size_inches(10.5, 10.5)
        self.__save_plot('linear_regression_best_case.png')
        self.__show_plot()
        plt.close(plt.gcf())

    def __get_pearson_r(self, data):
        return stats.pearsonr(self.dataframe["ADCIRC"], data)[0]

    def __get_bias(self, data):
        bias = list()
        for i, obs in enumerate(data):
            bias.append(self.dataframe["ADCIRC"][i]-obs)
        return sum(bias)/len(bias)

    def __get_RMS_demeaned(self, data, bias):
        rms_demeaned = list()
        for i, obs in enumerate(data):
            rms_demeaned.append((self.dataframe["ADCIRC"][i]-obs-bias)**2)
        return np.sqrt(sum(rms_demeaned)/(len(rms_demeaned)-1))

    def __get_scatter_index(self, data, rms_demeaned):
        return rms_demeaned/(sum(data)/len(data))

    def __get_90percent_index(self, data):
        return self._90_percent_index/len(data)

    def __show_plot(self):
        if self.args.show is True:
            plt.show()

    def __save_plot(self, filename):
        if self.args.save_path is not None:
            plt.savefig(self.args.save_path+'/'+filename, bbox_inches='tight')


def get_parser():
    parser = argparse.ArgumentParser(description='Program for validating \
                  tidal harmonic constituent outputs against COOPS data.')
    parser.add_argument('event_id', help='ATCF id of the hurricane')
    parser.add_argument('maxele', help='Path to the maxele file.')
    parser.add_argument('navd88_mesh', help='Path to the datum conversion \
                       mesh from mesh datum to NAVD88. Required since HWMs are \
                       in NAVD88 and usually the mesh is in LMSL based on \
                       VDatum.')
    parser.add_argument('--subtitle', help='Custom subtitle that applies to all the plots.')
    parser.add_argument('--units', default='meters', choices=['feet', 'meters'], help='Units for the plots.')
    parser.add_argument('--fort14', help='Path to mesh file. Required if maxele file is ASCII.')
    parser.add_argument('--save-path', '--save', help='Path where to save figures.')
    show = parser.add_mutually_exclusive_group(required=False)
    show.add_argument('--show', dest='show', action='store_true', help='Shows plots to screen as they are generated (default).')
    show.add_argument('--no-show', dest='show', action='store_false', help='Prevents the plots from showing to screen. Useful for only saving the plots without showing them.')
    parser.set_defaults(show=True)
    parser.add_argument('--method', choices=['nearest', 'nanmean', 'idw'], default='nearest', help='Pairing method. Nearest within closing element is default.')
    parser.add_argument('--csv', help='Use CSV file for High Water Marks instead of USGS REST server.')
    return parser


def main():
    HighWaterMarkValidation(get_parser().parse_args())


if __name__ == '__main__':
    main()
