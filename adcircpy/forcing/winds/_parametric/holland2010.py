import warnings

import matplotlib.pyplot as plt
import numpy as np
from pyschism.forcing.winds.atcf import Bdeck
from scipy import optimize


def holland_B(hurdat, coriolis=True):
    # air_density = 1.225
    air_density = 1.15

    def with_coriolis(Vmax, Rmax, Pn, Pc, eye_lat):
        f = 2.0 * 7.2921e-5 * np.sin(np.radians(np.abs(eye_lat)))
        return (Vmax ** 2 + Vmax * Rmax * f * air_density * np.exp(1)) / (Pn - Pc)

    def no_coriolis(Vmax, Pn, Pc):
        return (Vmax ** 2 * air_density * np.exp(1)) / (Pn - Pc)

    for data in hurdat.values():

        Pc = data['central_pressure']
        Pn = data['background_pressure']
        # avoid negative Holland B parameter as initial guess
        if Pn <= Pc:
            Pn = Pc + 1.0
        if coriolis:
            return with_coriolis(
                data['max_sustained_wind_speed'],
                data['radius_of_maximum_winds'],
                Pn,
                Pc,
                data['eye']['lat'],
            )
        else:
            return no_coriolis(data['max_sustained_wind_speed'], Pn, Pc)


def main():
    storm_id = 'AL152017'
    # storm_id = 'AL182012'
    hurdat = Bdeck(storm_id).data
    for time, data in hurdat.items():
        if len(data['isotachs'].keys()) != 4:
            continue
        # initial guesses
        Vmax = data['max_sustained_wind_speed']
        Rmax = data['radius_of_maximum_winds']
        # print(data)
        # exit()
        x = 1.0
        B = 1.0

        def holland2010(r, B, x):
            return Vmax * (((Rmax / r) ** B) * np.exp(1 - (Rmax / r) ** B)) ** x

        def V(B, x):
            def v(r):
                return holland2010(r, B, x)

            return v

        # B = holland_B(hurdat)
        for quad, isotachs in data['isotachs'].items():
            xdata = []
            ydata = []
            for y, x in isotachs.items():
                xdata.append(x)
                ydata.append(y)
            # xdata.append(Rmax)
            # ydata.append(Vmax)
            # add bounds
            bi = np.finfo(float).eps  # avoid divide by zero
            bf = data['radius_of_last_closed_isobar']
            bounds = (bi, bf)
            p0 = [B, x]
            # do curve fitting
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                popt, pcov = optimize.curve_fit(
                    holland2010,
                    xdata,
                    ydata,
                    p0=p0,
                    # bounds=bounds,
                    method='dogbox',
                )
                print(popt)
            v = V(*popt)
            radii = np.linspace(bi, bf, num=500)
            res = []
            for i in radii:
                res.append(v(i))
            results = np.array(res)
            plt.plot(radii, results, label=quad)
            # plt.gca().axis('scaled')
        plt.legend()
        plt.show()
        plt.close(plt.gcf())
        # Vmax
        # res = minimize(holland2010, x0=[Vmax, Rmax, B, 0.1] )
        # print(res)


def init():
    if __name__ == '__main__':
        try:
            import colored_traceback

            colored_traceback.add_hook(always=True)
        except ModuleNotFoundError:
            pass
        main()


init()
