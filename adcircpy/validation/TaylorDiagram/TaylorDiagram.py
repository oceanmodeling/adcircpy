"""
Taylor diagram (Taylor, 2001) test implementation.
http://www-pcmdi.llnl.gov/about/staff/Taylor/CV/Taylor_diagram_primer.htm
"""

__version__ = "Time-stamp: <2012-02-17 20:59:35 ycopin>"
__author__ = "Yannick Copin <yannick.copin@laposte.net>"

# Rev. 1.0  saeed Moghimi  moghimis@gmail.com

import unittest

import matplotlib.pyplot as plt
import numpy as np


class TaylorDiagram:
    """Taylor diagram: plot model standard deviation and correlation
    to reference (data) sample in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=111, label='_'):
        """Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using mpl_toolkits.axisartist.floating_axes. refstd is
        the reference standard deviation to be compared to.
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd  # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.concatenate((np.arange(10) / 10., [0.95, 0.99]))
        tlocs = np.arccos(rlocs)  # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)  # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent
        self.smin = 0
        self.smax = 1.5 * self.refstd

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(
                                               0, np.pi / 2,  # 1stquadrant
                                               self.smin, self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           )

        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")

        ax.axis["right"].set_axis_direction("top")  # "Y axis"
        ax.axis["right"].label.set_text("Standard deviation")

        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")

        ax.axis["bottom"].set_visible(False)  # Useless

        # Contours along standard deviations
        ax.grid(False)

        self._ax = ax  # Graphical axes
        self.ax = ax.get_aux_axes(tr)  # Polar coordinates

        # Add reference point and stddev contour
        # print("Reference std:", self.refstd)
        # l, = self.ax.plot([0], self.refstd, 'k*',     rem by saeed
        #                  ls='', ms=10, label=label)  rem by saeed
        # t = np.linspace(0, np.pi/2)                   rem by saeed
        # r = np.zeros_like(t) + self.refstd            rem by saeed
        # self.ax.plot(t,r, 'k--', label='_')           rem by saeed

        for ip in range(len(rlocs)):
            x = [self.smin, self.smax]
            y = [np.arccos(rlocs[ip]), np.arccos(rlocs[ip])]
            self.ax.plot(y, x, 'grey', linewidth=0.25)

        # Collect sample points for latter use (e.g. legend)
        # self.samplePoints = [l]                       rem by saeed
        self.samplePoints = []  # add by saeed

    def add_sample(self, stddev, corrcoef, ref, *args, **kwargs):
        """Add sample (stddev,corrcoeff) to the Taylor diagram. args
        and kwargs are directly propagated to the Figure.plot
        command."""

        l, = self.ax.plot(np.arccos(corrcoef), stddev,
                          *args, **kwargs)  # (theta,radius)
        self.samplePoints.append(l)

        if ref:
            t = np.linspace(0, np.pi / 2)  # add by saeed
            r = np.zeros_like(t) + stddev  # add by saeed
            self.ax.plot(t, r, 'grey', linewidth=0.25,
                         label='_')  # add by saeed
        return l

    def add_contours(self, levels, data_std, **kwargs):
        """Add constant centered RMS difference contours."""
        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax),
                             np.linspace(0, np.pi / 2))
        # Compute centered RMS difference
        rms = np.sqrt(data_std ** 2 + rs ** 2 - 2 * data_std * rs * np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        return contours


# Unittest
class TaylorDiagramTestCase(unittest.TestCase):
    def test_TaylorDiagram(self):
        # Reference dataset
        x = np.linspace(0, 4 * np.pi, 100)
        data = np.sin(x)
        refstd = data.std(ddof=1)  # Reference standard deviation

        # Models
        m1 = data + 0.2 * np.random.randn(len(x))  # Model 1
        m2 = 0.8 * data + .1 * np.random.randn(len(x))  # Model 2
        m3 = np.sin(x - np.pi / 10)  # Model 3

        # Compute stddev and correlation coefficient of models
        samples = np.array([[m.std(ddof=1), np.corrcoef(data, m)[0, 1]]
                            for m in (m1, m2, m3)])

        fig = plt.figure(1, figsize=(10, 4))
        fig.clf()
        ax1 = fig.add_subplot(1, 2, 1, xlabel='X', ylabel='Y')
        # Taylor diagram
        dia = TaylorDiagram(refstd, fig=fig, rect=122, label="Reference")
        colors = plt.matplotlib.cm.jet(np.linspace(0, 1, len(samples)))
        ax1.plot(x, data, 'ko', label='Data')
        for i, m in enumerate([m1, m2, m3]):
            ax1.plot(x, m, c=colors[i], label='Model %d' % (i + 1))
        ax1.legend(numpoints=1, prop=dict(size='small'), loc='best')
        # Add samples to Taylor diagram
        for i, (stddev, corrcoef) in enumerate(samples):
            dia.add_sample(stddev, corrcoef, ref=False,
                           marker='s', ls='', c=colors[i],
                           label="Model %d" % (i + 1))

        dia.add_sample(refstd, 1.0, ref=True, marker='*', ls='',
                       c='k', label="Ref.")

        # Add RMS contours, and label them
        contours = dia.add_contours(levels=5, data_std=refstd, colors='0.5')
        plt.clabel(contours, inline=1, fontsize=10)

        # Add a figure legend
        fig.legend(dia.samplePoints,
                   [p.get_label() for p in dia.samplePoints],
                   numpoints=1, prop=dict(size='small'), loc='upper right')
        plt.show(block=False)
