from .colors import *
from .SunSpatial import *

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.tri as tri
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

mpl.rcParams['font.size'] = 7


def draw_pie(dist, xpos, ypos, size, ax=None):
    angles = np.linspace(2*np.pi*0, 2*np.pi*dist)
    x = [0] + np.cos(angles).tolist()
    y = [0] + np.sin(angles).tolist()
    xy = np.column_stack([x, y])
    ax.scatter([xpos], [ypos], marker=xy, s=size, c="#000000", alpha=0.5, edgecolors='none')


def draw_box(dist, xpos, ypos, size, ax=None):
    y0 = 2
    xy = [[1, -y0], [1, y0], [-1, y0], [-1, -y0]]
    ax.scatter([xpos], [ypos], marker=xy, s=size, c="#cccccc", alpha=1, edgecolors='none', zorder=5)
    y = y0*(-1 + 2*dist)
    xy = [[1, -y0], [1, y], [-1, y], [-1, -y0]]
    ax.scatter([xpos], [ypos], marker=xy, s=size, c="#3a3a3a", alpha=1, edgecolors='none', zorder=5)


def draw_boxa(dist, xpos, ypos, size, ax=None):
    a = 1
    xy = [[a, -a], [a, a], [-a, a], [-a, -a]]
    ax.scatter([xpos], [ypos], marker=xy, s=size, c="#000000", alpha=0.2, edgecolors='none')
    a = math.sqrt(dist)
    xy = [[a, -a], [a, a], [-a, a], [-a, -a]]
    ax.scatter([xpos], [ypos], marker=xy, s=dist*size, c="#000000", alpha=0.5, edgecolors='none')


class SunPlotter:
    sunSpatialPlot: SunSpatial
    sunSpatialCalc: SunSpatial
    tris: np.ndarray
    levels: list
    labelPositionsE: list
    labelPositionsH: list

    def __init__(self, sunSpatialPlot: SunSpatial, sunSpatialCalc: SunSpatial, levels=[0., 1., 0.1],
                 labelPositionsE=[], labelPositionsH=[]
                 ):
        self.sunSpatialPlot = sunSpatialPlot
        self.sunSpatialCalc = sunSpatialCalc
        triangulation = tri.Triangulation(
            [e.hourAngle() for e in self.sunSpatialPlot.nodesE],
            [e.declination() for e in self.sunSpatialPlot.nodesE]
        )
        self.tris = triangulation.triangles
        self.levels = levels
        self.labelPositionsE = labelPositionsE
        self.labelPositionsH = labelPositionsH

    def showEquatorial(self, etaF: typing.Callable, **kwargs):
        etaMin, etaMax, etaStep = kwargs.get('levels', self.levels)
        colorMap = kwargs.get('cmap', colorMapL)
        labelsFormat = kwargs.get('labelsFormat', '%0.1f')
        labelPositions = kwargs.get('labelPositions', self.labelPositionsE)
        cbarTitle = kwargs.get('cbarTitle', r'$\eta$')
        showWeights = kwargs.get('showWeights', False)

        pointsPlot = np.array([
            (e.hourAngle(), e.declination()) for e in self.sunSpatialPlot.nodesE
        ])
        xm = pointsPlot[:, 0]/angleHour
        ym = pointsPlot[:, 1]/degree
        zm = np.array([etaF(v) for v in self.sunSpatialPlot.nodesV])

        points = np.array([
            (e.hourAngle(), e.declination()) for e in self.sunSpatialCalc.nodesE
        ])
        snws = np.array(self.sunSpatialCalc.nodesW.copy())
        if len(snws):
            snws = np.fabs(snws)
            snwsmax = np.max(snws)
            snws /= snwsmax
            # snws = 1 -snws
            # snws = snws**2

        levelsAll = np.arange(etaMin, etaMax + 1e-6, etaStep/10)
        ticksLocs = (0.2, etaStep)  # major, minor
        ticksMajor = np.arange(etaMin, etaMax + 1e-6, etaStep)

        # figure
        figure, axes = plt.subplots(constrained_layout=True)
        figure.set_size_inches(3.5, 1.5)
        figure.set_dpi(150)

        # parts
        CS = axes.tricontourf(xm, ym, zm, levels=levelsAll, cmap=colorMap, zorder=-2, linestyles='solid')
        CS2 = axes.tricontour(xm, ym, zm, levels=CS.levels[::1], colors='#00000020', linewidths=0.3, linestyles='solid')
        CS3 = axes.tricontour(xm, ym, zm, levels=CS.levels[::10], colors='#00000050', linewidths=0.3,
                              linestyles='solid')

        axes.clabel(CS3, fmt=labelsFormat, fontsize=7, colors=colorAxes, manual=labelPositions)

        if showWeights and len(snws) > 0:
            for q, w in zip(points, snws):
                draw_box(w, q[0]/angleHour, q[1]/degree, 2*32, ax=axes)
        else:
            axes.scatter(points[:, 0]/angleHour, points[:, 1]/degree, s=[15], c="#000000", alpha=0.5, edgecolors='none',
                         zorder=5)

        # frame
        axes.set_xlabel(r'Hour angle $\omega$, h')
        axes.set_xlim(-8, 8)
        axes.xaxis.set_major_locator(MultipleLocator(2))
        axes.xaxis.set_minor_locator(MultipleLocator(0.5))

        axes.set_ylabel(r'Declination $\delta$, deg')
        axes.set_ylim(-28, 28)
        if showWeights and len(snws) > 0:
            axes.set_ylim(-30, 30)
        axes.yaxis.set_major_locator(MultipleLocator(10))
        axes.yaxis.set_minor_locator(MultipleLocator(2))

        axes.tick_params(which='both', direction='in', top=True, right=True, width=0.5, color=colorAxes)
        axes.tick_params(which='major', length=3)
        axes.tick_params(which='minor', length=1.5)

        for spine in axes.spines.values():
            spine.set_edgecolor(colorAxes)
            spine.set_linewidth(0.5)

        for v in [-6, 0, 6]:
            axes.axvline(x=v, c=colorGrid, ls='--', lw=0.5, zorder=-1)
        for v in [0]:
            axes.axhline(y=v, c=colorGrid, ls='--', lw=0.5, zorder=-1)

        # colorbar
        cbar = SunPlotter.addColorBar(figure, CS, 20, ticksLocs, ticksMajor, cbarTitle)
        cbar.add_lines(CS3)

        # show
        # plt.show()
        # figure.savefig("plot.png", dpi=300)
        return figure

    def showHorizontal(self, etaF: typing.Callable, **kwargs):
        etaMin, etaMax, etaStep = kwargs.get('levels', self.levels)
        colorMap = kwargs.get('cmap', colorMapL)
        labelsFormat = kwargs.get('labelsFormat', '%0.1f')
        labelPositions = kwargs.get('labelPositions', self.labelPositionsH)
        cbarTitle = kwargs.get('cbarTitle', r'$\eta$')
        showWeights = kwargs.get('showWeights', False)

        pointsPlot = np.array([
            (h.azimuth(), h.elevation()) for h in self.sunSpatialPlot.nodesH
        ])
        azs = pointsPlot[:, 0]/degree
        els = pointsPlot[:, 1]/degree

        zm = np.array([etaF(v) for v in self.sunSpatialPlot.nodesV])

        points = np.array([
            (h.azimuth(), h.elevation()) for h in self.sunSpatialCalc.nodesH
        ])
        snws = np.array(self.sunSpatialCalc.nodesW.copy())
        if len(snws):
            snws = np.fabs(snws)
            snwsmax = np.max(snws)
            snws /= snwsmax
            # snws = 1 -snws
            # snws = snws**2

        levelsAll = np.arange(etaMin, etaMax + 1e-6, etaStep/10)
        ticksLocs = (0.2, etaStep)  # major, minor
        ticksMajor = np.arange(etaMin, etaMax + 1e-6, etaStep)

        # figure
        figure, axes = plt.subplots(constrained_layout=True)

        figure.set_size_inches(3.5, 1.5)
        figure.set_dpi(150)

        # parts
        CS = axes.tricontourf(azs, els, self.tris, zm, levels=levelsAll, cmap=colorMap, zorder=-2, linestyles='solid')
        CS2 = axes.tricontour(azs, els, self.tris, zm, levels=CS.levels[::1], colors='#00000020', linewidths=0.3,
                              linestyles='solid')
        CS3 = axes.tricontour(azs, els, self.tris, zm, levels=CS.levels[::10], colors='#00000050', linewidths=0.3,
                              linestyles='solid')

        axes.clabel(CS3, fmt=labelsFormat, fontsize=7, colors=colorAxes, manual=labelPositions)

        if showWeights and len(snws) > 0:
            for q, w in zip(points, snws):
                draw_box(w, q[0]/angleHour, q[1]/degree, 20, ax=axes)
        else:
            axes.scatter(points[:, 0]/degree, points[:, 1]/degree, s=[15], c="#000000", alpha=0.5, edgecolors='none',
                         zorder=5)

        # frame
        axes.set_xlabel(r'Azimuth $\gamma$, deg')
        axes.set_xlim(50, 360 - 50)
        axes.xaxis.set_major_locator(MultipleLocator(30))
        axes.xaxis.set_minor_locator(MultipleLocator(10))

        axes.set_ylabel(r'Elevation $\alpha$, deg')
        axes.set_ylim(-5, 95)
        axes.yaxis.set_major_locator(MultipleLocator(30))
        axes.yaxis.set_minor_locator(MultipleLocator(10))

        axes.tick_params(which='both', direction='in', top=True, right=True, width=0.5, color=colorAxes)
        axes.tick_params(which='major', length=3)
        axes.tick_params(which='minor', length=1.5)

        for spine in axes.spines.values():
            spine.set_edgecolor(colorAxes)
            spine.set_linewidth(0.5)

        colorD = '#B0B0B0FF'
        for v in np.arange(30, 360, 30):
            axes.axvline(x=v, c=colorD, ls='-', lw=0.25, zorder=-5)
        for v in np.arange(0, 90.1, 30):
            axes.axhline(y=v, c=colorD, ls='-', lw=0.25, zorder=-5)

        # colorbar
        cbar = SunPlotter.addColorBar(figure, CS, 20, ticksLocs, ticksMajor, cbarTitle)
        cbar.add_lines(CS3)

        return figure

    @staticmethod
    def addColorBar(figure, CS, aspect, ticksLocs, ticksMajor, cbarTitle):
        cbar = figure.colorbar(CS, aspect=aspect)
        cbar.ax.set_xlabel(cbarTitle, labelpad=7)
        # cbar.ax.set_xlabel(r'${\rm W}/{\rm m}^2$', labelpad=7)
        # cbar.ax.set_ylabel('$\eta$', rotation=0)

        # cbar.ax.yaxis.set_major_locator(MultipleLocator(ticksLocs[0]))
        cbar.ax.yaxis.set_minor_locator(MultipleLocator(ticksLocs[1]))
        cbar.set_ticks(ticksMajor)

        cbar.ax.tick_params(which='both', direction='in', top=True, right=True, width=0.5, color=colorAxes)
        cbar.ax.tick_params(which='major', length=3)
        cbar.ax.tick_params(which='minor', length=1.5)

        cbar.outline.set_linewidth(0.5)
        for spine in cbar.ax.spines.values():
            spine.set_edgecolor(colorAxes)
            spine.set_linewidth(0.5)

        return cbar
