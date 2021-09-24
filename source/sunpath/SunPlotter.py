from .colors import *
from .SunSpatial import *

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.tri as tri
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

mpl.rcParams['font.size'] = 7


def draw_box(f, x, y, size, ax=None):
    w = 1
    h = 2
    box_full = [[w, -h], [w, h], [-w, h], [-w, -h]]
    ax.scatter([x], [y], marker=box_full, s=size, c="#cccccc", alpha=1, edgecolors='none', zorder=5)
    t = h*(-1 + 2*f)
    box_part = [[w, -h], [w, t], [-w, t], [-w, -h]]
    ax.scatter([x], [y], marker=box_part, s=size, c="#3a3a3a", alpha=1, edgecolors='none', zorder=5)


class SunPlotter:
    sunSpatialPlot: SunSpatial
    sunSpatialCalc: SunSpatial
    tris: np.ndarray
    levels: list
    labelPositionsE: list
    labelPositionsH: list

    def __init__(self, sunSpatialPlot: SunSpatial, sunSpatialCalc: SunSpatial, levels=[0., 1., 0.1, 10], 
                 labelPositionsE=[], labelPositionsH=[]
                 ):
        self.sunSpatialPlot = sunSpatialPlot
        self.sunSpatialCalc = sunSpatialCalc
        triangulation = tri.Triangulation(
            [e.hourAngle() for e in self.sunSpatialPlot.nodesE],
            [e.declination() for e in self.sunSpatialPlot.nodesE]
        )
        self.tris = triangulation.triangles
        self.levels = levels # start, finish, step, substeps
        self.labelPositionsE = labelPositionsE
        self.labelPositionsH = labelPositionsH

    def showEquatorial(self, etaF: typing.Callable, **kwargs):
        etaMin, etaMax, etaStep, etaSubsteps = kwargs.get('levels', self.levels)
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
            (e.hourAngle()/angleHour, e.declination()/degree) for e in self.sunSpatialCalc.nodesE
        ])
        values = np.array(self.sunSpatialCalc.nodesW.copy())
        if showWeights and len(values) > 0:
            values = np.fabs(values)
            vMax = np.max(values)
            values /= vMax
        else:
            showWeights = False

        etaSubstep = etaStep/etaSubsteps
        iMin = math.floor(etaMin/etaSubstep)    
        iMax = math.ceil(etaMax/etaSubstep)  
        levelsAll = etaSubstep*np.arange(iMin, iMax + 1)
        iMin = math.ceil(etaMin/etaStep)    
        iMax = math.floor(etaMax/etaStep) 
        levelsMajor = etaStep*np.arange(iMin, iMax + 1)
        ticksMajor = np.arange(etaMin, etaMax + 1e-6, etaStep)

        # figure
        figure, axes = plt.subplots(constrained_layout=True)
        figure.set_size_inches(3.5, 1.5)
        figure.set_dpi(150)

        # parts
        CS = axes.tricontourf(xm, ym, zm, levels=levelsAll, cmap=colorMap, zorder=-2, linestyles='solid')
        CS2 = axes.tricontour(xm, ym, zm, levels=CS.levels[::1], colors='#00000020', linewidths=0.3, linestyles='solid')
        CS3 = axes.tricontour(xm, ym, zm, levels=levelsMajor, colors='#00000050', linewidths=0.3,
                              linestyles='solid')

        axes.clabel(CS3, fmt=labelsFormat, fontsize=7, colors=colorAxes, manual=labelPositions)

        if showWeights:
            for q, w in zip(points, values):
                draw_box(w, q[0], q[1], 2*32, ax=axes)
        else:
            axes.scatter(points[:, 0], points[:, 1], s=[15], c="#000000", alpha=0.5, edgecolors='none', zorder=5)

        # frame
        axes.set_xlabel(r'Hour angle $\omega$, h')
        axes.set_xlim(-8, 8)
        axes.xaxis.set_major_locator(MultipleLocator(2))
        axes.xaxis.set_minor_locator(MultipleLocator(0.5))

        axes.set_ylabel(r'Declination $\delta$, deg')
        axes.set_ylim(-28, 28)
        if showWeights:
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
        cbar = SunPlotter.addColorBar(figure, CS, 20, (), levelsMajor, cbarTitle)
        cbar.add_lines(CS3)

        return figure

    def showHorizontal(self, etaF: typing.Callable, **kwargs):
        etaMin, etaMax, etaStep, etaSubsteps = kwargs.get('levels', self.levels)
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
            (h.azimuth()/degree, h.elevation()/degree) for h in self.sunSpatialCalc.nodesH
        ])
        values = np.array(self.sunSpatialCalc.nodesW.copy())
        if showWeights and len(values) > 0:
            values = np.fabs(values)
            vMax = np.max(values)
            values /= vMax
        else:
            showWeights = False

        etaSubstep = etaStep/etaSubsteps
        iMin = math.floor(etaMin/etaSubstep)    
        iMax = math.ceil(etaMax/etaSubstep)  
        levelsAll = etaSubstep*np.arange(iMin, iMax + 1)
        iMin = math.ceil(etaMin/etaStep)    
        iMax = math.floor(etaMax/etaStep) 
        levelsMajor = etaStep*np.arange(iMin, iMax + 1)
        ticksMajor = np.arange(etaMin, etaMax + 1e-6, etaStep)

        # figure
        figure, axes = plt.subplots(constrained_layout=True)

        figure.set_size_inches(3.5, 1.5)
        figure.set_dpi(150)

        # parts
        CS = axes.tricontourf(azs, els, self.tris, zm, levels=levelsAll, cmap=colorMap, zorder=-2, linestyles='solid')
        CS2 = axes.tricontour(azs, els, self.tris, zm, levels=CS.levels[::1], colors='#00000020', linewidths=0.3,
                              linestyles='solid')
        CS3 = axes.tricontour(azs, els, self.tris, zm, levels=levelsMajor, colors='#00000050', linewidths=0.3,
                              linestyles='solid')

        axes.clabel(CS3, fmt=labelsFormat, fontsize=7, colors=colorAxes, manual=labelPositions)

        if showWeights:
            for q, w in zip(points, values):
                draw_box(w, q[0], q[1], 2*32, ax=axes)
        else:
            axes.scatter(points[:, 0], points[:, 1], s=[15], c="#000000", alpha=0.5, edgecolors='none', zorder=5)

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
        cbar = SunPlotter.addColorBar(figure, CS, 20, (), levelsMajor, cbarTitle)
        cbar.add_lines(CS3)

        return figure

    @staticmethod
    def addColorBar(figure, CS, aspect, ticksLocs, ticksMajor, cbarTitle):
        cbar = figure.colorbar(CS, aspect=aspect)
        cbar.ax.set_xlabel(cbarTitle, labelpad=7)
        # cbar.ax.set_xlabel(r'${\rm W}/{\rm m}^2$', labelpad=7)
        # cbar.ax.set_ylabel('$\eta$', rotation=0)

        # cbar.ax.yaxis.set_major_locator(MultipleLocator(ticksLocs[0]))
        # cbar.ax.yaxis.set_minor_locator(MultipleLocator(ticksLocs[1]))
        cbar.set_ticks(ticksMajor)

        cbar.ax.tick_params(which='both', direction='in', top=True, right=True, width=0.5, color=colorAxes)
        cbar.ax.tick_params(which='major', length=3)
        cbar.ax.tick_params(which='minor', length=1.5)

        cbar.outline.set_linewidth(0.5)
        for spine in cbar.ax.spines.values():
            spine.set_edgecolor(colorAxes)
            spine.set_linewidth(0.5)

        return cbar
