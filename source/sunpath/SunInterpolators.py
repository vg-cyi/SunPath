from .SunSpatial import *

from scipy import interpolate
from scipy.special import xlogy
from enum import Enum

import csv
import urllib.request
from io import StringIO


def readEtaRefCSV(filename: str) -> typing.Callable:
    """
    read simulation data from CSV
    :param filename:
    :return:
    """
    # read
    if filename.startswith("http"):
        response = urllib.request.urlopen(filename)
        file = StringIO(response.read().decode('utf-8'))
    else:
        file = open(filename, "r")

    reader = csv.reader(file, delimiter=',')  # quoting=csv.QUOTE_NONNUMERIC
    data = []
    for row in reader:
        data += [row]
    data = np.array(data)    

    # convert
    iAzimuth = np.where(data[0] == 'azimuth')[0][0]
    iElevation = np.where(data[0] == 'elevation')[0][0]
    iEfficiency = np.where(data[0] == 'efficiency')[0][0]
    assert iAzimuth >= 0
    assert iElevation >= 0
    assert iEfficiency >= 0
    
    ans = []
    for q in data[1:]:
        azimuth = float(q[iAzimuth])
        elevation = float(q[iElevation])
        eta = float(q[iEfficiency])
        ans.append([azimuth, elevation, eta])
    ans = np.array(ans)

    # interpolator
    gammas = np.unique(ans[:, 0])*degree
    alphas = np.unique(ans[:, 1])*degree
    etas = np.reshape(ans[:, 2], (len(gammas), len(alphas)))

    etaFuncTemp = interpolate.RectBivariateSpline(gammas, alphas, etas)

    def etaFunc(v: mu.Vector) -> float:
        h = SunCalculator.findHorizontalFromVector(v)
        return etaFuncTemp(h.azimuth(), h.elevation())[0, 0]

    return etaFunc

    
def readSolarPilot(filename: str) -> typing.Callable:
    """
    read simulation data from SolarPilot
    :param filename:
    :return:
    """
    # read
    if filename.startswith("http"):
        response = urllib.request.urlopen(filename)
        file = StringIO(response.read().decode('utf-8'))
    else:
        file = open(filename, "r")

    reader = csv.reader(file, delimiter=',')  # quoting=csv.QUOTE_NONNUMERIC
    data = []
    for row in reader:
        data += [row]
    data = np.array(data).T

    # convert
    iPowerIn = np.where(data[0] == 'Power incident on field')[0][0]
    iPowerOut = np.where(data[0] == 'Power absorbed by the receiver')[0][0]
    iAzimuth = np.where(data[0] == 'fluxsim: Solar azimuth angle (0=N)')[0][0]
    iElevation = np.where(data[0] == 'fluxsim: Solar elevation angle')[0][0]
    iAbsorption = np.where(data[0] == 'Absorption efficiency')[0][0]
    assert iPowerIn >= 0
    assert iPowerOut >= 0
    assert iAzimuth >= 0
    assert iElevation >= 0
    assert iAbsorption >= 0

    ans = []
    for q in data[2:-1]:
        azimuth = float(q[iAzimuth])
        elevation = float(q[iElevation])
        powerIn = float(q[iPowerIn])
        powerOut = float(q[iPowerOut])
        etaAbsorption = float(q[iAbsorption])/100
        eta = powerOut/etaAbsorption/powerIn
        ans.append([azimuth, elevation, eta])
    ans = np.array(ans)

    # interpolator
    gammas = np.unique(ans[:, 0])*degree
    alphas = np.unique(ans[:, 1])*degree
    etas = np.reshape(ans[:, 2], (len(gammas), len(alphas)))

    etaFuncTemp = interpolate.RectBivariateSpline(gammas, alphas, etas)

    def etaFunc(v: mu.Vector) -> float:
        h = SunCalculator.findHorizontalFromVector(v)
        return etaFuncTemp(h.azimuth(), h.elevation())[0, 0]

    return etaFunc


class SphericalRBF:
    nodes: list
    amplitudes: list
    sigma: float

    def kernel(self, ra: mu.Vector, rb: mu.Vector) -> float:
        arg = ra.dot(rb) - 1.
        arg /= self.sigma ** 2
        return math.exp(arg)

    def __init__(self, nodes: list, values: list, sigma: float = 20.*degree):
        self.sigma = sigma
        self.nodes = nodes
        matrix = np.array([[
            self.kernel(p, q)
            for q in self.nodes]
            for p in self.nodes])
        self.amplitudes = np.linalg.solve(matrix, values)
        # self.amplitudes = np.linalg.lstsq(matrix, values, rcond=None)[0]

    def __call__(self, v: mu.Vector) -> float:
        temp = [a*self.kernel(vn, v) for vn, a in zip(self.nodes, self.amplitudes)]
        return math.fsum(temp)


class InterpolationMethod:
    """
    method: str
        RBF-2D,
        RBF-3D
        triangular,
        B-spline

    kernel (for RBF-2D): str or callable
        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'thin_plate': r**2 * log(r)
        'cubic': r**3
        'quintic': r**5

    kernel (for RBF-3D): str or callable
        'polyharmonic':
        'gaussian'
        'gaussian-scipy'

    order:
        for triangular, polyharmonic and B-splines

    sigma: float
        shape parameter

    precondition:

    node: int
        for localized kernels
    """
    def __init__(self, method:str = 'RBF-3D', kernel='polyharmonic',
                 order=5, sigma=20*degree,
                 precondition=lambda s: 1., node=-1):
        self.method = method
        self.kernel = kernel
        self.order = order
        self.sigma = sigma
        self.precondition = precondition
        self.node = node

    def findWeights(self, sunTemporal: SunTemporal, sunSpatial: SunSpatial, etaRef: typing.Callable):
        overlaps = []
        n0 = self.node
        for n in range(len(sunSpatial.nodesV)):
            self.node = n
            etaL = MakeInterpolation(sunSpatial, etaRef, self)
            overlap = sunTemporal.integrate(etaL)
            overlaps.append(overlap)
        self.node = n0
        sunSpatial.nodesW = overlaps

    def checkDiff(self, sunTemporal: SunTemporal, sunSpatial: SunSpatial, etaRef: typing.Callable):
        integralRef = sunTemporal.integrate(etaRef)
        integralW = sunTemporal.integrate(lambda x: 1.)
        etaInterp = MakeInterpolation(sunSpatial, etaRef, self)
        integral = sunSpatial.integrate(etaInterp)
        return (integral - integralRef)/integralW


class InterpolationFormat(Enum):
    f_xy = 1
    f_xy_bound = 2
    f_xy_grid = 3
    f_xyz = 4
    f_v = 5


class MakeInterpolation:
    sc: SunCalculator
    method: InterpolationMethod
    interpolationFormat: InterpolationFormat
    etaInterpolated: typing.Callable

    def __init__(self, sunSpatial: SunSpatial, etaRef: typing.Callable,
                 method: InterpolationMethod = InterpolationMethod()):
        self.sc = sunSpatial.sc
        self.method = method
        im = self.method

        points = np.array([
            (e.hourAngle(), e.declination()) for e in sunSpatial.nodesE
        ])
        if im.node < 0:
            values = [etaRef(v)*im.precondition(v) for v in sunSpatial.nodesV]
        else:
            values = np.zeros(len(sunSpatial.nodesV))
            values[im.node] = im.precondition(sunSpatial.nodesV[im.node])

        if im.method == 'triangular' and im.order == 1:
            self.etaInterpolated = interpolate.LinearNDInterpolator(
                points, values, fill_value=0.
            )
            self.interpolationFormat = InterpolationFormat.f_xy_bound
        elif im.method == 'triangular' and im.order == 2:
            self.etaInterpolated = interpolate.CloughTocher2DInterpolator(
                points, values, fill_value=0.
            )
            self.interpolationFormat = InterpolationFormat.f_xy_bound
        elif im.method == 'RBF-2D':
            self.etaInterpolated = interpolate.Rbf(
                points[:, 0], points[:, 1], values, function=im.kernel
            )
            self.interpolationFormat = InterpolationFormat.f_xy
        elif im.method == 'RBF-3D':
            points = np.array([
                (v.x, v.y, v.z) for v in sunSpatial.nodesV
            ])
            if isinstance(im.kernel, str):
                if im.kernel == 'polyharmonic':
                    k = im.order
                    if k%2 == 0:
                        kernelt = lambda r: xlogy(r ** k, r)
                    else:
                        kernelt = lambda r: r ** k
                    self.etaInterpolated = interpolate.Rbf(
                        points[:, 0], points[:, 1], points[:, 2], values, function=kernelt
                    )
                    self.interpolationFormat = InterpolationFormat.f_xyz
                    return
                elif im.kernel == 'gaussian-scipy':
                    self.etaInterpolated = interpolate.Rbf(
                        points[:, 0], points[:, 1], points[:, 2], values,
                        norm=lambda ra, rb: np.dot(ra, rb) - 1.,
                        function=lambda s, r: np.exp(r/s.epsilon**2),
                        epsilon=im.sigma
                    )
                    self.interpolationFormat = InterpolationFormat.f_xyz
                    return
                elif im.kernel == 'gaussian':
                    self.etaInterpolated = SphericalRBF(sunSpatial.nodesV, values, im.sigma)
                    self.interpolationFormat = InterpolationFormat.f_v
                    return
            else:
                self.etaInterpolated = interpolate.Rbf(
                    points[:, 0], points[:, 1], points[:, 2], values, function=im.kernel
                )
                self.interpolationFormat = InterpolationFormat.f_xyz
                return
        elif im.method == 'B-spline':
            assert len(sunSpatial.dims) == 2
            omegas = np.linspace(-6*angleHour, 6*angleHour, sunSpatial.dims[0])
            deltas = np.linspace(-self.sc.obliquity, self.sc.obliquity, sunSpatial.dims[1])
            zs = np.reshape(np.array(values), (sunSpatial.dims[1], sunSpatial.dims[0])).T
            ks = im.order
            self.etaInterpolated = interpolate.RectBivariateSpline(omegas, deltas, zs, kx=ks[0], ky=ks[1])
            self.interpolationFormat = InterpolationFormat.f_xy_grid

    def __call__(self, v: mu.Vector) -> float:
        im = self.method
        pc = im.precondition(v)
        if self.interpolationFormat == InterpolationFormat.f_v:
            return self.etaInterpolated(v)/pc
        elif self.interpolationFormat == InterpolationFormat.f_xyz:
            return self.etaInterpolated(v.x, v.y, v.z)/pc
        elif self.interpolationFormat == InterpolationFormat.f_xy:
            e = self.sc.findEquatorialFromVector(v)
            return self.etaInterpolated(e.hourAngle(), e.declination())/pc
        elif self.interpolationFormat == InterpolationFormat.f_xy_bound:
            e = self.sc.findEquatorialFromVector(v)
            declination = np.clip(e.declination(), -self.sc.obliquity, self.sc.obliquity)
            return self.etaInterpolated(e.hourAngle(), declination)/pc
        elif self.interpolationFormat == InterpolationFormat.f_xy_grid:
            e = self.sc.findEquatorialFromVector(v)
            declination = np.clip(e.declination(), -self.sc.obliquity, self.sc.obliquity)
            s = 6*angleHour/self.sc.findHourAngleMax(declination)
            return self.etaInterpolated(s*e.hourAngle(), declination)[0][0]/pc

# Inverse Distance Weighted
# https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python/3119544#3119544
