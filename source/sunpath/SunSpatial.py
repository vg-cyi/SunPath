from .SunTemporal import *

def sampleInterval(a: float, b: float, rho: float, middle: bool = False) -> np.ndarray:
    """
    :param a: interval begin
    :param b: interval end
    :param rho: resolution
    :param middle: make point in the middle
    :return: points
    """
    if middle:
        n = 2*round((b - a)/(2*rho))
    else:
        n = round((b - a)/rho)
    return np.linspace(a, b, n + 1)


class SunSpatial:
    sc: SunCalculator
    nodesV: list  # vectors
    nodesH: list  # horizontal
    nodesE: list  # equatorial
    nodesW: list  # weights
    dims: list  # dimensions, used only for grid

    def __init__(self, sc: SunCalculator):
        self.sc = sc
        self.nodesW = []

    def sampleEquatorial(self, rho: float, grid: bool = False) -> None:
        """
        sample sun path in equatorial coordinates with a resolution
        :param rho: resolution
        :param grid: make stretched rectangular grid
        :return:
        """
        sc = self.sc

        ans = []
        deltas = sampleInterval(-sc.obliquity, sc.obliquity, rho)
        if grid:
            omegas = sampleInterval(-6.*angleHour, 6.*angleHour, rho, middle=True)
            self.dims = [len(omegas), len(deltas)]
            for delta in deltas:
                s = sc.findHourAngleMax(delta)/(6.*angleHour)
                for omega in omegas:
                    ans += [PointEquatorial(omega*s, delta)]
        else:
            self.dims = []
            for delta in deltas:
                omegaMax = sc.findHourAngleMax(delta)
                omegas = sampleInterval(-omegaMax, omegaMax, rho, middle=True)
                for omega in omegas:
                    ans += [PointEquatorial(omega, delta)]

        self.nodesE = ans
        self.nodesV = [sc.findVectorFromEquatorial(e) for e in self.nodesE]
        self.nodesH = [sc.findHorizontalFromVector(v) for v in self.nodesV]

    def info(self):
        n = len(self.nodesE)
        print('points: {:n}'.format(n))

    def checkAccuracy(self, etaF: typing.Callable, etaF_ref: typing.Callable) -> None:
        """
        quick check of interpolation accuracy
        :param etaF: interpolated function
        :param etaF_ref: reference function
        :return:
        """
        pointsF = np.array([etaF(v) for v in self.nodesV])
        pointsF_ref = np.array([etaF_ref(v) for v in self.nodesV])
        delta_rms = np.sqrt(np.mean((pointsF - pointsF_ref) ** 2))
        delta_m = np.mean(pointsF - pointsF_ref)
        print('delta_rms = {:.3f}%'.format(delta_rms*100))
        print('delta_m = {:.3f}%'.format(delta_m*100))

    def integrate(self, values) -> float:
        """
        find annual integral
        :param values: values at sampling points or a function
        :return: integral
        """
        if isinstance(values, list):
            vs = values
        elif callable(values):
            vs = [values(v) for v in self.nodesV]
        temp = [v*o for v, o in zip(vs, self.nodesW)]
        return sum(temp)