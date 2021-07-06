from .coordinates import *

import numpy as np
import typing


class SunCalculator:
    obliquity: float
    location: Location

    def __init__(self, location: Location = Location()):
        self.obliquity = 23.4*degree
        self.location = location

    @staticmethod
    def findVectorFromHorizontal(h: PointHorizontal) -> mu.Vector:
        gamma = h.azimuth()
        cosGamma = math.cos(gamma)
        sinGamma = math.sin(gamma)
        alpha = h.elevation()
        cosAlpha = math.cos(alpha)
        sinAlpha = math.sin(alpha)
        return mu.Vector((
            sinGamma*cosAlpha,
            cosGamma*cosAlpha,
            sinAlpha
        ))

    @staticmethod
    def findHorizontalFromVector(v: mu.Vector) -> PointHorizontal:
        gamma = math.atan2(v.x, v.y)
        alpha = math.asin(v.z)
        if gamma < 0.:
            gamma += 2.*math.pi
        return PointHorizontal(gamma, alpha)

    def findVectorFromEquatorial(self, e: PointEquatorial) -> mu.Vector:
        phi = self.location.point.latitude()
        cosPhi = math.cos(phi)
        sinPhi = math.sin(phi)
        omega = e.hourAngle()
        cosOmega = math.cos(omega)
        sinOmega = math.sin(omega)
        delta = e.declination()
        cosDelta = math.cos(delta)
        sinDelta = math.sin(delta)
        return mu.Vector((
            -sinOmega*cosDelta,
            cosPhi*sinDelta - sinPhi*cosOmega*cosDelta,
            sinPhi*sinDelta + cosPhi*cosOmega*cosDelta
        ))

    def findEquatorialFromVector(self, v: mu.Vector) -> PointEquatorial:
        phi = self.location.point.latitude()
        cosPhi = math.cos(phi)
        sinPhi = math.sin(phi)
        omega = math.atan2(-v.x, -v.y*sinPhi + v.z*cosPhi)
        temp = np.clip(v.y*cosPhi + v.z*sinPhi, -1., 1.)
        delta = math.asin(temp)
        return PointEquatorial(omega, delta)

    def findEquatorialFromHorizontal(self, h: PointHorizontal) -> PointEquatorial:
        v = self.findVectorFromHorizontal(h)
        return self.findEquatorialFromVector(v)

    def findHorizontalFromEquatorial(self, e: PointEquatorial) -> PointHorizontal:
        v = self.findVectorFromEquatorial(e)
        return self.findHorizontalFromVector(v)

    def findHourAngleMax(self, declination: float) -> float:
        phi = self.location.point.latitude()
        temp = -math.tan(phi)*math.tan(declination)
        temp = np.clip(temp, -1., 1.)
        hourAngle = math.acos(temp)
        return hourAngle

    # def findEfromER(self, e: PointEquatorial):
    #     (omegaR, delta) = e
    #     delta = e.declination()
    #     omegaMax = self.findHourAngleMax(delta)
    #     return (omegaR/(6*hour)*omegaMax, delta)
    #
    # def findERfromEquatorial(self, e: PointEquatorial):
    #     delta = e.declination()
    #     omegaMax = self.findHourAngleMax(delta)
    #     return (e.hourAngle()/omegaMax*6*hour, delta)
