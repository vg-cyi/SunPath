import math
import mathutils as mu
import datetime
from datetime import timedelta

degree = math.pi/180
angleHour = 360*degree/24
arcMinute = degree/60
arcSecond = arcMinute/60


class Point:
    x: float
    y: float

    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y


# class Vec2(mu.Vector):
#     @classmethod
#     def fromXY(cls, x: float, y: float):
#         return cls((x, y))
#
# Vec2.fromXY(1., 2.)


class PointHorizontal(Point):
    """
    horizontal coordinates
    """

    def __init__(self, azimuth: float = 0., elevation: float = 0.):
        Point.__init__(self, azimuth, elevation)

    def azimuth(self) -> float:
        return self.x

    def elevation(self) -> float:
        return self.y

    def __repr__(self):
        return 'azimuth = {:.3f}d, elevation = {:.3f}d'.format(self.azimuth()/degree, self.elevation()/degree)


class PointEquatorial(Point):
    """
    equatorial coordinates
    """

    def __init__(self, hourAngle: float = 0., declination: float = 0.):
        Point.__init__(self, hourAngle, declination)

    def hourAngle(self) -> float:
        return self.x

    def declination(self) -> float:
        return self.y

    def __repr__(self):
        return 'hour angle = {:.3f}d, declination = {:.3f}d'.format(self.hourAngle()/degree, self.declination()/degree)


class PointGeographic(Point):
    """
    geographic coordinates
    """

    def __init__(self, latitude: float = 0.*degree, longitude: float = 0.*degree):
        assert abs(latitude) <= 90.*degree
        assert abs(longitude) <= 180.*degree
        Point.__init__(self, longitude, latitude)

    def latitude(self) -> float:
        return self.y

    def longitude(self) -> float:
        return self.x

    def __rer__(self):
        return 'latitude = {:.3f}d, longitude = {:.3f}d'.format(self.latitude()/degree, self.longitude()/degree)


class Location:
    """
    geographic location with name and time zone
    """
    name: str
    point: PointGeographic
    timeDelta: datetime.timedelta

    def __init__(self,
                 name: str = "",
                 point: PointGeographic = PointGeographic(),
                 timeDelta: timedelta = timedelta(hours=0)
                 ):
        self.name = name
        self.point = point
        self.timeDelta = timeDelta
