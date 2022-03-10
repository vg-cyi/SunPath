from .SunCalculator import *

import suncalc

import csv
import urllib.request
from io import StringIO
import re


def irradianceIneichen(s: mu.Vector) -> float:
    """
    Ineichen model for clear sky irradiance (DNI)
    :param s:
    :return:
    """
    if s.z < 0.001:
        return 0.
    else:
        return 1618.*math.exp(-0.606/s.z ** 0.491)


class SunTemporal:
    sc: SunCalculator
    tStep: timedelta
    pointsV: list  # sun positions
    pointsW: list  # irradiance

    def __init__(self, sc: SunCalculator):
        self.sc = sc

    def sampleYear(self, year: int = 2021, tStep: timedelta = timedelta(minutes=60),
                   irradiance: typing.Callable = irradianceIneichen) -> None:
        """
        cache sun positions and irradiance with a fixed time step
        :param year: year
        :param tStep: time step
        :param irradiance: function for clear sky irradiance
        :return:
        """
        sc = self.sc
        self.tStep = tStep

        tz = datetime.timezone(sc.location.timeDelta)
        tA = datetime.datetime(year, 1, 1, 0, 0, 0, tzinfo=tz)
        tB = datetime.datetime(year + 1, 1, 1, 0, 0, 0, tzinfo=tz)

        self.pointsV = []
        latitude = sc.location.point.latitude()/degree
        longitude = sc.location.point.longitude()/degree
        t = tA
        while t < tB:
            sp = suncalc.get_position(t, longitude, latitude)
            h = PointHorizontal(sp['azimuth'] + math.pi, sp['altitude'])
            self.pointsV += [sc.findVectorFromHorizontal(h)]
            t += tStep

        self.pointsW = [irradiance(v) for v in self.pointsV]

    @staticmethod
    def compressTMY(filenameIn: str, filenameOut: str):
        """
        Compress 1-minute data to 1-hour data
        the average irradiance for preceeding hour is found
        """
        fileIn = open(filenameIn, "r")
        reader = csv.reader(file, delimiter=',')
        fileOut = open(filenameOut, 'w')

        header = next(reader)  # Source,Latitude,Longitude,Time Zone
        fileOut.write(header.join(',') + '\n')

        header = next(reader)  # TMY3,0.,0.,0.
        fileOut.write(header.join(',') + '\n')

        header = next(reader)  # Year,Month,Day,Hour,Minute,DNI
        fileOut.write(header.join(',') + '\n')
        nYear = -1
        nMonth = -1
        nDay = -1
        nHour = -1
        nMinute = -1
        nSecond = -1
        nDNI = -1
        for i, h in enumerate(header):
            if re.search('Year', h, re.IGNORECASE):
                nYear = i
            elif re.search('Month', h, re.IGNORECASE):
                nMonth = i
            elif re.search('Day', h, re.IGNORECASE):
                nDay = i
            elif re.search('Hour', h, re.IGNORECASE):
                nHour = i
            elif re.search('Minute', h, re.IGNORECASE):
                nMinute = i
            elif re.search('Second', h, re.IGNORECASE):
                nSecond = i
            elif re.search('(DNI|Beam)', h, re.IGNORECASE):
                nDNI = i
        assert nYear >= 0
        assert nMonth >= 0
        assert nDay >= 0
        assert nHour >= 0
        assert nDNI >= 0

        hourPrev = -1
        irradianceHour = 0.
        irradianceN = 0
        for row in reader:
            year = int(row[nYear])
            month = int(row[nMonth])
            day = int(row[nDay])
            hour = int(row[nHour])
            irradiance = float(row[nDNI])

            irradianceHour += irradiance
            irradianceN += 1
            if hourPrev != hour and irradianceN > 0:
                row[nDNI] = irradianceHour/irradianceN
                if nMinute > 0:
                    row[nMinute] = 0
                if nSecond > 0:
                    row[nSecond] = 0
                fileOut.write(row.join(',') + '\n')
                hourPrev = hour
                irradianceHour = 0.
                irradianceN = 0

    def readTMY(self, filename: str, tMid: float = -0.5) -> None:
        """
        read TMY weather file
        :param filename:
        :param tMid: add this offset (in units of time step) to get the center of interval relative to timeStamp
        :return:
        """
        sc = self.sc

        if filename.startswith("http"):
            response = urllib.request.urlopen(filename)
            file = StringIO(response.read().decode('utf-8'))
        else:
            file = open(filename, "r")

        reader = csv.reader(file, delimiter=',')

        header = next(reader)
        nLatitude = -1
        nLongitude = -1
        nTimeZone = -1
        for i, h in enumerate(header):
            if re.search('Latitude', h, re.IGNORECASE):
                nLatitude = i
            elif re.search('Longitude', h, re.IGNORECASE):
                nLongitude = i
            elif re.search('Time Zone', h, re.IGNORECASE):
                nTimeZone = i
        assert nLatitude >= 0
        assert nLongitude >= 0
        assert nTimeZone >= 0

        row = next(reader)
        latitude = float(row[nLatitude])
        longitude = float(row[nLongitude])
        timeZone = float(row[nTimeZone])

        sc.location.name = "TMY"
        sc.location.point = PointGeographic(latitude*degree, longitude*degree)
        sc.location.timeDelta = datetime.timedelta(hours=timeZone)

        header = next(reader)
        nYear = -1
        nMonth = -1
        nDay = -1
        nHour = -1
        nMinute = -1
        nSecond = -1
        nDNI = -1
        for i, h in enumerate(header):
            if re.search('Year', h, re.IGNORECASE):
                nYear = i
            elif re.search('Month', h, re.IGNORECASE):
                nMonth = i
            elif re.search('Day', h, re.IGNORECASE):
                nDay = i
            elif re.search('Hour', h, re.IGNORECASE):
                nHour = i
            elif re.search('Minute', h, re.IGNORECASE):
                nMinute = i
            elif re.search('Second', h, re.IGNORECASE):
                nSecond = i
            elif re.search('(DNI|Beam)', h, re.IGNORECASE):
                nDNI = i
        assert nYear >= 0
        assert nMonth >= 0
        assert nDay >= 0
        assert nHour >= 0
        assert nDNI >= 0

        listT = []
        listW = []
        tz = datetime.timezone(sc.location.timeDelta)
        for row in reader:
            year = int(row[nYear])
            month = int(row[nMonth])
            day = int(row[nDay])
            hour = int(row[nHour])
            if nMinute > 0:
                minute = int(row[nMinute])
            else:
                minute = 0
            if nSecond > 0:
                second = int(row[nSecond])
            else:
                second = 0
            irradiance = float(row[nDNI])

            t = datetime.datetime(year, month, day, hour, minute, second, tzinfo=tz)
            listT += [t]
            listW += [irradiance]

        self.tStep = listT[1] - listT[0]
        dt = tMid*self.tStep
        # print(dt)

        listV = []
        for t in listT:
            sp = suncalc.get_position(t + dt, longitude, latitude)
            h = PointHorizontal(sp['azimuth'] + math.pi, sp['altitude'])
            listV += [SunCalculator.findVectorFromHorizontal(h)]

        # a correct tMid should minimize wBelow
        wBelow = 0.
        wTotal = 0.
        for v, w in zip(listV, listW):
            wTotal += w
            if v.z < 0. and w > 0.:
                wBelow += w

        ts = self.tStep/datetime.timedelta(hours=1)
        wBelow *= ts
        wTotal *= ts
        print('annual insolation (all elevations): {:.3f} kWh/m2'.format(wTotal/1000))
        print('annual insolation (negative elevations): {:.3f} kWh/m2'.format(wBelow/1000))

        self.pointsV = listV
        self.pointsW = listW

    def checkAccuracy(self, etaF: typing.Callable, etaRefF: typing.Callable, verbose: bool = True) -> list:
        ansRMS = 0.
        ansM = 0.
        ansWM = 0.
        wTotal = 0.
        n = 0
        for v, w in zip(self.pointsV, self.pointsW):
            if v.z < 0.: continue
            diff = etaF(v) - etaRefF(v)
            ansRMS += diff ** 2
            ansM += diff
            ansWM += w*diff
            wTotal += w
            n += 1
        ansRMS = math.sqrt(ansRMS/n)
        ansM /= n
        ansWM /= wTotal
        if verbose:
            print('delta_rms = {:.3f}%'.format(ansRMS*100))
            print('delta_m = {:.3f}%'.format(ansM*100))  # mean
            print('delta_wm = {:.3f}%'.format(ansWM*100))  # weighted mean
        return [ansRMS, ansM, ansWM]

    def integrate(self, etaF: typing.Callable) -> float:
        ans = 0.
        for v, w in zip(self.pointsV, self.pointsW):
            if v.z <= 0.: continue
            eta = etaF(v)
            ans += w*eta

        ans *= self.tStep.seconds/3600.
        return ans

    def average(self, etaF: typing.Callable) -> float:
        a = self.integrate(etaF)
        b = self.integrate(lambda x: 1.)
        return a/b
