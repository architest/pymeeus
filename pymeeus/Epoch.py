# -*- coding: utf-8 -*-


# PyMeeus: Python module implementing astronomical algorithms.
# Copyright (C) 2018  Dagoberto Salazar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import calendar
import datetime
from math import floor
# from base import TOL


"""
.. module:: Epoch
   :synopsis: Class to handle time
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


DAY2SEC = 86400.0
"""Number of seconds per day"""

DAY2MIN = 1440.0
"""Number of minutes per day"""

DAY2HOURS = 24.0
"""Number of hours per day"""

LEAP_TABLE = {1972.5: 1, 1973.0: 2, 1974.0: 3, 1975.0: 4, 1976.0: 5,
              1977.0: 6, 1978.0: 7, 1979.0: 8, 1980.0: 9, 1981.5: 10,
              1982.5: 11, 1983.5: 12, 1985.5: 13, 1988.0: 14, 1990.0: 15,
              1991.0: 16, 1992.5: 17, 1993.5: 18, 1994.5: 19, 1996.0: 20,
              1997.5: 21, 1999.0: 22, 2006.0: 23, 2009.0: 24, 2012.5: 25,
              2015.5: 26, 2017.0: 27}
"""This table represents the point in time FROM WHERE the given number of leap
seconds is valid. Given that leap seconds are (so far) always added at
June 30th or December 31st, a leap second added in 1997/06/30 is represented
here as '1997.5', while a leap second added in 2005/12/31 appears here as
'2006.0'."""


class Epoch(object):
    """
    Class Epoch deals with the tasks related to time handling.

    The constructor takes either a single JDE value, or a series of values
    representing year, month, day, hours, minutes, seconds. This series of
    values is supposed to be in UTC time (civil time). It is also possible to
    provide another Epoch object as input of the constructor.

    When a UTC time is provided, it is converted to International Atomic Time
    (TAI) using an internal table of leap seconds, and from there, it is
    converted to (and stored as) Terrestrial Time (TT). Given that leap seconds
    are added or subtracted in an irregular basis, it is not possible to
    predict them in advance, and the internal leap second table will become
    outdated at some point in time. To counter this, you have two options:

    - Download an updated version of this Pymeeus package.
    - Use the argument 'leap_seconds' in the constructor or 'set()' method to
    provide the correct number of leap seconds (w.r.t. TAI) to be applied.

    For instance, if at some time in the future the TAI-UTC difference is 43
    seconds, you should set 'leap_seconds=43' if you don't have an updated
    version of this class.

    In order to know which is the most updated leap second value stored in this
    class, you may use the 'get_last_leap_second()' method.

    The UTC to TT correction is done by default, but you may disable it by
    setting 'lead_seconds=0'. In that case, it is supposed that the input data
    is already in TT scale.

    :note: Internally, time values are stored as a Julian Ephemeris Day (JDE),
    based on the uniform scale of Dynamical Time, or more specifically,
    Terrestial Time (TT) (itself the redefinition of Terrestrial Dynamical
    Time, TDT).

    :note: The UTC-TT conversion is composed of three parts: a) TT-TAI,
    composed of 32.184 s, b) TAI-UTC(1972), 10 s, and c) UTC(1972)-UTC(Now),
    which is the current amount of leap seconds. When you do 'leap_seconds=43',
    you modify the c) part, while when you do 'leap_seconds=0.0', you disable
    the three corrections.

    :note: Given that this class stores the epoch as JDE, if the JDE value is
    in the order of millions of days then, for a computer with 15-digit
    accuracy, the final time resolution is about 10 milliseconds. That is
    considered enough for most applications of this class.
    """

    def __init__(self, *args, **kwargs):
        """Epoch constructor.

        This constructor takes either a single JDE value, or a series of values
        representing year, month, day, hours, minutes, seconds. This series of
        values is supposed to be in UTC time (civil time).

        It is also possible to provide another Epoch object as input for the
        constructor, or the year, month, etc. arguments can be provided in a
        tuple or list. Moreover, it is also possible provide 'date' or
        'datetime' objects for initialization.

        The 'month' value can be provided as an integer (1 = January, 2 =
        February, etc), or it can be provided as short (Jan, Feb, ...) or long
        (January, February, ...) names. Also, hours, minutes, seconds can be
        provided separately, or as decimals of the day value.

        If 'leap_seconds' argument is set to a value different than zero, then
        that value will be used for the UTC->TAI conversion, and the internal
        leap seconds table will be bypassed. On the other hand, if it is set
        to zero, then the UTC to TT correction is disabled, and it is supposed
        that the input data is already in TT scale.

        :param \*args: Either JDE, Epoch, date, datetime or year, month, day,
        hours, minutes, seconds values, by themselves or inside a tuple or list
        :type \*args: int, float, Epoch, tuple, list, date, datetime
        :param leap_seconds: If different from zero, this is the value to be
        used in the UTC->TAI conversion. If equals to zero, conversion is
        disabled.
        :type leap_seconds: int, float

        :returns: Epoch object.
        :rtype: Epoch
        :raises: ValueError if input values are in the wrong range.
        :raises: TypeError if input values are of wrong type.

        >>> e = Epoch(1987, 6, 19.5, leap_seconds=0.0)
        >>> print(e)
        2446966.0
        """
        # Initialize field
        self._jde = 0.0
        self.set(*args, **kwargs)   # Use 'set()' method to handle the setup

    def set(self, *args, **kwargs):
        """Method used to set the value of this object.

        This method takes either a single JDE value, or a series of values
        representing year, month, day, hours, minutes, seconds. This series of
        values is supposed to be in UTC time (civil time).

        It is also possible to provide another Epoch object as input for the
        'set()' method, or the year, month, etc arguments can be provided in a
        tuple or list. Moreover, it is also possible provide 'date' or
        'datetime' objects for initialization.

        The 'month' value can be provided as an integer (1 = January, 2 =
        February, etc), or it can be provided as short (Jan, Feb, ...) or long
        (January, February, ...) names. Also, hours, minutes, seconds can be
        provided separately, or as decimals of the day value.

        If 'leap_seconds' argument is set to a value different than zero, then
        that value will be used for the UTC->TAI conversion, and the internal
        leap seconds table will be bypassed. On the other hand, if it is set
        to zero, then the UTC to TT correction is disabled, and it is supposed
        that the input data is already in TT scale.

        :param \*args: Either JDE, Epoch, date, datetime or year, month, day,
        hours, minutes, seconds values, by themselves or inside a tuple or list
        :type \*args: int, float, Epoch, tuple, list, date, datetime
        :param leap_seconds: If different from zero, this is the value to be
        used in the UTC->TAI conversion. If equals to zero, conversion is
        disabled. If not given, UTC to TT conversion is carried out (default).
        :type leap_seconds: int, float

        :returns: None.
        :rtype: None
        :raises: ValueError if input values are in the wrong range.
        :raises: TypeError if input values are of wrong type.

        >>> e = Epoch()
        >>> e.set(1987, 6, 19.5, leap_seconds=0.0)
        >>> print(e)
        2446966.0
        >>> e.set(1977, 'Apr', 26.4, leap_seconds=0.0)
        >>> print(e)
        2443259.9
        >>> e.set(1957, 'October', 4.81, leap_seconds=0.0)
        >>> print(e)
        2436116.31
        >>> e.set(333, 'Jan', 27, 12, leap_seconds=0.0)
        >>> print(e)
        1842713.0
        >>> e.set(1900, 'Jan', 1, leap_seconds=0.0)
        >>> print(e)
        2415020.5
        >>> e.set(-1001, 'august', 17.9, leap_seconds=0.0)
        >>> print(e)
        1355671.4
        >>> e.set(-4712, 1, 1.5, leap_seconds=0.0)
        >>> print(e)
        0.0
        >>> e.set((1600, 12, 31), leap_seconds=0.0)
        >>> print(e)
        2305812.5
        >>> e.set([1988, 'JUN', 19, 12], leap_seconds=0.0)
        >>> print(e)
        2447332.0
        >>> d = datetime.date(2000, 1, 1)
        >>> e.set(d, leap_seconds=0.0)
        >>> print(e)
        2451544.5
        >>> e.set(837, 'Apr', 10, 7, 12, leap_seconds=0.0)
        >>> print(e)
        2026871.8
        >>> d = datetime.datetime(837, 4, 10, 7, 12, 0, 0)
        >>> e.set(d, leap_seconds=0.0)
        >>> print(e)
        2026871.8
        """
        # Clean up the internal parameters
        self._jde = 0.0
        # If no arguments are given, return. Internal values are 0.0
        if len(args) == 0:
            return
        # If we have only one argument, it can be a JDE or another Epoch object
        elif len(args) == 1:
            if isinstance(args[0], Epoch):
                self._jde = args[0]._jde
                return
            elif isinstance(args[0], (int, float)):
                self._jde = args[0]
                return
            elif isinstance(args[0], (tuple, list)):
                year, month, day, hours, minutes, sec = \
                        self._check_values(*args[0])
            elif isinstance(args[0], datetime.datetime):
                d = args[0]
                year, month, day, hours, minutes, sec = \
                    self._check_values(d.year, d.month, d.day, d.hour,
                                       d.minute, d.second + d.microsecond/1e6)
            elif isinstance(args[0], datetime.date):
                d = args[0]
                year, month, day, hours, minutes, sec = \
                    self._check_values(d.year, d.month, d.day)
            else:
                raise TypeError("Invalid input value")
        elif len(args) == 2:
            # Insuficient data to set the Epoch
            raise ValueError("Invalid number of input values")
        elif len(args) >= 3:        # Year, month, day
            year, month, day, hours, minutes, sec = self._check_values(*args)
        day += hours/DAY2HOURS + minutes/DAY2MIN + sec/DAY2SEC
        # Handle the 'leap_seconds' argument, if pressent
        if 'leap_seconds' in kwargs:
            # Compute JDE
            self._jde = self._compute_jde(year, month, day, utc2tt=False,
                                          leap_seconds=kwargs['leap_seconds'])
        else:
            self._jde = self._compute_jde(year, month, day)

    def _compute_jde(self, y, m, d, utc2tt=True, leap_seconds=0.0):
        """Method to compute the Julian Ephemeris Day (JDE).

        :param y: Year
        :type y: int
        :param m: Month
        :type m: int
        :param d: Day
        :type d: float
        :param utc2tt: Whether correction UTC to TT is done automatically.
        :type utc2tt: bool
        :param leap_seconds: Number of leap seconds to apply
        :type leap_seconds: float

        :returns: Julian Ephemeris Day (JDE)
        :rtype: float
        """

        # The best approach here is first convert to JDE, and then adjust secs
        if m <= 2:
            y -= 1
            m += 12
        a = floor(y/100.0)
        b = 0.0
        if not Epoch.is_julian(y, m, floor(d)):
            b = 2.0 - a + floor(a/4.0)
        jde = floor(365.25*(y + 4716.0)) + floor(30.6001*(m + 1.0)) + \
            d + b - 1524.5
        # If enabled, let's convert from UTC to TT, adding the needed seconds
        deltasec = 0.0
        # In this case, UTC to TT correction is applied automatically
        if utc2tt:
            deltasec = 32.184       # Difference between TT and TAI
            if y > 1972 or (y == 1972 and m >= 7):
                deltasec += 10.0    # Difference between UTC and TAI in 1972
                deltasec += Epoch.leap_seconds(y, m)
        else:                           # Correction is NOT automatic
            if leap_seconds != 0.0:     # We apply provided leap seconds
                deltasec = 32.184       # Difference between TT and TAI
                if y > 1972 or (y == 1972 and m >= 7):
                    deltasec += 10.0    # Difference between UTC-TAI in 1972
                    deltasec += leap_seconds
        return jde + deltasec/DAY2SEC

    def _check_values(self, *args):
        """This method takes the input arguments to 'set()' method (year,
        month, day, etc) and carries out some sanity checks on them.

        It returns a tuple containing those values separately, assigning zeros
        to those arguments which were not provided.

        :param \*args: Year, month, day, hours, minutes, seconds values.
        :type \*args: int, float

        :returns: Tuple with year, month, day, hours, minutes, seconds values.
        :rtype: tuple
        :raises: ValueError if input values are in the wrong range, or too few
        arguments given as input.
        """

        # This list holds the maximum amount of days a given month can have
        maxdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # Initialize some variables
        year = -9999
        month = -9999
        day = -9999
        hours = 0.0
        minutes = 0.0
        sec = 0.0

        # Carry out some basic checks
        if len(args) < 3:
            raise ValueError("Invalid number of input values")
        elif len(args) >= 3:        # Year, month, day
            year = args[0]
            month = args[1]
            day = args[2]
        if len(args) >= 4:          # Year, month, day, hour
            hours = args[3]
        if len(args) >= 5:          # Year, month, day, hour, minutes
            minutes = args[4]
        if len(args) >= 6:          # Year, month, day, hour, minutes, seconds
            sec = args[5]
        if year < -4712:        # No negative JDE will be allowed
            raise ValueError("Invalid value for the input year")
        if day < 1 or day > 31:
            raise ValueError("Invalid value for the input day")
        if hours < 0 or hours > 23:
            raise ValueError("Invalid value for the input hours")
        if minutes < 0 or minutes > 59:
            raise ValueError("Invalid value for the input minutes")
        if sec < 0 or sec > 59:
            raise ValueError("Invalid value for the input seconds")

        # Test the days according to the month
        month = Epoch.get_month(month)
        limit_day = maxdays[month - 1]
        # We need extra tests if month is '2' (February)
        if month == 2:
            if Epoch.is_leap(year):
                limit_day = 29
        if day > limit_day:
            raise ValueError("Invalid value for the input day")

        # We are ready to return the parameters
        return year, month, day, hours, minutes, sec

    @staticmethod
    def is_julian(year, month, day):
        """This method returns True if date given is in the Julian calendar.

        :param year: Year
        :type y: int
        :param month: Month
        :type m: int
        :param day: Day
        :type day: int

        :returns: Whether the provided date belongs to Julian calendar or not.
        :rtype: bool

        >>> Epoch.is_julian(1997, 5, 27.1)
        False
        >>> Epoch.is_julian(1397, 7, 7.0)
        True
        """
        if year < 1582:
            return True
        elif year == 1582 and month < 10:
            return True
        elif year == 1582 and month == 10 and day < 5.0:
            return True
        else:
            return False

    @staticmethod
    def get_month(month):
        """Method to get the month as a integer in the [1, 12] range.

        :param month: Month, in numeric, short name or long name format
        :type DD: int, float, str

        :returns: Month as integer in the [1, 12] range.
        :rtype: int
        :raises: ValueError if input month value is invalid.

        >>> Epoch.get_month(4.0)
        4
        >>> Epoch.get_month('Oct')
        10
        >>> Epoch.get_month('FEB')
        2
        >>> Epoch.get_month('August')
        8
        >>> Epoch.get_month('NOVEMBER')
        11
        """

        months_mmm = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

        months_full = ["January", "February", "March", "April", "May", "June",
                       "July", "August", "September", "October", "November",
                       "December"]

        if isinstance(month, (int, float)):
            month = floor(month)            # Truncate if it has decimals
            if month >= 1 and month <= 12:
                return int(month)
            else:
                raise ValueError("Invalid value for the input month")
        elif isinstance(month, str):
            month = month.strip().capitalize()
            if len(month) == 3:
                if month in months_mmm:
                    return (months_mmm.index(month) + 1)
                else:
                    raise ValueError("Invalid value for the input month")
            else:
                if month in months_full:
                    return (months_full.index(month) + 1)
                else:
                    raise ValueError("Invalid value for the input month")

    @staticmethod
    def is_leap(year):
        """Method to check if a given year is a leap year.

        :param year: Year to be checked.
        :type year: int, float

        :returns: Whether or not year is a leap year.
        :rtype: bool
        :raises: ValueError if input year value is invalid.

        >>> Epoch.is_leap(2003)
        False
        >>> Epoch.is_leap(2012)
        True
        >>> Epoch.is_leap(-1000)
        True
        >>> Epoch.is_leap(1000)
        True
        """
        if isinstance(year, (int, float)):
            # Mind the difference between Julian and Gregorian calendars
            if year >= 1582:
                year = floor(year)
                return calendar.isleap(year)
            else:
                return (abs(year) % 4) == 0
        else:
            raise ValueError("Invalid value for the input year")

    @staticmethod
    def getDOY(YYYY, MM, DD):
        """This method returns the Day Of Year (DOY) for the given date.

        :param YYYY: Year, in four digits format
        :type YYYY: int, float
        :param MM: Month, in numeric format (1 = January, 2 = February, etc)
        :type MM: int, float
        :param DD: Day, in numeric format
        :type DD: int, float

        :returns: Day Of Year (DOY).
        :rtype: float
        :raises: ValueError if input values correspond to a wrong date.

        >>> Epoch.getDOY(1999, 1, 29)
        29.0
        >>> Epoch.getDOY(1978, 11, 14)
        318.0
        >>> Epoch.getDOY(2017, 12, 31.7)
        365.7
        >>> Epoch.getDOY(2012, 3, 3.1)
        63.1
        >>> Epoch.getDOY(-400, 2, 29.9)
        60.9
        """
        # Let's carry out first some basic checks
        if DD < 1 or DD >= 32 or MM < 1 or MM > 12:
            raise ValueError("Invalid input data")
        day = int(DD)
        frac = DD % 1
        if YYYY >= 1:                           # datetime's minimum year is 1
            try:
                d = datetime.date(YYYY, MM, day)
            except ValueError:
                raise ValueError("Invalid input date")
            doy = d.timetuple().tm_yday
        else:
            k = 2 if Epoch.is_leap(YYYY) else 1
            doy = floor((275.0*MM)/9.0) - k*floor((MM + 9.0)/12.0) + day - 30.0
        return float(doy + frac)

    @staticmethod
    def doy2date(year, doy):
        """This method takes a year and a Day Of Year values, and returns the
        corresponding date.

        :param year: Year, in four digits format
        :type year: int, float
        :param doy: Day of Year number
        :type doy: int, float

        :returns: Year, month, day.
        :rtype: tuple
        :raises: ValueError if either input year or doy values are invalid.

        >>> t = Epoch.doy2date(1999, 29)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        1999/1/29.0
        >>> t = Epoch.doy2date(2017, 365.7)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        2017/12/31.7
        >>> t = Epoch.doy2date(2012, 63.1)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        2012/3/3.1
        >>> t = Epoch.doy2date(-1004, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -1004/2/29.0
        >>> t = Epoch.doy2date(0, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        0/2/29.0
        >>> t = Epoch.doy2date(1, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        1/3/1.0
        >>> t = Epoch.doy2date(-1, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -1/3/1.0
        >>> t = Epoch.doy2date(-2, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -2/3/1.0
        >>> t = Epoch.doy2date(-3, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -3/3/1.0
        >>> t = Epoch.doy2date(-4, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -4/2/29.0
        >>> t = Epoch.doy2date(-5, 60)
        >>> print("{}/{}/{}".format(t[0], t[1], round(t[2], 1)))
        -5/3/1.0
        """
        if isinstance(year, (int, float)) and isinstance(doy, (int, float)):
            frac = float(doy % 1)
            doy = int(doy)
            if year >= 1:                   # datetime's minimum year is 1
                ref = datetime.date(year, 1, 1)
                mydate = datetime.date.fromordinal(ref.toordinal() + doy - 1)
                return year, mydate.month, mydate.day + frac
            else:
                # The algorithm provided by Meeus doesn't work for years below
                # +1. This little hack solves that problem (the 'if' result is
                # inverted here).
                k = 1 if Epoch.is_leap(year) else 2
                if doy < 32:
                    m = 1
                else:
                    m = floor((9.0*(k + doy))/275.0 + 0.98)
                d = doy - floor((275.0*m)/9.0) + k*floor((m + 9.0)/12.0) + 30
                return year, int(m), d + frac
        else:
            raise ValueError("Invalid input values")

    @staticmethod
    def leap_seconds(year, month):
        """Returns the leap seconds accumulated for the given year and month.

        :param year: Year
        :type year: int
        :param month: Month, in numeric format ([1:12] range)
        :type month: int

        :returns: Leap seconds accumulated for given year and month.
        :rtype: int

        >>> Epoch.leap_seconds(1972, 4)
        0
        >>> Epoch.leap_seconds(1972, 6)
        0
        >>> Epoch.leap_seconds(1972, 7)
        1
        >>> Epoch.leap_seconds(1983, 6)
        11
        >>> Epoch.leap_seconds(1983, 7)
        12
        >>> Epoch.leap_seconds(1985, 8)
        13
        >>> Epoch.leap_seconds(2016, 11)
        26
        >>> Epoch.leap_seconds(2017, 1)
        27
        >>> Epoch.leap_seconds(2018, 7)
        27
        """

        list_years = sorted(LEAP_TABLE.keys())
        # First test the extremes of the table
        if (year + month/12.0) <= list_years[0]:
            return 0
        if (year + month/12.0) >= list_years[-1]:
            return LEAP_TABLE[list_years[-1]]
        lyear = (year + 0.25) if month <= 6 else (year + 0.75)
        idx = 0
        while lyear > list_years[idx]:
            idx += 1
        return LEAP_TABLE[list_years[idx - 1]]

    @staticmethod
    def get_last_leap_second():
        """Method to get the date and value of the last leap second added to
        the table

        :returns: Tuple with year, month, day, leap second value.
        :rtype: tuple
        """
        list_years = sorted(LEAP_TABLE.keys())
        lyear = list_years[-1]
        lseconds = LEAP_TABLE[lyear]
        year = floor(lyear)
        # So far, leap seconds are added either on June 30th or December 31th
        if lyear % 1 == 0.0:
            year -= 1
            month = 12
            day = 31.0
        else:
            month = 6
            day = 30.0
        return year, month, day, lseconds

    def __str__(self):
        """Method used when trying to print the object.

        :returns: Internal JDE value as a string.
        :rtype: string

        >>> e = Epoch(1987, 6, 19.5, leap_seconds=0.0)
        >>> print(e)
        2446966.0
        """
        return str(self._jde)

    def get_date(self, **kwargs):
        """This method converts the internal JDE value back to a date.

        Use 'leap_seconds=0.0' to disable the automatic TT to UTC conversion
        mechanism, or provide a non zero value to 'leap_seconds' to apply a
        specific leap seconds value.

        :param leap_seconds: Optional value for leap seconds.
        :type leap_seconds: int, float

        :returns: Year, month, day in a tuple
        :rtype: tuple

        >>> e = Epoch(2436116.31)
        >>> y, m, d = e.get_date(leap_seconds=0.0)
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        1957/10/4.81
        >>> e = Epoch(1988, 1, 27)
        >>> y, m, d = e.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        1988/1/27.0
        >>> e = Epoch(1842713.0)
        >>> y, m, d = e.get_date(leap_seconds=0.0)
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        333/1/27.5
        >>> e = Epoch(1507900.13)
        >>> y, m, d = e.get_date(leap_seconds=0.0)
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        -584/5/28.63
        """
        jd = self._jde + 0.5
        z = floor(jd)
        f = jd % 1
        if z < 2299161:
            a = z
        else:
            alpha = floor((z - 1867216.25)/36524.25)
            a = z + 1 + alpha - floor(alpha/4.0)
        b = a + 1524
        c = floor((b - 122.1)/365.25)
        d = floor(365.25*c)
        e = floor((b - d)/30.6001)
        day = b - d - floor(30.6001*e) + f
        if e < 14:
            month = e - 1
        elif e == 14 or e == 15:
            month = e - 13
        if month > 2:
            year = c - 4716
        elif month == 1 or month == 2:
            year = c - 4715
        year = int(year)
        month = int(month)

        tt2utc = True
        if 'leap_seconds' in kwargs:
            tt2utc = False
            leap_seconds = kwargs['leap_seconds']
        # If enabled, let's convert from TT to UTC, subtracting needed seconds
        deltasec = 0.0
        # In this case, TT to UTC correction is applied automatically
        if tt2utc:
            deltasec = 32.184       # Difference between TT and TAI
            if year > 1972 or (year == 1972 and month >= 7):
                deltasec += 10.0    # Difference between UTC and TAI in 1972
                deltasec += Epoch.leap_seconds(year, month)
        else:                           # Correction is NOT automatic
            if leap_seconds != 0.0:     # We apply provided leap seconds
                deltasec = 32.184       # Difference between TT and TAI
                if year > 1972 or (year == 1972 and month >= 7):
                    deltasec += 10.0    # Difference between UTC-TAI in 1972
                    deltasec += leap_seconds

        if deltasec != 0.0:
            doy = Epoch.getDOY(year, month, day)
            doy -= deltasec/DAY2SEC
            # Check that we didn't change year
            if doy < 1.0:
                year -= 1
                doy = 366.0 + doy if Epoch.is_leap(year) else 365.0 + doy
            year, month, day = Epoch.doy2date(year, doy)
        return year, month, day

    def dow(self, as_string=False):
        """Method to return the day of week corresponding to this Epoch.

        By default, this method returns an integer value: 0 for Sunday, 1 for
        Monday, etc. However, when 'as_string=True' is passed, the names of the
        days are returned.

        :param as_string: Whether result will be given as a integer or as a
        string. False by default.
        :type as_string: bool

        :returns: Day of the week, as a integer or as a string.
        :rtype: int, str

        >>> e = Epoch(1954, 'June', 30)
        >>> e.dow()
        3
        >>> e = Epoch(2018, 'Feb', 14.9)
        >>> e.dow(as_string=True)
        'Wednesday'
        >>> e = Epoch(2018, 'Feb', 15)
        >>> e.dow(as_string=True)
        'Thursday'
        >>> e = Epoch(2018, 'Feb', 15.99)
        >>> e.dow(as_string=True)
        'Thursday'
        >>> e.set(2018, 'Jul', 15.4)
        >>> e.dow(as_string=True)
        'Sunday'
        >>> e.set(2018, 'Jul', 15.9)
        >>> e.dow(as_string=True)
        'Sunday'
        """
        jd = floor(self._jde - 0.5) + 2.0
        doy = int(jd % 7)
        if not as_string:
            return doy
        else:
            day_names = ['Sunday', 'Monday', 'Tuesday', 'Wednesday',
                         'Thursday', 'Friday', 'Saturday']
            return day_names[doy]

    def mjd(self):
        """This method returns the Modified Julian Day (MJD).

        :returns: Modified Julian Day (MJD).
        :rtype: float

        >>> e = Epoch(1858, 'NOVEMBER', 17, leap_seconds=0.0)
        >>> e.mjd()
        0.0
        """
        return self._jde - 2400000.5

    def jde(self):
        """Method to return the internal value of the Julian Ephemeris Day.

        :returns: The internal value of the Julian Ephemeris Day.
        :rtype: float

        >>> a = Epoch(-1000, 2, 29.0, leap_seconds=0.0)
        >>> print(a.jde())
        1355866.5
        """
        return self()

    def __call__(self):
        """Method used when Epoch is called only with parenthesis.

        :returns: The internal value of the Julian Ephemeris Day.
        :rtype: float

        >>> a = Epoch(-122, 1, 1.0, leap_seconds=0.0)
        >>> print(a())
        1676497.5
        """
        return self._jde

    def __add__(self, b):
        """This method defines the addition between an Epoch and some days.

        :param b: Value to be added, in days.
        :type b: int, float

        :returns: A new Epoch object.
        :rtype: Epoch
        :raises: TypeError if operand is of wrong type.

        >>> a = Epoch(1991, 7, 11)
        >>> b = a + 10000
        >>> y, m, d = b.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        2018/11/26.0
        """
        if isinstance(b, (int, float)):
            return Epoch(self._jde + float(b))
        else:
            raise TypeError("Wrong operand type")

    def __sub__(self, b):
        """This method defines the subtraction between Epochs or between an
        Epoch and a given number of days.

        :param b: Value to be subtracted, either an Epoch or days.
        :type b: Epoch, int, float

        :returns: A new Epoch object if parameter 'b' is in days, or the
        difference between provided Epochs, in days.
        :rtype: Epoch, float
        :raises: TypeError if operand is of wrong type.

        >>> a = Epoch(1986, 2, 9.0)
        >>> print(round(a(), 2))
        2446470.5
        >>> b = Epoch(1910, 4, 20.0)
        >>> print(round(b(), 2))
        2418781.5
        >>> c = a - b
        >>> print(round(c, 2))
        27689.0
        >>> a = Epoch(2003, 12, 31.0)
        >>> b = a - 365.5
        >>> y, m, d = b.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        2002/12/30.5
        """
        if isinstance(b, (int, float)):
            return Epoch(self._jde - b)
        elif isinstance(b, Epoch):
            return float(self._jde - b._jde)
        else:
            raise TypeError("Invalid operand type")

    def __iadd__(self, b):
        """This method defines the accumulative addition to this Epoch.

        :param b: Value to be added, in days.
        :type b: int, float

        :returns: This Epoch.
        :rtype: Epoch
        :raises: TypeError if operand is of wrong type.

        >>> a = Epoch(2003, 12, 31.0)
        >>> a += 32.5
        >>> y, m, d = a.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        2004/2/1.5
        """
        if isinstance(b, (int, float)):
            self = self + b
            return self
        else:
            raise TypeError("Wrong operand type")

    def __isub__(self, b):
        """This method defines the accumulative subtraction to this Epoch.

        :param b: Value to be subtracted, in days.
        :type b: int, float

        :returns: This Epoch.
        :rtype: Epoch
        :raises: TypeError if operand is of wrong type.

        >>> a = Epoch(2001, 12, 31.0)
        >>> a -= 2*365
        >>> y, m, d = a.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        2000/1/1.0
        """
        if isinstance(b, (int, float)):
            self = self - b
            return self
        else:
            raise TypeError("Wrong operand type")

    def __radd__(self, b):
        """This method defines the addition to a Epoch by the right

        :param b: Value to be added, in days.
        :type b: int, float

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if operand is of wrong type.

        >>> a = Epoch(2004, 2, 27.8)
        >>> b = 2.2 + a
        >>> y, m, d = b.get_date()
        >>> print("{}/{}/{}".format(y, m, round(d, 2)))
        2004/3/1.0
        """
        if isinstance(b, (int, float)):
            return self.__add__(b)              # It is the same as by the left
        else:
            raise TypeError("Wrong operand type")

    def __int__(self):
        """This method returns the internal JDE value as an int.

        :returns: Internal JDE value as an int.
        :rtype: int

        >>> a = Epoch(2434923.85)
        >>> int(a)
        2434923
        """
        return int(self._jde)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's do some work with the Epoch class
    print('\n' + 35*'*')
    print("*** Use of Epoch class")
    print(35*'*' + '\n')

    e = Epoch(1987, 6, 19.5)
    print_me("JDE for 1987/6/19.5", e)

    # Redefine the Epoch object
    e.set(333, 'Jan', 27, 12, leap_seconds=0.0)
    print_me("JDE for 333/1/27.5", e)

    # We can create an Epoch from a 'date' or 'datetime' object
    d = datetime.datetime(837, 4, 10, 7, 12, 0, 0)
    f = Epoch(d, leap_seconds=0.0)
    print_me("JDE for 837/4/10.3", f)


if __name__ == '__main__':

    main()
