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


from math import sqrt, sin, cos, tan, atan, atan2, asin, acos
from base import TOL
from Angle import Angle
from Epoch import Epoch, JDE2000


"""
.. module:: Coordinates
   :synopsis: Module including different functions to handle coordinates
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


NUTATION_ARG_TABLE = [
    [0, 0, 0, 0, 1], [-2, 0, 0, 2, 2], [0, 0, 0, 2, 2], [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [-2, 1, 0, 2, 2], [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2], [-2, -1, 0, 2, 2], [-2, 0, 1, 0, 0], [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2], [2, 0, 0, 0, 0], [0, 0, 1, 0, 1], [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1], [0, 0, 1, 2, 1], [-2, 0, 2, 0, 0], [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2], [0, 0, 2, 2, 2], [0, 0, 2, 0, 0], [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0], [-2, 0, 0, 2, 0], [0, 0, -1, 2, 1], [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1], [-2, 2, 0, 2, 2], [0, 1, 0, 0, 1], [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1], [0, 0, 2, -2, 0], [2, 0, -1, 2, 1], [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2], [-2, 1, 1, 0, 0], [0, -1, 0, 2, 2], [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0], [-2, 0, 2, 2, 2], [-2, 0, 1, 2, 1], [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1], [0, -1, 1, 0, 0], [-2, -1, 0, 2, 1], [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1], [-2, 0, 2, 0, 1], [-2, 1, 0, 2, 1], [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0], [-2, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2], [-1, -1, 1, 0, 0], [0, 1, 1, 0, 0], [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2], [0, 0, 3, 2, 2], [2, -1, 0, 2, 2]]
"""This table contains the periodic terms for the argument of the nutation. In
Meeus' book this is Table 22.A and can be found in pages 145-146."""

NUTATION_SINE_COEF_TABLE = [
    [-171996.0, -174.2], [-13187.0, -1.6], [-2274.0, -0.2], [2062.0, 0.2],
    [1426.0, -3.4], [712.0, 0.1], [-517.0, 1.2], [-386.0, -0.4], [-301.0, 0.0],
    [217.0, -0.5], [-158.0, 0.0], [129.0, 0.1], [123.0, 0.0], [63.0, 0.0],
    [63.0, 0.1], [-59.0, 0.0], [-58.0, -0.1], [-51.0, 0.0], [48.0, 0.0],
    [46.0, 0.0], [-38.0, 0.0], [-31.0, 0.0], [29.0, 0.0], [29.0, 0.0],
    [26.0, 0.0], [-22.0, 0.0], [21.0, 0.0], [17.0, -0.1], [16.0, 0.0],
    [-16.0, 0.1], [-15.0, 0.0], [-13.0, 0.0], [-12.0, 0.0], [11.0, 0.0],
    [-10.0, 0.0], [-8.0, 0.0], [7.0, 0.0], [-7.0, 0.0], [-7.0, 0.0],
    [-7.0, 0.0], [6.0, 0.0], [6.0, 0.0], [6.0, 0.0], [-6.0, 0.0], [-6.0, 0.0],
    [5.0, 0.0], [-5.0, 0.0], [-5.0, 0.0], [-5.0, 0.0], [4.0, 0.0], [4.0, 0.0],
    [4.0, 0.0], [-4.0, 0.0], [-4.0, 0.0], [-4.0, 0.0], [3.0, 0.0], [-3.0, 0.0],
    [-3.0, 0.0], [-3.0, 0.0], [-3.0, 0.0], [-3.0, 0.0], [-3.0, 0.0],
    [-3.0, 0.0]]
"""This table contains the periodic terms for the coefficients of the sine of
the argument of the nutation, and they are used to compute Delta psi. Units are
in 0.0001''. In Meeus' book this is Table 22.A and can be found in pages
145-146."""

NUTATION_COSINE_COEF_TABLE = [
    [92025.0, 8.9], [5736.0, -3.1], [977.0, -0.5], [-895.0, 0.5], [54.0, -0.1],
    [-7.0, 0.0], [224.0, -0.6], [200.0, 0.0], [129.0, -0.1], [-95.0, 0.3],
    [0.0, 0.0], [-70.0, 0.0], [-53.0, 0.0], [0.0, 0.0], [-33.0, 0.0],
    [26.0, 0.0], [32.0, 0.0], [27.0, 0.0], [0.0, 0.0], [-24.0, 0.0],
    [16.0, 0.0], [13.0, 0.0], [0.0, 0.0], [-12.0, 0.0], [0.0, 0.0], [0.0, 0.0],
    [-10.0, 0.0], [0.0, 0.0], [-8.0, 0.0], [7.0, 0.0], [9.0, 0.0], [7.0, 0.0],
    [6.0, 0.0], [0.0, 0.0], [5.0, 0.0], [3.0, 0.0], [-3.0, 0.0], [0.0, 0.0],
    [3.0, 0.0], [3.0, 0.0], [0.0, 0.0], [-3.0, 0.0], [-3.0, 0.0], [3.0, 0.0],
    [3.0, 0.0], [0.0, 0.0], [3.0, 0.0], [3.0, 0.0], [3.0, 0.0]]
"""This table contains the periodic terms for the coefficients of the cosine of
the argument of the nutation, and they are used to compute Delta epsilon. Units
are in 0.0001''. In Meeus' book this is Table 22.A and can be found in pages
145-146."""


def mean_obliquity(*args, **kwargs):
    """This function computes the mean obliquity (epsilon0) at the provided
    date.

    This function internally uses an :class:`Epoch` object, and the **utc**
    argument then controls the way the UTC->TT conversion is handled for that
    object. If **leap_seconds** argument is set to a value different than zero,
    then that value will be used for the UTC->TAI conversion, and the internal
    leap seconds table will be bypassed.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param utc: Whether the provided epoch is a civil time (UTC) or TT
    :type utc: bool
    :param leap_seconds: This is the value to be used in the UTC->TAI
        conversion, instead of taking it from internal leap seconds table.
    :type leap_seconds: int, float

    :returns: The mean obliquity of the ecliptic, as an :class:`Angle`
    :rtype: :class:`Angle`
    :raises: ValueError if input values are in the wrong range.
    :raises: TypeError if input values are of wrong type.

    >>> e0 = mean_obliquity(1987, 4, 10)
    >>> a = e0.dms_tuple()
    >>> a[0]
    23
    >>> a[1]
    26
    >>> round(a[2], 3)
    27.407
    """

    # Get the Epoch object corresponding to input parameters
    t = Epoch.check_input_date(*args, **kwargs)
    # Let's redefine u in units of 100 Julian centuries from Epoch J2000.0
    u = (t.jde() - 2451545.0)/3652500.0
    epsilon0 = Angle(23, 26, 21.448)
    delta = u*(-4680.93 + u*(-1.55 + u*(1999.25 + u*(-51.38 + u*(-249.67 +
               u*(-39.05 + u*(7.12 + u*(27.87 + u*(5.79 + u*2.45)))))))))
    delta = Angle(0, 0, delta)
    epsilon0 += delta
    return epsilon0


def true_obliquity(*args, **kwargs):
    """This function computes the true obliquity (epsilon) at the provided
    date. The true obliquity is the mean obliquity (epsilon0) plus the
    correction provided by the nutation in obliquity (Delta epsilon).

    This function internally uses an :class:`Epoch` object, and the **utc**
    argument then controls the way the UTC->TT conversion is handled for that
    object. If **leap_seconds** argument is set to a value different than zero,
    then that value will be used for the UTC->TAI conversion, and the internal
    leap seconds table will be bypassed.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param utc: Whether the provided epoch is a civil time (UTC) or TT
    :type utc: bool
    :param leap_seconds: This is the value to be used in the UTC->TAI
        conversion, instead of taking it from internal leap seconds table.
    :type leap_seconds: int, float

    :returns: The true obliquity of the ecliptic, as an Angle
    :rtype: :class:`Angle`
    :raises: ValueError if input values are in the wrong range.
    :raises: TypeError if input values are of wrong type.

    >>> epsilon = true_obliquity(1987, 4, 10)
    >>> a = epsilon.dms_tuple()
    >>> a[0]
    23
    >>> a[1]
    26
    >>> round(a[2], 3)
    36.849
    """

    epsilon0 = mean_obliquity(*args, **kwargs)
    delta_epsilon = nutation_obliquity(*args, **kwargs)
    return (epsilon0 + delta_epsilon)


def nutation_longitude(*args, **kwargs):
    """This function computes the nutation in longitude (Delta psi) at the
    provided date.

    This function internally uses an :class:`Epoch` object, and the **utc**
    argument then controls the way the UTC->TT conversion is handled for that
    object. If **leap_seconds** argument is set to a value different than zero,
    then that value will be used for the UTC->TAI conversion, and the internal
    leap seconds table will be bypassed.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param utc: Whether the provided epoch is a civil time (UTC) or TT
    :type utc: bool
    :param leap_seconds: This is the value to be used in the UTC->TAI
        conversion, instead of taking it from internal leap seconds table.
    :type leap_seconds: int, float

    :returns: The nutation in longitude (Delta psi), as an Angle
    :rtype: :class:`Angle`
    :raises: ValueError if input values are in the wrong range.
    :raises: TypeError if input values are of wrong type.

    >>> dpsi = nutation_longitude(1987, 4, 10)
    >>> a = dpsi.dms_tuple()
    >>> a[0]
    0
    >>> a[1]
    0
    >>> round(a[2], 3)
    3.788
    >>> a[3]
    -1.0
    """

    # Get the Epoch object corresponding to input parameters
    t = Epoch.check_input_date(*args, **kwargs)
    # Let's redefine t in units of Julian centuries from Epoch J2000.0
    t = (t.jde() - 2451545.0)/36525.0
    # Let's compute the mean elongation of the Moon from the Sun
    d = 297.85036 + t*(445267.111480 + t*(-0.0019142 + t/189474.0))
    d = Angle(d)            # Convert into an Angle: It is easier to handle
    # Compute the mean anomaly of the Sun (from Earth)
    m = 357.52772 + t*(35999.050340 + t*(-0.0001603 - t/300000.0))
    m = Angle(m)
    # Compute the mean anomaly of the Moon
    mprime = 134.96298 + t*(477198.867398 + t*(0.0086972 + t/56250.0))
    mprime = Angle(mprime)
    # Now, let's compute the Moon's argument of latitude
    f = 93.27191 + t*(483202.017538 + t*(-0.0036825 + t/327270.0))
    f = Angle(f)
    # And finally, the longitude of the ascending node of the Moon's mean
    # orbit on the ecliptic, measured from the mean equinox of date
    omega = 125.04452 + t*(-1934.136261 + t*(0.0020708 + t/450000.0))
    omega = Angle(omega)
    # Let's store this results in a list, in preparation for using tables
    arguments = [d, m, mprime, f, omega]
    # Now is time of using the nutation tables
    deltapsi = 0.0
    for i in range(len(NUTATION_SINE_COEF_TABLE)):
        argument = Angle()
        coeff = 0.0
        for j in range(5):
            if NUTATION_ARG_TABLE[i][j]:    # Avoid multiplications by zero
                argument += NUTATION_ARG_TABLE[i][j] * arguments[j]
        coeff = NUTATION_SINE_COEF_TABLE[i][0]
        if NUTATION_SINE_COEF_TABLE[i][1]:
            coeff += NUTATION_SINE_COEF_TABLE[i][1] * t
        deltapsi += (coeff*sin(argument.rad()))/10000.0
    return Angle(0, 0, deltapsi)


def nutation_obliquity(*args, **kwargs):
    """This function computes the nutation in obliquity (Delta epsilon) at
    the provided date.

    This function internally uses an :class:`Epoch` object, and the **utc**
    argument then controls the way the UTC->TT conversion is handled for that
    object. If **leap_seconds** argument is set to a value different than zero,
    then that value will be used for the UTC->TAI conversion, and the internal
    leap seconds table will be bypassed.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param utc: Whether the provided epoch is a civil time (UTC) or TT
    :type utc: bool
    :param leap_seconds: This is the value to be used in the UTC->TAI
        conversion, instead of taking it from internal leap seconds table.
    :type leap_seconds: int, float

    :returns: The nutation in obliquity (Delta epsilon), as an
        :class:`Angle`
    :rtype: :class:`Angle`
    :raises: ValueError if input values are in the wrong range.
    :raises: TypeError if input values are of wrong type.

    >>> depsilon = nutation_obliquity(1987, 4, 10)
    >>> a = depsilon.dms_tuple()
    >>> a[0]
    0
    >>> a[1]
    0
    >>> round(a[2], 3)
    9.443
    >>> a[3]
    1.0
    """

    # Get the Epoch object corresponding to input parameters
    t = Epoch.check_input_date(*args, **kwargs)
    # Let's redefine t in units of Julian centuries from Epoch J2000.0
    t = (t.jde() - 2451545.0)/36525.0
    # Let's compute the mean elongation of the Moon from the Sun
    d = 297.85036 + t*(445267.111480 + t*(-0.0019142 + t/189474.0))
    d = Angle(d)            # Convert into an Angle: It is easier to handle
    # Compute the mean anomaly of the Sun (from Earth)
    m = 357.52772 + t*(35999.050340 + t*(-0.0001603 - t/300000.0))
    m = Angle(m)
    # Compute the mean anomaly of the Moon
    mprime = 134.96298 + t*(477198.867398 + t*(0.0086972 + t/56250.0))
    mprime = Angle(mprime)
    # Now, let's compute the Moon's argument of latitude
    f = 93.27191 + t*(483202.017538 + t*(-0.0036825 + t/327270.0))
    f = Angle(f)
    # And finally, the longitude of the ascending node of the Moon's mean
    # orbit on the ecliptic, measured from the mean equinox of date
    omega = 125.04452 + t*(-1934.136261 + t*(0.0020708 + t/450000.0))
    omega = Angle(omega)
    # Let's store this results in a list, in preparation for using tables
    arguments = [d, m, mprime, f, omega]
    # Now is time of using the nutation tables
    deltaepsilon = 0.0
    for i in range(len(NUTATION_COSINE_COEF_TABLE)):
        argument = Angle()
        coeff = 0.0
        for j in range(5):
            if NUTATION_ARG_TABLE[i][j]:    # Avoid multiplications by zero
                argument += NUTATION_ARG_TABLE[i][j] * arguments[j]
        coeff = NUTATION_COSINE_COEF_TABLE[i][0]
        if NUTATION_COSINE_COEF_TABLE[i][1]:
            coeff += NUTATION_COSINE_COEF_TABLE[i][1] * t
        deltaepsilon += (coeff*cos(argument.rad()))/10000.0
    return Angle(0, 0, deltaepsilon)


def precession_equatorial(start_epoch, final_epoch, start_ra, start_dec,
                          p_motion_ra=0.0, p_motion_dec=0.0):
    """This function converts the equatorial coordinates (right ascension and
    declination) given for an epoch and a equinox, to the corresponding
    values for another epoch and equinox. Only the **mean** positions, i.e.
    the effects of precession and proper motion, are considered here.

    :param start_epoch: Initial epoch when initial coordinates are given
    :type start_epoch: :py:class:`Epoch`
    :param final_epoch: Final epoch for when coordinates are going to be
        computed
    :type final_epoch: :py:class:`Epoch`
    :param start_ra: Initial right ascension
    :type start_ra: :py:class:`Angle`
    :param start_dec: Initial declination
    :type start_dec: :py:class:`Angle`
    :param p_motion_ra: Proper motion in right ascension, in degrees per
        year. Zero by default.
    :type p_motion_ra: :py:class:`Angle`
    :param p_motion_dec: Proper motion in declination, in degrees per year.
        Zero by default.
    :type p_motion_dec: :py:class:`Angle`

    :returns: Equatorial coordinates (right ascension, declination, in that
        order) corresponding to the final epoch, given as two objects
        :class:`Angle` inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> start_epoch = JDE2000
    >>> final_epoch = Epoch(2028, 11, 13.19)
    >>> alpha0 = Angle(2, 44, 11.986, ra=True)
    >>> delta0 = Angle(49, 13, 42.48)
    >>> pm_ra = Angle(0, 0, 0.03425, ra=True)
    >>> pm_dec = Angle(0, 0, -0.0895)
    >>> alpha, delta = precession_equatorial(start_epoch, final_epoch, alpha0,
    ...                                      delta0, pm_ra, pm_dec)
    >>> print(alpha.ra_str(False, 3))
    2:46:11.331
    >>> print(delta.dms_str(False, 2))
    49:20:54.54
    """

    # First check that input values are of correct types
    if not(isinstance(start_epoch, Epoch) and
           isinstance(final_epoch, Epoch) and
           isinstance(start_ra, Angle) and
           isinstance(start_dec, Angle)):
        raise TypeError("Invalid input types")
    if isinstance(p_motion_ra, (int, float)):
        p_motion_ra = Angle(p_motion_ra)
    if isinstance(p_motion_dec, (int, float)):
        p_motion_dec = Angle(p_motion_dec)
    if not (isinstance(p_motion_ra, Angle) and
            isinstance(p_motion_dec, Angle)):
        raise TypeError("Invalid input types")
    tt = (start_epoch - JDE2000)/36525.0
    t = (final_epoch - start_epoch)/36525.0
    # Correct starting coordinates by proper motion
    start_ra += p_motion_ra*t*100.0
    start_dec += p_motion_dec*t*100.0
    # Compute the conversion parameters
    zeta = t*((2306.2181 + tt*(1.39656 - 0.000139*tt)) +
              t*((0.30188 - 0.000344*tt) + 0.017998*t))
    z = t*((2306.2181 + tt*(1.39656 - 0.000139*tt)) +
           t*((1.09468 + 0.000066*tt) + 0.018203*t))
    theta = t*(2004.3109 + tt*(-0.85330 - 0.000217*tt) +
               t*(-(0.42665 + 0.000217*tt) - 0.041833*t))
    # Redefine the former values as Angles
    zeta = Angle(0, 0, zeta)
    z = Angle(0, 0, z)
    theta = Angle(0, 0, theta)
    a = cos(start_dec.rad())*sin(start_ra.rad() + zeta.rad())
    b = cos(theta.rad()) * cos(start_dec.rad()) * \
        cos(start_ra.rad() + zeta.rad()) - \
        sin(theta.rad()) * sin(start_dec.rad())
    c = sin(theta.rad()) * cos(start_dec.rad()) * \
        cos(start_ra.rad() + zeta.rad()) + \
        cos(theta.rad()) * sin(start_dec.rad())
    final_ra = atan2(a, b) + z.rad()
    if start_dec > 85.0:        # Coordinates are close to the pole
        final_dec = sqrt(a*a + b*b)
    else:
        final_dec = asin(c)
    # Convert results to Angles. Please note results are in radians
    final_ra = Angle(final_ra, radians=True)
    final_dec = Angle(final_dec, radians=True)
    return (final_ra, final_dec)


def precession_ecliptical(start_epoch, final_epoch, start_lon, start_lat,
                          p_motion_lon=0.0, p_motion_lat=0.0):
    """This function converts the ecliptical coordinates (longitude and
    latitude) given for an epoch and a equinox, to the corresponding
    values for another epoch and equinox. Only the **mean** positions, i.e.
    the effects of precession and proper motion, are considered here.

    :param start_epoch: Initial epoch when initial coordinates are given
    :type start_epoch: :py:class:`Epoch`
    :param final_epoch: Final epoch for when coordinates are going to be
        computed
    :type final_epoch: :py:class:`Epoch`
    :param start_lon: Initial longitude
    :type start_lon: :py:class:`Angle`
    :param start_lat: Initial latitude
    :type start_lat: :py:class:`Angle`
    :param p_motion_lon: Proper motion in longitude, in degrees per year.
        Zero by default.
    :type p_motion_lon: :py:class:`Angle`
    :param p_motion_lat: Proper motion in latitude, in degrees per year.
        Zero by default.
    :type p_motion_lat: :py:class:`Angle`

    :returns: Ecliptical coordinates (longitude, latitude, in that order)
        corresponding to the final epoch, given as two :class:`Angle`
        objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> start_epoch = JDE2000
    >>> final_epoch = Epoch(-214, 6, 30.0)
    >>> lon0 = Angle(149.48194)
    >>> lat0 = Angle(1.76549)
    >>> lon, lat = precession_ecliptical(start_epoch, final_epoch, lon0, lat0)
    >>> print(round(lon(), 3))
    118.704
    >>> print(round(lat(), 3))
    1.615
    """

    # First check that input values are of correct types
    if not(isinstance(start_epoch, Epoch) and
           isinstance(final_epoch, Epoch) and
           isinstance(start_lon, Angle) and
           isinstance(start_lat, Angle)):
        raise TypeError("Invalid input types")
    if isinstance(p_motion_lon, (int, float)):
        p_motion_lon = Angle(p_motion_lon)
    if isinstance(p_motion_lat, (int, float)):
        p_motion_lat = Angle(p_motion_lat)
    if not (isinstance(p_motion_lon, Angle) and
            isinstance(p_motion_lat, Angle)):
        raise TypeError("Invalid input types")
    tt = (start_epoch - JDE2000)/36525.0
    t = (final_epoch - start_epoch)/36525.0
    # Correct starting coordinates by proper motion
    start_lon += p_motion_lon*t*100.0
    start_lat += p_motion_lat*t*100.0
    # Compute the conversion parameters
    eta = t*((47.0029 + tt*(-0.06603 + 0.000598*tt))
             + t*((-0.03302 + 0.000598*tt) + 0.00006*t))
    pi = tt*(3289.4789 + 0.60622*tt) + t*(-(869.8089 + 0.50491*tt)
                                          + 0.03536*t)
    p = t*(5029.0966 + tt*(2.22226 - 0.000042*tt)
           + t*(1.11113 - 0.000042*tt - 0.000006*t))
    eta = Angle(0, 0, eta)
    pi = Angle(0, 0, pi)
    p = Angle(0, 0, p)
    # But beware!: There is still a missing constant for pi. We didn't add
    # it before because of the mismatch between degrees and seconds
    pi += 174.876384
    a = cos(eta.rad()) * cos(start_lat.rad()) \
        * sin(pi.rad() - start_lon.rad()) \
        - sin(eta.rad()) * sin(start_lat.rad())
    b = cos(start_lat.rad()) * cos(pi.rad() - start_lon.rad())
    c = cos(eta.rad()) * sin(start_lat.rad()) + sin(eta.rad()) \
        * cos(start_lat.rad()) * sin(pi.rad() - start_lon.rad())
    final_lon = p.rad() + pi.rad() - atan2(a, b)
    final_lat = asin(c)
    # Convert results to Angles. Please note results are in radians
    final_lon = Angle(final_lon, radians=True)
    final_lat = Angle(final_lat, radians=True)
    return (final_lon, final_lat)


def p_motion_equa2eclip(p_motion_ra, p_motion_dec, ra, dec, lat, epsilon):
    """It is usual that proper motions are given in equatorial coordinates,
    not in ecliptical ones. Therefore, this function converts the provided
    proper motions in equatorial coordinates to the corresponding ones in
    ecliptical coordinates.

    :param p_motion_ra: Proper motion in right ascension, in degrees per
        year, as an :class:`Angle` object
    :type p_motion_ra: :py:class:`Angle`
    :param p_motion_dec: Proper motion in declination, in degrees per year,
        as an :class:`Angle` object
    :type p_motion_dec: :py:class:`Angle`
    :param ra: Right ascension of the astronomical object, as degrees in an
        :class:`Angle` object
    :type ra: :py:class:`Angle`
    :param dec: Declination of the astronomical object, as degrees in an
        :class:`Angle` object
    :type dec: :py:class:`Angle`
    :param lat: Ecliptical latitude of the astronomical object, as degrees
        in an :class:`Angle` object
    :type lat: :py:class:`Angle`
    :param epsilon: Obliquity of the ecliptic
    :type epsilon: :py:class:`Angle`

    :returns: Proper motions in ecliptical longitude and latitude (in that
        order), given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.
    """

    # First check that input values are of correct types
    if not(isinstance(p_motion_ra, Angle) and
           isinstance(p_motion_dec, Angle) and
           isinstance(ra, Angle) and isinstance(dec, Angle) and
           isinstance(lat, Angle) and isinstance(epsilon, Angle)):
        raise TypeError("Invalid input types")
    pm_ra = p_motion_ra.rad()
    pm_dec = p_motion_dec.rad()
    s_eps = sin(epsilon.rad())
    c_eps = cos(epsilon.rad())
    s_ra = sin(ra.rad())
    c_ra = cos(ra.rad())
    s_dec = sin(dec.rad())
    c_dec = cos(dec.rad())
    c_lat = cos(lat.rad())
    se_ca = s_eps*c_ra
    se_sd_sa = s_eps*s_dec*s_ra
    pa_cd = pm_ra*c_dec
    ce_cd = c_eps*c_dec
    cl2 = c_lat*c_lat
    p_motion_lon = (pm_dec*se_ca + pa_cd*(ce_cd + se_sd_sa)) / cl2
    p_motion_lat = (pm_dec*(ce_cd + se_sd_sa) - pa_cd*se_ca) / c_lat
    return (p_motion_lon, p_motion_lat)


def precession_newcomb(start_epoch, final_epoch, start_ra, start_dec,
                       p_motion_ra=0.0, p_motion_dec=0.0):
    """This function implements the Newcomb precessional equations used in
    the old FK4 system. It takes equatorial coordinates (right ascension
    and declination) given for an epoch and a equinox, and converts them to
    the corresponding values for another epoch and equinox. Only the
    **mean** positions, i.e. the effects of precession and proper motion,
    are considered here.

    :param start_epoch: Initial epoch when initial coordinates are given
    :type start_epoch: :py:class:`Epoch`
    :param final_epoch: Final epoch for when coordinates are going to be
        computed
    :type final_epoch: :py:class:`Epoch`
    :param start_ra: Initial right ascension
    :type start_ra: :py:class:`Angle`
    :param start_dec: Initial declination
    :type start_dec: :py:class:`Angle`
    :param p_motion_ra: Proper motion in right ascension, in degrees per
        year. Zero by default.
    :type p_motion_ra: :py:class:`Angle`
    :param p_motion_dec: Proper motion in declination, in degrees per year.
        Zero by default.
    :type p_motion_dec: :py:class:`Angle`

    :returns: Equatorial coordinates (right ascension, declination, in that
        order) corresponding to the final epoch, given as two objects
        :class:`Angle` inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.
    """

    # First check that input values are of correct types
    if not(isinstance(start_epoch, Epoch) and
           isinstance(final_epoch, Epoch) and
           isinstance(start_ra, Angle) and
           isinstance(start_dec, Angle)):
        raise TypeError("Invalid input types")
    if isinstance(p_motion_ra, (int, float)):
        p_motion_ra = Angle(p_motion_ra)
    if isinstance(p_motion_dec, (int, float)):
        p_motion_dec = Angle(p_motion_dec)
    if not (isinstance(p_motion_ra, Angle) and
            isinstance(p_motion_dec, Angle)):
        raise TypeError("Invalid input types")
    tt = (start_epoch - 2415020.3135)/36524.2199
    t = (final_epoch - start_epoch)/36524.2199
    # Correct starting coordinates by proper motion
    start_ra += p_motion_ra*t*100.0
    start_dec += p_motion_dec*t*100.0
    # Compute the conversion parameters
    zeta = t*(2304.25 + 1.396*tt + t*(0.302 + 0.018*t))
    z = zeta + t*t*(0.791 + 0.001*t)
    theta = t*(2004.682 - 0.853*tt - t*(0.426 + 0.042*t))
    # Redefine the former values as Angles
    zeta = Angle(0, 0, zeta)
    z = Angle(0, 0, z)
    theta = Angle(0, 0, theta)
    a = cos(start_dec.rad())*sin(start_ra.rad() + zeta.rad())
    b = cos(theta.rad()) * cos(start_dec.rad()) * \
        cos(start_ra.rad() + zeta.rad()) - \
        sin(theta.rad()) * sin(start_dec.rad())
    c = sin(theta.rad()) * cos(start_dec.rad()) * \
        cos(start_ra.rad() + zeta.rad()) + \
        cos(theta.rad()) * sin(start_dec.rad())
    final_ra = atan2(a, b) + z.rad()
    if start_dec > 85.0:        # Coordinates are close to the pole
        final_dec = sqrt(a*a + b*b)
    else:
        final_dec = asin(c)
    # Convert results to Angles. Please note results are in radians
    final_ra = Angle(final_ra, radians=True)
    final_dec = Angle(final_dec, radians=True)
    return (final_ra, final_dec)


def motion_in_space(start_ra, start_dec, distance, velocity,
                    p_motion_ra, p_motion_dec, time):
    """This function computes the star's true motion through space relative
    to the Sun, allowing to compute the start proper motion at a given
    time.

    :param start_ra: Initial right ascension
    :type start_ra: :py:class:`Angle`
    :param start_dec: Initial declination
    :type start_dec: :py:class:`Angle`
    :param distance: Star's distance to the Sun, in parsecs. If distance is
        given in light-years, multipy it by 0.3066. If the star's parallax
        **pi** (in arcseconds) is given, use (1.0/pi).
    :type distance: float
    :param velocity: Radial velocity in km/s
    :type velocity: float
    :param p_motion_ra: Proper motion in right ascension, in degrees per
        year.
    :type p_motion_ra: :py:class:`Angle`
    :param p_motion_dec: Proper motion in declination, in degrees per year.
    :type p_motion_dec: :py:class:`Angle`
    :param time: Number of years since starting epoch, positive in the
        future, negative in the past
    :type time: float

    :returns: Equatorial coordinates (right ascension, declination, in that
        order) corresponding to the final epoch, given as two objects
        :class:`Angle` inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> ra = Angle(6, 45, 8.871, ra=True)
    >>> dec = Angle(-16.716108)
    >>> pm_ra = Angle(0, 0, -0.03847, ra=True)
    >>> pm_dec = Angle(0, 0, -1.2053)
    >>> dist = 2.64
    >>> vel = -7.6
    >>> alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec,
    ...                                -1000.0)
    >>> print(alpha.ra_str(False, 2))
    6:45:47.16
    >>> print(delta.dms_str(False, 1))
    -16:22:56.0
    >>> alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec,
    ...                                -4000.0)
    >>> print(alpha.ra_str(False, 2))
    6:47:39.91
    >>> print(delta.dms_str(False, 1))
    -15:23:30.6
    """
    # >>> ra = Angle(101.286962)

    # First check that input values are of correct types
    if not(isinstance(start_ra, Angle) and
           isinstance(start_dec, Angle)):
        raise TypeError("Invalid input types")
    if isinstance(p_motion_ra, (int, float)):
        p_motion_ra = Angle(p_motion_ra)
    if isinstance(p_motion_dec, (int, float)):
        p_motion_dec = Angle(p_motion_dec)
    if not (isinstance(p_motion_ra, Angle) and
            isinstance(p_motion_dec, Angle)):
        raise TypeError("Invalid input types")
    if not(isinstance(distance, (int, float)) and
           isinstance(velocity, (int, float)) and
           isinstance(time, (int, float))):
        raise TypeError("Invalid input types")
    dr = velocity/977792.0
    x = distance*cos(start_dec.rad())*cos(start_ra.rad())
    y = distance*cos(start_dec.rad())*sin(start_ra.rad())
    z = distance*sin(start_dec.rad())
    dx = (x/distance)*dr - z*p_motion_dec.rad()*cos(start_ra.rad()) \
        - y*p_motion_ra.rad()
    dy = (y/distance)*dr - z*p_motion_dec.rad()*sin(start_ra.rad()) \
        + x*p_motion_ra.rad()
    dz = (z/distance)*dr + distance*p_motion_dec.rad()*cos(start_dec.rad())
    xp = x + time*dx
    yp = y + time*dy
    zp = z + time*dz
    final_ra = atan2(yp, xp)
    final_dec = atan(zp/sqrt(xp*xp + yp*yp))
    # Convert results to Angles. Please note results are in radians
    final_ra = Angle(final_ra, radians=True)
    final_dec = Angle(final_dec, radians=True)
    return (final_ra, final_dec)


def equatorial2ecliptical(right_ascension, declination, obliquity):
    """This function converts from equatorial coordinated (right ascension and
    declination) to ecliptical coordinates (longitude and latitude).

    :param right_ascension: Right ascension, as an Angle object
    :type start_epoch: :py:class:`Angle`
    :param declination: Declination, as an Angle object
    :type start_epoch: :py:class:`Angle`
    :param obliquity: Obliquity of the ecliptic, as an Angle object
    :type obliquity: :py:class:`Angle`

    :returns: Ecliptical coordinates (longitude, latitude, in that order),
        given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> ra = Angle(7, 45, 18.946, ra=True)
    >>> dec = Angle(28, 1, 34.26)
    >>> epsilon = Angle(23.4392911)
    >>> lon, lat = equatorial2ecliptical(ra, dec, epsilon)
    >>> print(round(lon(), 5))
    113.21563
    >>> print(round(lat(), 5))
    6.68417
    """

    # First check that input values are of correct types
    if not(isinstance(right_ascension, Angle) and
           isinstance(declination, Angle) and
           isinstance(obliquity, Angle)):
        raise TypeError("Invalid input types")
    ra = right_ascension.rad()
    dec = declination.rad()
    eps = obliquity.rad()
    lon = atan2((sin(ra)*cos(eps) + tan(dec)*sin(eps)), cos(ra))
    lat = asin(sin(dec)*cos(eps) - cos(dec)*sin(eps)*sin(ra))
    lon = Angle(lon, radians=True)
    lat = Angle(lat, radians=True)
    return (lon, lat)


def ecliptical2equatorial(longitude, latitude, obliquity):
    """This function converts from ecliptical coordinates (longitude and
    latitude) to equatorial coordinated (right ascension and declination).

    :param longitude: Ecliptical longitude, as an Angle object
    :type longitude: :py:class:`Angle`
    :param latitude: Ecliptical latitude, as an Angle object
    :type latitude: :py:class:`Angle`
    :param obliquity: Obliquity of the ecliptic, as an Angle object
    :type obliquity: :py:class:`Angle`

    :returns: Equatorial coordinates (right ascension, declination, in that
        order), given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> lon = Angle(113.21563)
    >>> lat = Angle(6.68417)
    >>> epsilon = Angle(23.4392911)
    >>> ra, dec = ecliptical2equatorial(lon, lat, epsilon)
    >>> print(ra.ra_str(n_dec=3))
    7h 45' 18.946''
    >>> print(dec.dms_str(n_dec=2))
    28d 1' 34.26''
    """

    # First check that input values are of correct types
    if not(isinstance(longitude, Angle) and
           isinstance(latitude, Angle) and
           isinstance(obliquity, Angle)):
        raise TypeError("Invalid input types")
    lon = longitude.rad()
    lat = latitude.rad()
    eps = obliquity.rad()
    ra = atan2((sin(lon)*cos(eps) - tan(lat)*sin(eps)), cos(lon))
    dec = asin(sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon))
    ra = Angle(ra, radians=True)
    dec = Angle(dec, radians=True)
    return (ra, dec)


def equatorial2horizontal(hour_angle, declination, geo_latitude):
    """This function converts from equatorial coordinates (right ascension and
    declination) to local horizontal coordinates (azimuth and elevation).

    Following Meeus' convention, the azimuth is measured westward from the
    SOUTH. If you want the azimuth to be measured from the north (common custom
    between navigators and meteorologits), you should add 180 degrees.

    The hour angle (H) comprises information about the sidereal time, the
    observer's geodetic longitude (positive west from Greenwich) and the right
    ascension. If theta is the local sidereal time, theta0 the sidereal time at
    Greenwich, lon the observer's longitude and ra the right ascension, the
    following expressions hold:

        H = theta - ra
        H = theta0 - lon - ra

    :param hour_angle: Hour angle, as an Angle object
    :type hour_angle: :py:class:`Angle`
    :param declination: Declination, as an Angle object
    :type declination: :py:class:`Angle`
    :param geo_latitude: Geodetic latitude of the observer, as an Angle object
    :type geo_latitude: :py:class:`Angle`

    :returns: Local horizontal coordinates (azimuth, elevation, in that order),
        given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> lon = Angle(77, 3, 56)
    >>> lat = Angle(38, 55, 17)
    >>> ra = Angle(23, 9, 16.641, ra=True)
    >>> dec = Angle(-6, 43, 11.61)
    >>> theta0 = Angle(8, 34, 57.0896, ra=True)
    >>> eps = Angle(23, 26, 36.87)
    >>> delta = Angle(0, 0, ((-3.868*cos(eps.rad()))/15.0), ra=True)
    >>> theta0 += delta
    >>> h = theta0 - lon - ra
    >>> azi, ele = equatorial2horizontal(h, dec, lat)
    >>> print(round(azi, 3))
    68.034
    >>> print(round(ele, 3))
    15.125
    """

    # First check that input values are of correct types
    if not(isinstance(hour_angle, Angle) and
           isinstance(declination, Angle) and
           isinstance(geo_latitude, Angle)):
        raise TypeError("Invalid input types")
    h = hour_angle.rad()
    dec = declination.rad()
    lat = geo_latitude.rad()
    azi = atan2(sin(h), (cos(h)*sin(lat) - tan(dec)*cos(lat)))
    ele = asin(sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(h))
    azi = Angle(azi, radians=True)
    ele = Angle(ele, radians=True)
    return (azi, ele)


def horizontal2equatorial(azimuth, elevation, geo_latitude):
    """This function converts from local horizontal coordinates (azimuth and
    elevation) to equatorial coordinates (right ascension and declination).

    Following Meeus' convention, the azimuth is measured westward from the
    SOUTH.

    This function returns the hour angle and the declination. The hour angle
    (H) comprises information about the sidereal time, the observer's geodetic
    longitude (positive west from Greenwich) and the right ascension. If theta
    is the local sidereal time, theta0 the sidereal time at Greenwich, lon the
    observer's longitude and ra the right ascension, the following expressions
    hold:

        H = theta - ra
        H = theta0 - lon - ra

    :param azimuth: Azimuth, measured westward from south, as an Angle object
    :type azimuth: :py:class:`Angle`
    :param elevation: Elevation from the horizon, as an Angle object
    :type elevation: :py:class:`Angle`
    :param geo_latitude: Geodetic latitude of the observer, as an Angle object
    :type geo_latitude: :py:class:`Angle`

    :returns: Equatorial coordinates (as hour angle and declination, in that
        order), given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> azi = Angle(68.0337)
    >>> ele = Angle(15.1249)
    >>> lat = Angle(38, 55, 17)
    >>> h, dec = horizontal2equatorial(azi, ele, lat)
    >>> print(round(h, 4))
    64.3521
    >>> print(dec.dms_str(n_dec=0))
    -6d 43' 12.0''
    """

    # First check that input values are of correct types
    if not(isinstance(azimuth, Angle) and
           isinstance(elevation, Angle) and
           isinstance(geo_latitude, Angle)):
        raise TypeError("Invalid input types")
    azi = azimuth.rad()
    ele = elevation.rad()
    lat = geo_latitude.rad()
    h = atan2(sin(azi), (cos(azi)*sin(lat) + tan(ele)*cos(lat)))
    dec = asin(sin(lat)*sin(ele) - cos(lat)*cos(ele)*cos(azi))
    h = Angle(h, radians=True)
    dec = Angle(dec, radians=True)
    return (h, dec)


def equatorial2galactic(right_ascension, declination):
    """This function converts from equatorial coordinates (right ascension and
    declination) to galactic coordinates (longitude and latitude).

    The current galactic system of coordinates was defined by the International
    Astronomical Union in 1959, using the standard equatorial system of epoch
    B1950.0.

    :param right_ascension: Right ascension, as an Angle object
    :type right_ascension: :py:class:`Angle`
    :param declination: Declination, as an Angle object
    :type declination: :py:class:`Angle`

    :returns: Galactic coordinates (longitude and latitude, in that order),
        given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> ra = Angle(17, 48, 59.74, ra=True)
    >>> dec = Angle(-14, 43, 8.2)
    >>> lon, lat = equatorial2galactic(ra, dec)
    >>> print(round(lon, 4))
    12.9593
    >>> print(round(lat, 4))
    6.0463
    """

    # First check that input values are of correct types
    if not(isinstance(right_ascension, Angle) and
           isinstance(declination, Angle)):
        raise TypeError("Invalid input types")
    ra = right_ascension.rad()
    dec = declination.rad()
    c1 = Angle(192.25)
    c1 = c1.rad()
    c1ra = c1 - ra
    c2 = Angle(27.4)
    c2 = c2.rad()
    x = atan2(sin(c1ra), (cos(c1ra)*sin(c2) - tan(dec)*cos(c2)))
    lon = Angle(-x, radians=True)
    lon = 303.0 + lon
    lat = asin(sin(dec)*sin(c2) + cos(dec)*cos(c2)*cos(c1ra))
    lat = Angle(lat, radians=True)
    return (lon, lat)


def galactic2equatorial(longitude, latitude):
    """This function converts from galactic coordinates (longitude and latitude)
    to equatorial coordinates (right ascension and declination).

    The current galactic system of coordinates was defined by the International
    Astronomical Union in 1959, using the standard equatorial system of epoch
    B1950.0.

    :param longitude: Longitude, as an Angle object
    :type longitude: :py:class:`Angle`
    :param latitude: Latitude, as an Angle object
    :type latitude: :py:class:`Angle`

    :returns: Equatorial coordinates (right ascension and declination, in that
        order), given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> lon = Angle(12.9593)
    >>> lat = Angle(6.0463)
    >>> ra, dec = galactic2equatorial(lon, lat)
    >>> print(ra.ra_str(n_dec=1))
    17h 48' 59.7''
    >>> print(dec.dms_str(n_dec=0))
    -14d 43' 8.0''
    """

    # First check that input values are of correct types
    if not(isinstance(longitude, Angle) and
           isinstance(latitude, Angle)):
        raise TypeError("Invalid input types")
    lon = longitude.rad()
    lat = latitude.rad()
    c1 = Angle(123.0)
    c1 = c1.rad()
    c2 = Angle(27.4)
    c2 = c2.rad()
    lc1 = lon - c1
    y = atan2(sin(lc1), (cos(lc1)*sin(c2) - tan(lat)*cos(c2)))
    y = Angle(y, radians=True)
    ra = y + 12.25
    ra.to_positive()
    dec = asin(sin(lat)*sin(c2) + cos(lat)*cos(c2)*cos(lc1))
    dec = Angle(dec, radians=True)
    return (ra, dec)


def parallactic_angle(hour_angle, declination, geo_latitude):
    """This function computes the parallactic angle, an apparent rotation that
    appears because celestial bodies move along parallel circles. By
    convention, the parallactic angle is negative before the passage through
    the southern meridian (in the north hemisphere), and positive afterwards.
    Exactly on the meridian, its value is zero.

    Please note that when the celestial body is exactly at the zenith, the
    parallactic angle is not defined, and this function will return 'None'.

    The hour angle (H) comprises information about the sidereal time, the
    observer's geodetic longitude (positive west from Greenwich) and the right
    ascension. If theta is the local sidereal time, theta0 the sidereal time at
    Greenwich, lon the observer's longitude and ra the right ascension, the
    following expressions hold:

        H = theta - ra
        H = theta0 - lon - ra

    :param hour_angle: Hour angle, as an Angle object
    :type hour_angle: :py:class:`Angle`
    :param declination: Declination, as an Angle object
    :type declination: :py:class:`Angle`
    :param geo_latitude: Geodetic latitude of the observer, as an Angle object
    :type geo_latitude: :py:class:`Angle`

    :returns: Local horizontal coordinates (azimuth, elevation, in that order),
        given as two :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.
    """

    # First check that input values are of correct types
    if not(isinstance(hour_angle, Angle) and
           isinstance(declination, Angle) and
           isinstance(geo_latitude, Angle)):
        raise TypeError("Invalid input types")
    h = hour_angle.rad()
    dec = declination.rad()
    lat = geo_latitude.rad()
    den = tan(lat)*cos(dec) - sin(dec)*cos(h)
    if abs(den) < TOL:
        return None
    q = atan2(sin(h), den)
    q = Angle(q, radians=True)
    return q


def ecliptic_horizon(local_sidereal_time, geo_latitude, obliquity):
    """This function returns the longitudes of the two points of the ecliptic
    which are on the horizon, as well as the angle between the ecliptic and the
    horizon.

    :param local_sidereal_time: Local sidereal time, as an Angle object
    :type local_sidereal_time: :py:class:`Angle`
    :param geo_latitude: Geodetic latitude, as an Angle object
    :type geo_latitude: :py:class:`Angle`
    :param obliquity: Obliquity of the ecliptic, as an Angle object
    :type obliquity: :py:class:`Angle`

    :returns: Longitudes of the two points of the ecliptic which are on the
        horizon, and the angle between the ecliptic and the horizon (in that
        order), given as three :class:`Angle` objects inside a tuple
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> sidereal_time = Angle(5.0, ra=True)
    >>> lat = Angle(51.0)
    >>> epsilon = Angle(23.44)
    >>> lon1, lon2, i = ecliptic_horizon(sidereal_time, lat, epsilon)
    >>> print(lon1.dms_str(n_dec=1))
    169d 21' 29.9''
    >>> print(lon2.dms_str(n_dec=1))
    349d 21' 29.9''
    >>> print(round(i, 0))
    62.0
    """

    # First check that input values are of correct types
    if not(isinstance(local_sidereal_time, Angle) and
           isinstance(geo_latitude, Angle) and
           isinstance(obliquity, Angle)):
        raise TypeError("Invalid input types")
    theta = local_sidereal_time.rad()
    lat = geo_latitude.rad()
    eps = obliquity.rad()
    # First, let's compute the longitudes of the ecliptic points on the horizon
    lon1 = atan2(-cos(theta), (sin(eps)*tan(lat) + cos(eps)*sin(theta)))
    lon1 = Angle(lon1, radians=True)
    lon1.to_positive()
    # Get the second point, which is 180 degrees apart
    if lon1 < 180.0:
        lon2 = lon1 + 180.0
    else:
        lon2 = lon1
        lon1 = lon2 - 180.0
    # Now, compute the angle between the ecliptic and the horizon
    i = acos(cos(eps)*sin(lat) - sin(eps)*cos(lat)*sin(theta))
    i = Angle(i, radians=True)
    return (lon1, lon2, i)


def ecliptic_equator(longitude, latitude, obliquity):
    """This function returns the angle between the direction of the northern
    celestial pole and the direction of the north pole of the ecliptic, taking
    as reference the point whose ecliptic longitude and latitude are given.

    Please note that if we make latitude=0, the result is the angle between the
    ecliptic (at the given ecliptical longitude) and the east-west direction on
    the celestial sphere.

    :param longitude: Ecliptical longitude, as an Angle object
    :type longitude: :py:class:`Angle`
    :param latitude: Ecliptical latitude, as an Angle object
    :type latitude: :py:class:`Angle`
    :param obliquity: Obliquity of the ecliptic, as an Angle object
    :type obliquity: :py:class:`Angle`

    :returns: Angle between the direction of the northern celestial pole and
        the direction of the north pole of the ecliptic, given as one
        :class:`Angle` object
    :rtype: :class:`Angle`
    :raises: TypeError if input values are of wrong type.
    """

    # First check that input values are of correct types
    if not(isinstance(longitude, Angle) and
           isinstance(latitude, Angle) and
           isinstance(obliquity, Angle)):
        raise TypeError("Invalid input types")
    lon = longitude.rad()
    lat = latitude.rad()
    eps = obliquity.rad()
    q = atan2((cos(lon)*tan(eps)), (sin(lat)*sin(lon)*tan(eps) - cos(lat)))
    q = Angle(q, radians=True)
    return q


def diurnal_path_horizon(declination, geo_latitude):
    """This function returns the angle of the diurnal path of a celestial body
    relative to the horizon at the time of its rising or setting.

    :param declination: Declination, as an Angle object
    :type declination: :py:class:`Angle`
    :param geo_latitude: Geodetic latitude, as an Angle object
    :type geo_latitude: :py:class:`Angle`

    :returns: Angle of the diurnal path of the celestial body relative to the
        horizon at the time of rising or setting, given as one
        :class:`Angle` object
    :rtype: :class:`Angle`
    :raises: TypeError if input values are of wrong type.
    """

    # First check that input values are of correct types
    if not(isinstance(declination, Angle) and
           isinstance(geo_latitude, Angle)):
        raise TypeError("Invalid input types")
    dec = declination.rad()
    lat = geo_latitude.rad()
    b = tan(dec)*tan(lat)
    c = sqrt(1.0 - b*b)
    j = atan2(c*cos(dec), tan(lat))
    j = Angle(j, radians=True)
    return j


def times_rise_transit_set(longitude, latitude, alpha1, delta1, alpha2, delta2,
                           alpha3, delta3, h0, delta_t, theta0):
    """This function computes the times (in Universal Time UT) of rising,
    transit and setting of a given celestial body.

    .. note:: If the body is circumpolar there are no rising, transit nor
        setting times. In such a case a tuple with None's is returned

    .. note:: Care must be taken when interpreting the results. For instance,
        if the setting time is **smaller** than the rising time, it means that
        it belongs to the **following** day. Also, if the rising time is
        **bigger** than the setting time, it belong to the **previous** day.
        The same applies to the transit time.

    :param longitude: Geodetic longitude, as an Angle object. It is measured
        positively west from Greenwich, and negatively to the east.
    :type longitude: :py:class:`Angle`
    :param latitude: Geodetic latitude, as an Angle object
    :type latitude: :py:class:`Angle`
    :param alpha1: Apparent right ascension the previous day at 0h TT, as an
        Angle object
    :type alpha1: :py:class:`Angle`
    :param delta1: Apparent declination the previous day at 0h TT, as an Angle
        object
    :type delta1: :py:class:`Angle`
    :param alpha2: Apparent right ascension the current day at 0h TT, as an
        Angle object
    :type alpha2: :py:class:`Angle`
    :param delta2: Apparent declination the current day at 0h TT, as an Angle
        object
    :type delta2: :py:class:`Angle`
    :param alpha3: Apparent right ascension the following day at 0h TT, as an
        Angle object
    :type alpha3: :py:class:`Angle`
    :param delta3: Apparent declination the following day at 0h TT, as an Angle
        object
    :type delta3: :py:class:`Angle`
    :param h0: 'Standard' altitude: the geometric altitude of the center of the
        body at the time of apparent rising or setting, as degrees in an Angle
        object. It should be -0.5667 deg for stars and planets, -0.8333 deg
        for the Sun, and 0.125 deg for the Moon.
    :type h0: :py:class:`Angle`
    :param delta_t: The difference between Terrestrial Time and Universal Time
        (TT - UT) in seconds of time
    :type delta_t: float
    :param theta0: Apparent sidereal time at 0h TT on the current day for the
        meridian of Greenwich, as degrees in an Angle object
    :type theta0: :py:class:`Angle`

    :returns: A tuple with the times of rising, transit and setting, in that
        order, as hours in UT.
    :rtype: tuple
    :raises: TypeError if input values are of wrong type.

    >>> longitude = Angle(71, 5, 0.0)
    >>> latitude = Angle(42, 20, 0.0)
    >>> alpha1 = Angle(2, 42, 43.25, ra=True)
    >>> delta1 = Angle(18, 2, 51.4)
    >>> alpha2 = Angle(2, 46, 55.51, ra=True)
    >>> delta2 = Angle(18, 26, 27.3)
    >>> alpha3 = Angle(2, 51, 7.69, ra=True)
    >>> delta3 = Angle(18, 49, 38.7)
    >>> h0 = Angle(-0.5667)
    >>> delta_t = 56.0
    >>> theta0 = Angle(11, 50, 58.1, ra=True)
    >>> rising, transit, setting = times_rise_transit_set(longitude, latitude,\
                                                          alpha1, delta1, \
                                                          alpha2, delta2, \
                                                          alpha3, delta3, h0, \
                                                          delta_t, theta0)
    >>> print(round(rising, 4))
    12.4238
    >>> print(round(transit, 3))
    19.675
    >>> print(round(setting, 3))
    2.911
    """

    def check_value(m):
        while m < 0 or m > 1.0:
            if m < 0.0:
                m += 1
            elif m > 1.0:
                m -= 1
        return m

    def interpol(n, y1, y2, y3):
        a = y2 - y1
        b = y3 - y2
        c = b - a
        return y2 + n*(a + b + n*c)/2.0

    # First check that input values are of correct types
    if not(isinstance(longitude, Angle) and isinstance(latitude, Angle) and
           isinstance(alpha1, Angle) and isinstance(delta1, Angle) and
           isinstance(alpha2, Angle) and isinstance(delta2, Angle) and
           isinstance(alpha3, Angle) and isinstance(delta3, Angle) and
           isinstance(h0, Angle) and isinstance(theta0, Angle) and
           isinstance(delta_t, (int, float))):
        raise TypeError("Invalid input types")
    # Let's start computing approximate times
    h = h0.rad()
    lat = latitude.rad()
    d2 = delta2.rad()
    hh0 = (sin(h) - sin(lat)*sin(d2))/(cos(lat)*cos(d2))
    # Check if the body is circumpolar. In such case, there are no rising,
    # transit nor setting times, and a tuple with None's is returned
    if abs(hh0) > 1.0:
        return (None, None, None)
    hh0 = acos(hh0)
    hh0 = Angle(hh0, radians=True)
    hh0.to_positive()
    m0 = (alpha2 + longitude - theta0)/360.0
    m0 = m0()                               # m0 is an Angle. Convert to float
    m1 = m0 - hh0()/360.0
    m2 = m0 + hh0()/360.0
    m0 = check_value(m0)
    m1 = check_value(m1)
    m2 = check_value(m2)
    # Carry out this procedure twice
    for _ in range(2):
        # Interpolate alpha and delta values for each (m0, m1, m2)
        n = m0 + delta_t/86400.0
        transit_alpha = interpol(n, alpha1, alpha2, alpha3)
        n = m1 + delta_t/86400.0
        rise_alpha = interpol(n, alpha1, alpha2, alpha3)
        rise_delta = interpol(n, delta1, delta2, delta3)
        n = m2 + delta_t/86400.0
        set_alpha = interpol(n, alpha1, alpha2, alpha3)
        set_delta = interpol(n, delta1, delta2, delta3)
        # Compute the hour angles
        theta = theta0 + 360.985647*m0
        transit_ha = theta - longitude - transit_alpha
        delta_transit = transit_ha/(-360.0)
        theta = theta0 + 360.985647*m1
        rise_ha = theta - longitude - rise_alpha
        theta = theta0 + 360.985647*m2
        set_ha = theta - longitude - set_alpha
        # We need the elevations
        azi, rise_ele = equatorial2horizontal(rise_ha, rise_delta, latitude)
        azi, set_ele = equatorial2horizontal(set_ha, set_delta, latitude)
        delta_rise = (rise_ele - h0)/(360.0 * cos(rise_delta.rad()) *
                                      cos(lat) * sin(rise_ha.rad()))
        delta_set = (set_ele - h0)/(360.0 * cos(set_delta.rad()) *
                                    cos(lat) * sin(set_ha.rad()))
        m0 += delta_transit()
        m1 += delta_rise()
        m2 += delta_set()
    return (m1*24.0, m0*24.0, m2*24.0)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Earth class
    print('\n' + 35*'*')
    print("*** Use of Coordinate functions")
    print(35*'*' + '\n')

    # Here follows a series of important parameters related to the angle
    # between Earth's rotation axis and the ecliptic
    e0 = mean_obliquity(1987, 4, 10)
    print("The mean angle between Earth rotation axis and ecliptic axis for " +
          "1987/4/10 is:")
    print_me("Mean obliquity", e0.dms_str(n_dec=3))         # 23d 26' 27.407''
    epsilon = true_obliquity(1987, 4, 10)
    print("'True' (instantaneous) angle between those axes for 1987/4/10 is:")
    print_me("True obliquity", epsilon.dms_str(n_dec=3))    # 23d 26' 36.849''
    epsilon = true_obliquity(2018, 7, 29)
    print("'True' (instantaneous) angle between those axes for 2018/7/29 is:")
    print_me("True obliquity", epsilon.dms_str(True, 4))    # 23d 26' 7.2157''

    # The nutation effect is separated in two components: One parallel to the
    # ecliptic (nutation in longitude) and other perpendicular to the ecliptic
    # (nutation in obliquity)
    print("Nutation correction in longitude for 1987/4/10:")
    dpsi = nutation_longitude(1987, 4, 10)
    print_me("Nutation in longitude", dpsi.dms_str(n_dec=3))   # 0d 0' -3.788''
    print("Nutation correction in obliquity for 1987/4/10:")
    depsilon = nutation_obliquity(1987, 4, 10)            # 0d 0' 9.443''
    print_me("Nutation in obliquity", depsilon.dms_str(n_dec=3))

    print("")

    # We can compute the effects of precession on the equatorial coordinates of
    # a given star, taking also into account its proper motion

    start_epoch = JDE2000
    final_epoch = Epoch(2028, 11, 13.19)
    alpha0 = Angle(2, 44, 11.986, ra=True)
    delta0 = Angle(49, 13, 42.48)                             # 2h 44' 11.986''
    print_me("Initial right ascension", alpha0.ra_str(n_dec=3))
    print_me("Initial declination", delta0.dms_str(n_dec=2))  # 49d 13' 42.48''
    pm_ra = Angle(0, 0, 0.03425, ra=True)
    pm_dec = Angle(0, 0, -0.0895)
    alpha, delta = precession_equatorial(start_epoch, final_epoch, alpha0,
                                         delta0, pm_ra, pm_dec)
    print_me("Final right ascension", alpha.ra_str(n_dec=3))  # 2h 46' 11.331''
    print_me("Final declination", delta.dms_str(n_dec=2))     # 49d 20' 54.54''

    print("")

    # Something similar can also be done with the ecliptical coordinates
    start_epoch = JDE2000
    final_epoch = Epoch(-214, 6, 30.0)
    lon0 = Angle(149.48194)
    lat0 = Angle(1.76549)
    print_me("Initial ecliptical longitude", round(lon0(), 5))      # 149.48194
    print_me("Initial ecliptical latitude", round(lat0(), 5))       # 1.76549
    lon, lat = precession_ecliptical(start_epoch, final_epoch, lon0, lat0)
    print_me("Final ecliptical longitude", round(lon(), 3))         # 118.704
    print_me("Final ecliptical latitude", round(lat(), 3))          # 1.615

    print("")

    # It is possible to compute with relative accuracy the proper motion of the
    # stars, taking into account their distance to Sun and relative velocity
    ra = Angle(6, 45, 8.871, ra=True)
    dec = Angle(-16.716108)
    pm_ra = Angle(0, 0, -0.03847, ra=True)
    pm_dec = Angle(0, 0, -1.2053)
    dist = 2.64
    vel = -7.6
    alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -1000.0)

    print_me("Right ascension, year 2000", ra.ra_str(True, 2))
    print_me("Right ascension, year 1000", alpha.ra_str(True, 2))
    # 6h 45' 47.16''
    print_me("Declination, year 2000", dec.dms_str(True, 1))
    print_me("Declination, year 1000", delta.dms_str(True, 1))
    # -16d 22' 56.0''

    print("")

    # This module provides a series of functions to convert between equatorial,
    # ecliptical, horizontal and galactic coordinates

    ra = Angle(7, 45, 18.946, ra=True)
    dec = Angle(28, 1, 34.26)
    epsilon = Angle(23.4392911)
    lon, lat = equatorial2ecliptical(ra, dec, epsilon)
    print_me("Equatorial to ecliptical. Longitude", round(lon(), 5))
    # 113.21563
    print_me("Equatorial to ecliptical. Latitude", round(lat(), 5))
    # 6.68417

    print("")

    lon = Angle(113.21563)
    lat = Angle(6.68417)
    epsilon = Angle(23.4392911)
    ra, dec = ecliptical2equatorial(lon, lat, epsilon)
    print_me("Ecliptical to equatorial. Right ascension", ra.ra_str(n_dec=3))
    # 7h 45' 18.946''
    print_me("Ecliptical to equatorial. Declination", dec.dms_str(n_dec=2))
    # 28d 1' 34.26''

    print("")

    lon = Angle(77, 3, 56)
    lat = Angle(38, 55, 17)
    ra = Angle(23, 9, 16.641, ra=True)
    dec = Angle(-6, 43, 11.61)
    theta0 = Angle(8, 34, 57.0896, ra=True)
    eps = Angle(23, 26, 36.87)
    # Compute correction to convert from mean to apparent sidereal time
    delta = Angle(0, 0, ((-3.868*cos(eps.rad()))/15.0), ra=True)
    theta0 += delta
    h = theta0 - lon - ra
    azi, ele = equatorial2horizontal(h, dec, lat)
    print_me("Equatorial to horizontal: Azimuth", round(azi, 3))    # 68.034
    print_me("Equatorial to horizontal: Elevation", round(ele, 3))  # 15.125

    print("")

    azi = Angle(68.0337)
    ele = Angle(15.1249)
    lat = Angle(38, 55, 17)
    h, dec = horizontal2equatorial(azi, ele, lat)
    print_me("Horizontal to equatorial. Hour angle", round(h, 4))   # 64.3521
    print_me("Horizontal to equatorial. Declination", dec.dms_str(n_dec=0))
    # -6d 43' 12.0''

    print("")

    ra = Angle(17, 48, 59.74, ra=True)
    dec = Angle(-14, 43, 8.2)
    lon, lat = equatorial2galactic(ra, dec)
    print_me("Equatorial to galactic. Longitude", round(lon, 4))    # 12.9593
    print_me("Equatorial to galactic. Latitude", round(lat, 4))     # 6.0463

    print("")

    lon = Angle(12.9593)
    lat = Angle(6.0463)
    ra, dec = galactic2equatorial(lon, lat)
    print_me("Galactic to equatorial. Right ascension", ra.ra_str(n_dec=1))
    # 17h 48' 59.7''
    print_me("Galactic to equatorial. Declination", dec.dms_str(n_dec=0))
    # -14d 43' 8.0''

    print("")

    # Get the ecliptic longitudes of the two points of the ecliptic which are
    # on the horizon, as well as the angle between the ecliptic and the horizon
    sidereal_time = Angle(5.0, ra=True)
    lat = Angle(51.0)
    epsilon = Angle(23.44)
    lon1, lon2, i = ecliptic_horizon(sidereal_time, lat, epsilon)
    print_me("Longitude of ecliptic point #1 on the horizon",
             lon1.dms_str(n_dec=1))                         # 169d 21' 29.9''
    print_me("Longitude of ecliptic point #2 on the horizon",
             lon2.dms_str(n_dec=1))                         # 349d 21' 29.9''
    print_me("Angle between the ecliptic and the horizon", round(i, 0))  # 62.0

    print("")

    # Let's compute the angle of the diurnal path of a celestial body relative
    # to the horizon at the time of rising and setting
    dec = Angle(23.44)
    lat = Angle(40.0)
    j = diurnal_path_horizon(dec, lat)
    print_me("Diurnal path vs. horizon angle at time of rising and setting",
             j.dms_str(n_dec=1))                            # 45d 31' 28.4''


if __name__ == '__main__':

    main()
