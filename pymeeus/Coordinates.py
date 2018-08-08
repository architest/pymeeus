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


from math import sqrt, sin, cos, atan, atan2, asin
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
    :param utc: Whether the provided epoch is a civil time (UTC)
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
    """This method computes the true obliquity (epsilon) at the provided
    date. The true obliquity is the mean obliquity (epsilon0) plus the
    correction provided by the nutation in obliquity (Delta epsilon).

    This method internally uses an :class:`Epoch` object, and the
    **leap_seconds** argument then controls the way the UTC->TT conversion
    is handled for that object. If **leap_seconds** argument is set to a
    value different than zero, then that value will be used for the
    UTC->TAI conversion, and the internal leap seconds table will be
    bypassed. On the other hand, if it is set to zero, then the UTC to TT
    correction is disabled, and it is supposed that the input data is
    already in TT scale.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param leap_seconds: If different from zero, this is the value to be
       used in the UTC->TAI conversion. If equals to zero, conversion is
       disabled. If not given, UTC to TT conversion is carried out
       (default).
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
    """This method computes the nutation in longitude (Delta psi) at the
    provided date.

    This method internally uses an :class:`Epoch` object, and the
    **leap_seconds** argument then controls the way the UTC->TT conversion
    is handled for that object. If **leap_seconds** argument is set to a
    value different than zero, then that value will be used for the
    UTC->TAI conversion, and the internal leap seconds table will be
    bypassed. On the other hand, if it is set to zero, then the UTC to TT
    correction is disabled, and it is supposed that the input data is
    already in TT scale.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param leap_seconds: If different from zero, this is the value to be
       used in the UTC->TAI conversion. If equals to zero, conversion is
       disabled. If not given, UTC to TT conversion is carried out
       (default).
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
    """This method computes the nutation in obliquity (Delta epsilon) at
    the provided date.

    This method internally uses an :class:`Epoch` object, and the
    **leap_seconds** argument then controls the way the UTC->TT conversion
    is handled for that object. If **leap_seconds** argument is set to a
    value different than zero, then that value will be used for the
    UTC->TAI conversion, and the internal leap seconds table will be
    bypassed. On the other hand, if it is set to zero, then the UTC to TT
    correction is disabled, and it is supposed that the input data is
    already in TT scale.

    :param \*args: Either :class:`Epoch`, date, datetime or year, month,
        day values, by themselves or inside a tuple or list
    :type \*args: int, float, :py:class:`Epoch`, datetime, date, tuple,
        list
    :param leap_seconds: If different from zero, this is the value to be
       used in the UTC->TAI conversion. If equals to zero, conversion is
       disabled. If not given, UTC to TT conversion is carried out
       (default).
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
    """This method converts the equatorial coordinates (right ascension and
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
    """This method converts the ecliptical coordinates (longitude and
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
    not in ecliptical ones. Therefore, this method converts the provided
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
    """This method implements the Newcomb precessional equations used in
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
    """This method computes the star's true motion through space relative
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


if __name__ == '__main__':

    main()
