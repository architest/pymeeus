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


from math import sin, cos, atan2, asin

from Angle import Angle
from Epoch import Epoch, JDE2000
from Coordinates import mean_obliquity


"""
.. module:: Sun
   :synopsis: Module including functions regarding Sun position
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


def sun_true_longitude_coarse(epoch):
    """This function provides the Sun's true longitude with a relatively low
    accuracy of about 0.01 degree.

    :param epoch: Epoch to compute the position of the Sun
    :type epoch: :py:class:`Epoch`

    :returns: A tuple containing the true (ecliptical) longitude (as an Angle
        object) and the radius vector in astronomical units.
    :rtype: tuple
    :raises: TypeError if input value is of wrong type.

    >>> epoch = Epoch(1992, 10, 13)
    >>> true_lon, r = sun_true_longitude_coarse(epoch)
    >>> print(true_lon.dms_str(n_dec=0))
    199d 54' 36.0''
    >>> print(round(r, 5))
    0.99766
    """

    # First check that input values are of correct types
    if not(isinstance(epoch, Epoch)):
        raise TypeError("Invalid input type")
    # Compute the time in Julian centuries
    t = (epoch - JDE2000)/36525.0
    # Compute the geometric mean longitude of the Sun
    l0 = 280.46646 + t*(36000.76983 + t*0.0003032)
    l0 = Angle(l0)
    l0.to_positive()
    # Now, compute the mean anomaly of the Sun
    m = 357.52911 + t*(35999.05029 - t*0.0001537)
    m = Angle(m)
    mrad = m.rad()
    # The eccentricity of the Earth's orbit
    e = 0.016708634 - t*(0.000042037 + t*0.0000001267)
    # Equation of the center
    c = (1.914602 - t*(0.004817 + t*0.000014))*sin(mrad) + \
        (0.019993 - t*0.000101)*sin(2.0*mrad) + \
        0.000289*sin(3.0*mrad)
    c = Angle(c)
    true_lon = l0 + c
    true_anom = m + c
    # Sun's radius vector
    r = (1.000001018*(1.0 - e*e))/(1.0 + e*cos(true_anom.rad()))
    return (true_lon, r)


def sun_apparent_longitude_coarse(epoch):
    """This function provides the Sun's apparent longitude with a relatively
    low accuracy of about 0.01 degree.

    :param epoch: Epoch to compute the position of the Sun
    :type epoch: :py:class:`Epoch`

    :returns: A tuple containing the sun_apparent (ecliptical) longitude (as an
        Angle object) and the radius vector in astronomical units.
    :rtype: tuple
    :raises: TypeError if input value is of wrong type.

    >>> epoch = Epoch(1992, 10, 13)
    >>> app_lon, r = sun_apparent_longitude_coarse(epoch)
    >>> print(app_lon.dms_str(n_dec=0))
    199d 54' 32.0''
    >>> print(round(r, 5))
    0.99766
    """

    # First find the true longitude
    true_lon, r = sun_true_longitude_coarse(epoch)
    # Compute the time in Julian centuries
    t = (epoch - JDE2000)/36525.0
    # Then correct for nutation and aberration
    omega = 125.04 - 1934.136*t
    omega = Angle(omega)
    lambd = true_lon - 0.00569 - 0.00478*sin(omega.rad())
    return (lambd, r)


def sun_apparent_rightascension_declination_coarse(epoch):
    """This function provides the Sun's apparent right ascension and
    declination with a relatively low accuracy of about 0.01 degree.

    :param epoch: Epoch to compute the position of the Sun
    :type epoch: :py:class:`Epoch`

    :returns: A tuple containing the right ascension and the declination (as
        Angle objects) and the radius vector in astronomical units.
    :rtype: tuple
    :raises: TypeError if input value is of wrong type.

    >>> epoch = Epoch(1992, 10, 13)
    >>> ra, delta, r = sun_apparent_rightascension_declination_coarse(epoch)
    >>> print(ra.ra_str(n_dec=1))
    13h 13' 31.4''
    >>> print(delta.dms_str(n_dec=0))
    -7d 47' 6.0''
    >>> print(round(r, 5))
    0.99766
    """

    # First find the apparent longitude
    app_lon, r = sun_apparent_longitude_coarse(epoch)
    # Compute the obliquity of the ecliptic
    e0 = mean_obliquity(epoch)
    # Compute the time in Julian centuries
    t = (epoch - JDE2000)/36525.0
    # Then correct for nutation and aberration
    omega = 125.04 - 1934.136*t
    omega = Angle(omega)
    # Correct the obliquity
    e = e0 + 0.00256*cos(omega.rad())
    alpha = atan2(cos(e.rad())*sin(app_lon.rad()), cos(app_lon.rad()))
    alpha = Angle(alpha, radians=True)
    alpha.to_positive()
    delta = asin(sin(e.rad())*sin(app_lon.rad()))
    delta = Angle(delta, radians=True)
    return (alpha, delta, r)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Sun functions
    print('\n' + 35*'*')
    print("*** Use of Sun functions")
    print(35*'*' + '\n')

    # Compute an approximation of the Sun's true longitude
    epoch = Epoch(1992, 10, 13)
    true_lon, r = sun_true_longitude_coarse(epoch)
    print_me("Sun's approximate true longitude", true_lon.dms_str(n_dec=0))
    # 199d 54' 36.0''
    print_me("Sun's radius vector", round(r, 5))                    # 0.99766

    print("")

    # Now let's compute the Sun's approximate apparent longitude
    app_lon, r = sun_apparent_longitude_coarse(epoch)
    print_me("Sun's approximate apparent longitude", app_lon.dms_str(n_dec=0))
    # 199d 54' 32.0''

    print("")

    # And now is the turn for the apparent right ascension and declination
    ra, delta, r = sun_apparent_rightascension_declination_coarse(epoch)
    print_me("Sun's apparent right ascension", ra.ra_str(n_dec=1))
    # 13h 13' 31.4''
    print_me("Sun's apparent declination", delta.dms_str(n_dec=0))
    # -7d 47' 6.0''


if __name__ == '__main__':

    main()
