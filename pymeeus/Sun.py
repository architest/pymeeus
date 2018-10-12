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
from Earth import Earth


"""
.. module:: Sun
   :synopsis: Module including functions regarding Sun position
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Sun(object):
    """
    Class Sun handles the parameters related to the Sun.
    """

    def __init__(self):
        """Sun constructor.

        :returns: Sun object.
        :rtype: :py:class:`Sun`
        """

    @staticmethod
    def true_longitude_coarse(epoch):
        """This method provides the Sun's true longitude with a relatively low
        accuracy of about 0.01 degree.

        :param epoch: Epoch to compute the position of the Sun
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the true (ecliptical) longitude (as an
            Angle object) and the radius vector in astronomical units.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 10, 13)
        >>> true_lon, r = Sun.true_longitude_coarse(epoch)
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

    @staticmethod
    def apparent_longitude_coarse(epoch):
        """This method provides the Sun's apparent longitude with a relatively
        low accuracy of about 0.01 degree.

        :param epoch: Epoch to compute the position of the Sun
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the sun_apparent (ecliptical) longitude
            (as an Angle object) and the radius vector in astronomical units.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 10, 13)
        >>> app_lon, r = Sun.apparent_longitude_coarse(epoch)
        >>> print(app_lon.dms_str(n_dec=0))
        199d 54' 32.0''
        >>> print(round(r, 5))
        0.99766
        """

        # First find the true longitude
        sun = Sun()
        true_lon, r = sun.true_longitude_coarse(epoch)
        # Compute the time in Julian centuries
        t = (epoch - JDE2000)/36525.0
        # Then correct for nutation and aberration
        omega = 125.04 - 1934.136*t
        omega = Angle(omega)
        lambd = true_lon - 0.00569 - 0.00478*sin(omega.rad())
        return (lambd, r)

    @staticmethod
    def apparent_rightascension_declination_coarse(epoch):
        """This method provides the Sun's apparent right ascension and
        declination with a relatively low accuracy of about 0.01 degree.

        :param epoch: Epoch to compute the position of the Sun
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the right ascension and the declination
            (as Angle objects) and the radius vector in astronomical units.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epo = Epoch(1992, 10, 13)
        >>> ra, delta, r = Sun.apparent_rightascension_declination_coarse(epo)
        >>> print(ra.ra_str(n_dec=1))
        13h 13' 31.4''
        >>> print(delta.dms_str(n_dec=0))
        -7d 47' 6.0''
        >>> print(round(r, 5))
        0.99766
        """

        # First find the apparent longitude
        sun = Sun()
        app_lon, r = sun.apparent_longitude_coarse(epoch)
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

    @staticmethod
    def geometric_geocentric_position(epoch, toFK5=True):
        """"This method computes the geometric geocentric position of the Sun
        for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Venus position, as an Epoch object
        :type epoch: :py:class:`Epoch`
        :param toFK5: Whether or not the small correction to convert to the FK5
            system will be applied or not
        :type toFK5: bool

        :returns: A tuple with the geocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> l, b, r = Sun.geometric_geocentric_position(epoch, toFK5=False)
        >>> print(round(l.to_positive(), 6))
        199.906016
        >>> print(b.dms_str(n_dec=3))
        0.644''
        >>> print(round(r, 8))
        0.99760775
        """

        # NOTE: In page 169, Meeus gives a different value for the LONGITUDE
        # (199.907372 degrees) as the one presented above (199.906016 degrees).
        # After many checks and tests, I came to the conclusion that the result
        # above is the right one, and Meeus' result is wrong

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Use Earth heliocentric position to compute Sun's geocentric position
        lon, lat, r = Earth.geometric_heliocentric_position(epoch, toFK5)
        lon = lon.to_positive() + 180.0
        lat = -lat
        return lon, lat, r


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Sun functions
    print('\n' + 35*'*')
    print("*** Use of Sun class")
    print(35*'*' + '\n')

    # Compute an approximation of the Sun's true longitude
    epoch = Epoch(1992, 10, 13)
    true_lon, r = Sun.true_longitude_coarse(epoch)
    print_me("Sun's approximate true longitude", true_lon.dms_str(n_dec=0))
    # 199d 54' 36.0''
    print_me("Sun's radius vector", round(r, 5))                    # 0.99766

    print("")

    # Now let's compute the Sun's approximate apparent longitude
    app_lon, r = Sun.apparent_longitude_coarse(epoch)
    print_me("Sun's approximate apparent longitude", app_lon.dms_str(n_dec=0))
    # 199d 54' 32.0''

    print("")

    # And now is the turn for the apparent right ascension and declination
    ra, delta, r = Sun.apparent_rightascension_declination_coarse(epoch)
    print_me("Sun's apparent right ascension", ra.ra_str(n_dec=1))
    # 13h 13' 31.4''
    print_me("Sun's apparent declination", delta.dms_str(n_dec=0))
    # -7d 47' 6.0''


if __name__ == '__main__':

    main()
