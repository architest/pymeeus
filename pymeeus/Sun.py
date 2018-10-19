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

        :param epoch: Epoch to compute Sun position, as an Epoch object
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
        199.907297
        >>> print(b.dms_str(n_dec=3))
        0.744''
        >>> print(round(r, 8))
        0.99760852
        """

        # NOTE: In page 169, Meeus gives a different value for the LONGITUDE
        # (199.907372 degrees) as the one presented above (199.907297 degrees).
        # After many checks and tests, I came to the conclusion that the result
        # above is the right one, and Meeus' result is wrong.
        # On the other hand, the difference in LATITUDE may be due to the fact
        # that this software uses the complete set of VSOP87C terms, instead of
        # the abridged version in Meeus' book.

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Use Earth heliocentric position to compute Sun's geocentric position
        lon, lat, r = Earth.geometric_heliocentric_position(epoch, toFK5)
        lon = lon.to_positive() + 180.0
        lat = -lat
        return lon, lat, r

    @staticmethod
    def apparent_geocentric_position(epoch):
        """"This method computes the apparent geocentric position of the Sun
        for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Sun position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> lon, lat, r = Sun.apparent_geocentric_position(epoch)
        >>> print(lon.to_positive().dms_str(n_dec=3))
        199d 54' 21.548''
        >>> print(lat.dms_str(n_dec=3))
        0.721''
        >>> print(round(r, 8))
        0.99760852
        """

        # NOTE: In page 169, Meeus gives a different value for the LONGITUDE
        # (199d 54' 21.818'') as the one presented above (199d 54' 21.548'').
        # After many checks and tests, I came to the conclusion that the result
        # above is the right one, and Meeus' result is wrong.
        # On the other hand, the difference in LATITUDE may be due to the fact
        # that this software uses the complete set of VSOP87C terms, instead of
        # the abridged version in Meeus' book.

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Use Earth heliocentric position to compute Sun's geocentric position
        lon, lat, r = Earth.apparent_heliocentric_position(epoch)
        lon = lon.to_positive() + 180.0
        lat = -lat
        return lon, lat, r

    @staticmethod
    def rectangular_coordinates_mean_equinox(epoch):
        """"This method computes the rectangular geocentric equatorial
        coordinates (X, Y, Z) of the Sun, referred to the mean equinox of the
        date. The X axis is directed towards the vernal equinox (longitude 0),
        the Y axis lies in the plane of the equator and is directed towards
        longitude 90, and the Z axis is directed towards the north celestial
        pole.

        :param epoch: Epoch to compute Sun position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the X, Y, Z values in astronomical units
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> x, y, z = Sun.rectangular_coordinates_mean_equinox(epoch)
        >>> print(round(x, 7))
        -0.9379963
        >>> print(round(y, 6))
        -0.311654
        >>> print(round(z, 7))
        -0.1351207
        """

        # NOTE: In page 172, Meeus gives slightly different values for x, y, z
        # as the ones internally computed in the example above. After many
        # checks and tests, I came to the conclusion that the results above are
        # the right ones, and Meeus' results are wrong.

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Get Sun's geocentric position with reduction to FK5 system
        lon, lat, r = Sun.geometric_geocentric_position(epoch)
        # Get the mean obliquity
        epsilon0 = mean_obliquity(epoch)
        # Change to radians
        ll = lon.rad()
        b = lat.rad()
        e = epsilon0.rad()
        # Compute the results
        x = r*cos(ll)
        y = r*(sin(ll)*cos(e) - sin(b)*sin(e))
        z = r*(sin(ll)*sin(e) + sin(b)*cos(e))
        return x, y, z

    @staticmethod
    def rectangular_coordinates_J2000(epoch):
        """"This method computes the rectangular geocentric equatorial
        coordinates (X, Y, Z) of the Sun, referred to the standard equinox of
        J2000.0. The X axis is directed towards the vernal equinox (longitude
        0), the Y axis lies in the plane of the equator and is directed towards
        longitude 90, and the Z axis is directed towards the north celestial
        pole.

        :param epoch: Epoch to compute Sun position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the X, Y, Z values in astronomical units
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> x, y, z = Sun.rectangular_coordinates_J2000(epoch)
        >>> print(round(x, 8))
        -0.93740485
        >>> print(round(y, 8))
        -0.3131474
        >>> print(round(z, 8))
        -0.12456646
        """

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Second, compute Earth heliocentric position referred to J2000.0
        lon, lat, r = Earth.geometric_heliocentric_position_J2000(epoch)
        # Third, convert from Earth's heliocentric to Sun's geocentric
        lon = lon.to_positive() + 180.0
        lat = -lat
        x = r*cos(lat.rad())*cos(lon.rad())
        y = r*cos(lat.rad())*sin(lon.rad())
        z = r*sin(lat.rad())
        x += 0.00000044036*y - 0.000000190919*z
        y = -0.000000479966*x + 0.917482137087*y - 0.397776982902*z
        z = 0.397776982902*y + 0.917482137087*z
        return x, y, z

    @staticmethod
    def rectangular_coordinates_B1950(epoch):
        """"This method computes the rectangular geocentric equatorial
        coordinates (X, Y, Z) of the Sun, referred to the mean equinox of
        B1950.0. The X axis is directed towards the vernal equinox (longitude
        0), the Y axis lies in the plane of the equator and is directed towards
        longitude 90, and the Z axis is directed towards the north celestial
        pole.

        :param epoch: Epoch to compute Sun position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the X, Y, Z values in astronomical units
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> x, y, z = Sun.rectangular_coordinates_B1950(epoch)
        >>> print(round(x, 8))
        -0.94149557
        >>> print(round(y, 8))
        -0.30259922
        >>> print(round(z, 8))
        -0.11578695
        """

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Second, compute Earth heliocentric position referred to J2000.0
        lon, lat, r = Earth.geometric_heliocentric_position_J2000(epoch)
        # Third, convert from Earth's heliocentric to Sun's geocentric
        lon = lon.to_positive() + 180.0
        lat = -lat
        x = r*cos(lat.rad())*cos(lon.rad())
        y = r*cos(lat.rad())*sin(lon.rad())
        z = r*sin(lat.rad())
        x = 0.999925702634*x + 0.012189716217*y + 0.000011134016*z
        y = -0.011179418036*x + 0.917413998946*y - 0.397777041885*z
        z = -0.004859003787*x + 0.397747363646*y + 0.917482111428*z
        return x, y, z

    @staticmethod
    def rectangular_coordinates_equinox(epoch, equinox_epoch):
        """"This method computes the rectangular geocentric equatorial
        coordinates (X, Y, Z) of the Sun, referred to an arbitrary mean
        equinox. The X axis is directed towards the vernal equinox (longitude
        0), the Y axis lies in the plane of the equator and is directed towards
        longitude 90, and the Z axis is directed towards the north celestial
        pole.

        :param epoch: Epoch to compute Sun position, as an Epoch object
        :type epoch: :py:class:`Epoch`
        :param equinox_epoch: Epoch corresponding to the mean equinox
        :type equinox_epoch: :py:class:`Epoch`

        :returns: A tuple with the X, Y, Z values in astronomical units
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> e_equinox = Epoch(2467616.0)
        >>> x, y, z = Sun.rectangular_coordinates_equinox(epoch, e_equinox)
        >>> print(round(x, 8))
        -0.93373777
        >>> print(round(y, 8))
        -0.32235109
        >>> print(round(z, 8))
        -0.12856709
        """

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Second, compute Sun's rectangular coordinates w.r.t. J2000.0
        x0, y0, z0 = Sun.rectangular_coordinates_J2000(epoch)
        # Third, computed auxiliary angles
        t = (equinox_epoch - JDE2000)/36525.0
        tt = (epoch - equinox_epoch)/36525.0
        # Compute the conversion parameters
        zeta = t*((2306.2181 + tt*(1.39656 - 0.000139*tt)) +
                  t*((0.30188 - 0.000344*tt) + 0.017998*t))
        z = t*((2306.2181 + tt*(1.39656 - 0.000139*tt)) +
               t*((1.09468 + 0.000066*tt) + 0.018203*t))
        theta = t*(2004.3109 + tt*(-0.85330 - 0.000217*tt) +
                   t*(-(0.42665 + 0.000217*tt) - 0.041833*t))
        # Redefine the former values as Angles, and compute them in radians
        zeta = Angle(0, 0, zeta)
        zetar = zeta.rad()
        z = Angle(0, 0, z)
        zr = z.rad()
        theta = Angle(0, 0, theta)
        thetar = theta.rad()
        xx = cos(zetar)*cos(zr)*cos(thetar) - sin(zetar)*sin(zr)
        xy = sin(zetar)*cos(zr) + cos(zetar)*sin(zr)*cos(thetar)
        xz = cos(zetar)*sin(thetar)
        yx = -cos(zetar)*sin(zr) - sin(zetar)*cos(zr)*cos(thetar)
        yy = cos(zetar)*cos(zr) - sin(zetar)*sin(zr)*cos(thetar)
        yz = -sin(zetar)*sin(thetar)
        zx = -cos(zr)*sin(thetar)
        zy = -sin(zr)*sin(thetar)
        zz = cos(thetar)
        xp = xx*x0 + yx*y0 + zx*z0
        yp = xy*x0 + yy*y0 + zy*z0
        zp = xz*x0 + yz*y0 + zz*z0
        return xp, yp, zp


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

    print("")

    # Let's compute Sun's postion, but more accurately
    epoch = Epoch(1992, 10, 13.0)
    l, b, r = Sun.geometric_geocentric_position(epoch, toFK5=False)
    print_me("Geometric Geocentric Longitude", round(l.to_positive(), 6))
    # 199.906016
    print_me("Geometric Geocentric Latitude", b.dms_str(n_dec=3))
    # 0.644''
    print_me("Radius vector", round(r, 8))
    # 0.99760775

    print("")

    # Compute Sun's apparent postion
    l, b, r = Sun.apparent_geocentric_position(epoch)
    print_me("Apparent Geocentric Longitude", l.to_positive().dms_str(n_dec=3))
    # 199d 54' 16.937''
    print_me("Apparent Geocentric Latitude", b.dms_str(n_dec=3))
    # 0.621''
    print_me("Radius vector", round(r, 8))
    # 0.99760775

    print("")

    # We can compute rectangular coordinates referred to mean equinox of date
    x, y, z = Sun.rectangular_coordinates_mean_equinox(epoch)
    print("Rectangular coordinates referred to mean equinox of date:")
    print_me("X", round(x, 7))                                  # -0.9379963
    print_me("Y", round(y, 6))                                  # -0.311654
    print_me("Z", round(z, 7))                                  # -0.1351207

    print("")

    # Now, compute rectangular coordinates w.r.t. standard equinox J2000.0
    x, y, z = Sun.rectangular_coordinates_J2000(epoch)
    print("Rectangular coordinates w.r.t. standard equinox J2000.0:")
    print_me("X", round(x, 8))                                  # -0.93740485
    print_me("Y", round(y, 8))                                  # -0.3131474
    print_me("Z", round(z, 8))                                  # -0.12456646

    print("")

    # Compute rectangular coordinates w.r.t. mean equinox of B1950.0
    x, y, z = Sun.rectangular_coordinates_B1950(epoch)
    print("Rectangular coordinates w.r.t. mean equinox of B1950.0:")
    print_me("X", round(x, 8))                                  # -0.94149557
    print_me("Y", round(y, 8))                                  # -0.30259922
    print_me("Z", round(z, 8))                                  # -0.11578695

    print("")

    # Compute rectangular coordinates w.r.t. an arbitrary mean equinox
    e_equinox = Epoch(2467616.0)
    x, y, z = Sun.rectangular_coordinates_equinox(epoch, e_equinox)
    print("Rectangular coordinates w.r.t. an arbitrary mean equinox:")
    print_me("X", round(x, 8))                                  # -0.93373777
    print_me("Y", round(y, 8))                                  # -0.32235109
    print_me("Z", round(z, 8))                                  # -0.12856709


if __name__ == '__main__':

    main()
