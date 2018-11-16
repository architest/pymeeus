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


from math import sin, cos, acos, atan2, sqrt

from Angle import Angle
from Epoch import Epoch
from Coordinates import kepler_equation
from Sun import Sun


"""
.. module:: Minor
   :synopsis: Class to model celestial bodies like comets and minor planets
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Minor(object):
    """
    Class Minor models minor celestial bodies.
    """

    @staticmethod
    def geocentric_position(a, e, i, omega, w, t, epoch):
        """This method computes the geocentric position of a minor celestial
        body (right ascension and declination) for the given epoch, and
        referred to the standard equinox J2000.0. Additionally, it also
        computes the elongation angle to the Sun.

        :param a: Semi-major axis of the orbit, in Astronomical Units
        :type a: float
        :param e: Eccentricity of the orbit
        :type e: float
        :param i: Inclination of the orbit, as an Angle object
        :type i: :py:class:`Angle`
        :param omega: Longitude of the ascending node, as an Angle object
        :type omega: :py:class:`Angle`
        :param w: Argument of the perihelion, as an Angle object
        :type w: :py:class:`Angle`
        :param t: Epoch of passage by perihelion, as an Epoch object
        :type t: :py:class:`Epoch`
        :param epoch: Epoch to compute geocentric position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the right ascension, the declination and
            the elongation angle to the Sun, as Angle objects
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> a = 2.2091404
        >>> e = 0.8502196
        >>> i = Angle(11.94524)
        >>> omega = Angle(334.75006)
        >>> w = Angle(186.23352)
        >>> t = Epoch(1990, 10, 28.54502)
        >>> epoch = Epoch(1990, 10, 6.0)
        >>> ra, dec, p = Minor.geocentric_position(a, e, i, omega, w, t, epoch)
        >>> print(ra.ra_str(n_dec=1))
        10h 34' 13.7''
        >>> print(dec.dms_str(n_dec=0))
        19d 9' 32.0''
        >>> print(round(p, 2))
        40.51
        """

        # First check that input value is of correct types
        if not (isinstance(epoch, Epoch) and isinstance(a, float) and
                isinstance(e, float) and isinstance(i, Angle) and
                isinstance(omega, Angle) and isinstance(w, Angle) and
                isinstance(t, Epoch)):
            raise TypeError("Invalid input types")
        # Compute auxiliary quantities
        se = 0.397777156
        ce = 0.917482062
        omer = omega.rad()
        ir = i.rad()
        f = cos(omer)
        g = sin(omer) * ce
        h = sin(omer) * se
        p = -sin(omer) * cos(ir)
        q = cos(omer) * cos(ir) * ce - sin(ir) * se
        r = cos(omer) * cos(ir) * se + sin(ir) * ce
        aa = atan2(f, p)
        bb = atan2(g, q)
        cc = atan2(h, r)
        am = sqrt(f * f + p * p)
        bm = sqrt(g * g + q * q)
        cm = sqrt(h * h + r * r)
        # Compute the mean motion from the semi-major axis (degrees/day)
        n = 0.9856076686/(a * sqrt(a))
        # Time since perihelion
        t_peri = epoch - t
        # Now, compute the mean anomaly, in degrees
        m = t_peri * n
        m = Angle(m)
        # With the mean anomaly, use Kepler's equation to find E and v
        ee, v = kepler_equation(e, m)
        ee = Angle(ee).to_positive()
        # Get r
        er = ee.rad()
        rr = a * (1.0 - e * cos(er))
        # Compute the heliocentric rectangular equatorial coordinates
        wr = w.rad()
        vr = Angle(v).rad()
        x = rr * am * sin(aa + wr + vr)
        y = rr * bm * sin(bb + wr + vr)
        z = rr * cm * sin(cc + wr + vr)
        # Now let's compute Sun's rectangular coordinates
        xs, ys, zs = Sun.rectangular_coordinates_J2000(epoch)
        xi = x + xs
        eta = y + ys
        zeta = z + zs
        delta = sqrt(xi * xi + eta * eta + zeta * zeta)
        # We need to correct for the effect of light-time. Compute delay tau
        tau = 0.0057755183 * delta
        # Recompute some critical parameters
        t_peri = epoch - t - tau
        # Now, compute the mean anomaly, in degrees
        m = t_peri * n
        m = Angle(m)
        # With the mean anomaly, use Kepler's equation to find E and v
        ee, v = kepler_equation(e, m)
        ee = Angle(ee).to_positive()
        # Get r
        er = ee.rad()
        rr = a * (1.0 - e * cos(er))
        # Compute the heliocentric rectangular equatorial coordinates
        wr = w.rad()
        vr = Angle(v).rad()
        x = rr * am * sin(aa + wr + vr)
        y = rr * bm * sin(bb + wr + vr)
        z = rr * cm * sin(cc + wr + vr)
        xi = x + xs
        eta = y + ys
        zeta = z + zs
        ra = Angle(atan2(eta, xi), radians=True)
        dec = Angle(atan2(zeta, sqrt(xi * xi + eta * eta)), radians=True)
        r_sun = sqrt(xs * xs + ys * ys + zs * zs)
        psi = acos((xi * xs + eta * ys + zeta * zs) / (r_sun * delta))
        psi = Angle(psi, radians=True)
        return ra, dec, psi

    @staticmethod
    def heliocentric_ecliptical_position(a, e, i, omega, w, t, epoch):
        """This method computes the heliocentric position of a minor celestial
        body, providing the result in ecliptical coordinates.

        :param a: Semi-major axis of the orbit, in Astronomical Units
        :type a: float
        :param e: Eccentricity of the orbit
        :type e: float
        :param i: Inclination of the orbit, as an Angle object
        :type i: :py:class:`Angle`
        :param omega: Longitude of the ascending node, as an Angle object
        :type omega: :py:class:`Angle`
        :param w: Argument of the perihelion, as an Angle object
        :type w: :py:class:`Angle`
        :param t: Epoch of passage by perihelion, as an Epoch object
        :type t: :py:class:`Epoch`
        :param epoch: Epoch to compute geocentric position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing longitude and latitude, as Angle objects
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> a = 2.2091404
        >>> e = 0.8502196
        >>> i = Angle(11.94524)
        >>> omega = Angle(334.75006)
        >>> w = Angle(186.23352)
        >>> t = Epoch(1990, 10, 28.54502)
        >>> epoch = Epoch(1990, 10, 6.0)
        >>> lon, lat = Minor.heliocentric_ecliptical_position(a, e, i, omega, \
                                                              w, t, epoch)
        >>> print(lon.dms_str(n_dec=1))
        66d 51' 57.8''
        >>> print(lat.dms_str(n_dec=1))
        11d 56' 14.4''
        """

        # First check that input value is of correct types
        if not (isinstance(epoch, Epoch) and isinstance(a, float) and
                isinstance(e, float) and isinstance(i, Angle) and
                isinstance(omega, Angle) and isinstance(w, Angle) and
                isinstance(t, Epoch)):
            raise TypeError("Invalid input types")
        # Compute the mean motion from the semi-major axis (degrees/day)
        n = 0.9856076686/(a * sqrt(a))
        # Time since perihelion
        t_peri = epoch - t
        # Now, compute the mean anomaly, in degrees
        m = t_peri * n
        m = Angle(m)
        # With the mean anomaly, use Kepler's equation to find E and v
        ee, v = kepler_equation(e, m)
        ee = Angle(ee).to_positive()
        # Get r
        er = ee.rad()
        r = a * (1.0 - e * cos(er))
        # Compute the heliocentric rectangular ecliptical coordinates
        wr = w.rad()
        vr = Angle(v).rad()
        ur = wr + vr
        omer = omega.rad()
        ir = i.rad()
        x = r * (cos(omer) * cos(ur) - sin(omer) * sin(ur) * cos(ir))
        y = r * (sin(omer) * cos(ur) + cos(omer) * sin(ur) * cos(ir))
        z = r * sin(ir) * sin(ur)
        lon = atan2(y, x)
        lat = atan2(z, sqrt(x * x + y * y))
        return Angle(lon, radians=True), Angle(lat, radians=True)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Minor class
    print("\n" + 35 * "*")
    print("*** Use of Minor class")
    print(35 * "*" + "\n")

    # Let's compute the equatorial coordinates of comet Encke
    a = 2.2091404
    e = 0.8502196
    i = Angle(11.94524)
    omega = Angle(334.75006)
    w = Angle(186.23352)
    t = Epoch(1990, 10, 28.54502)
    epoch = Epoch(1990, 10, 6.0)
    ra, dec, elong = Minor.geocentric_position(a, e, i, omega, w, t, epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))     # 10h 34' 13.7''
    print_me("Declination", dec.dms_str(n_dec=0))       # 19d 9' 32.0''
    print_me("Elongation", round(elong, 2))             # 40.51

    print("")

    # Now compute the heliocentric ecliptical coordinates
    a = 2.2091404
    e = 0.8502196
    i = Angle(11.94524)
    omega = Angle(334.75006)
    w = Angle(186.23352)
    t = Epoch(1990, 10, 28.54502)
    epoch = Epoch(1990, 10, 6.0)
    lon, lat = Minor.heliocentric_ecliptical_position(a, e, i, omega,
                                                      w, t, epoch)
    print_me("Heliocentric ecliptical longitude", lon.dms_str(n_dec=1))
    # 66d 51' 57.8''
    print_me("Heliocentric ecliptical latitude", lat.dms_str(n_dec=1))
    # 11d 56' 14.4''


if __name__ == "__main__":

    main()
