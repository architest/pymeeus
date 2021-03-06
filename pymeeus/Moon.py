
# -*- coding: utf-8 -*-


# PyMeeus: Python module implementing astronomical algorithms.
# Copyright (C) 2021  Dagoberto Salazar
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


from math import sin, cos, asin, atan2
from pymeeus.base import iint
from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch, JDE2000
from pymeeus.Sun import Sun
from pymeeus.Coordinates import (
    nutation_longitude, true_obliquity, ecliptical2equatorial
)

"""
.. module:: Moon
   :synopsis: Class to model the Moon
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


PERIODIC_TERMS_LR_TABLE = [
    [0,  0,  1,  0, 6288774.0, -20905355.0],
    [2,  0, -1,  0, 1274027.0,  -3699111.0],
    [2,  0,  0,  0,  658314.0,  -2955968.0],
    [0,  0,  2,  0,  213618.0,   -569925.0],
    [0,  1,  0,  0, -185116.0,     48888.0],
    [0,  0,  0,  2, -114332.0,     -3149.0],
    [2,  0, -2,  0,   58793.0,    246158.0],
    [2, -1, -1,  0,   57066.0,   -152138.0],
    [2,  0,  1,  0,   53322.0,   -170733.0],
    [2, -1,  0,  0,   45758.0,   -204586.0],
    [0,  1, -1,  0,  -40923.0,   -129620.0],
    [1,  0,  0,  0,  -34720.0,    108743.0],
    [0,  1,  1,  0,  -30383.0,    104755.0],
    [2,  0,  0, -2,   15327.0,     10321.0],
    [0,  0,  1,  2,  -12528.0,         0.0],
    [0,  0,  1, -2,   10980.0,     79661.0],
    [4,  0, -1,  0,   10675.0,    -34782.0],
    [0,  0,  3,  0,   10034.0,    -23210.0],
    [4,  0, -2,  0,    8548.0,    -21636.0],
    [2,  1, -1,  0,   -7888.0,     24208.0],
    [2,  1,  0,  0,   -6766.0,     30824.0],
    [1,  0, -1,  0,   -5163.0,     -8379.0],
    [1,  1,  0,  0,    4987.0,    -16675.0],
    [2, -1,  1,  0,    4036.0,    -12831.0],
    [2,  0,  2,  0,    3994.0,    -10445.0],
    [4,  0,  0,  0,    3861.0,    -11650.0],
    [2,  0, -3,  0,    3665.0,     14403.0],
    [0,  1, -2,  0,   -2689.0,     -7003.0],
    [2,  0, -1,  2,   -2602.0,         0.0],
    [2, -1, -2,  0,    2390.0,     10056.0],
    [1,  0,  1,  0,   -2348.0,      6322.0],
    [2, -2,  0,  0,    2236.0,     -9884.0],
    [0,  1,  2,  0,   -2120.0,      5751.0],
    [0,  2,  0,  0,   -2069.0,         0.0],
    [2, -2, -1,  0,    2048.0,     -4950.0],
    [2,  0,  1, -2,   -1773.0,      4130.0],
    [2,  0,  0,  2,   -1595.0,         0.0],
    [4, -1, -1,  0,    1215.0,     -3958.0],
    [0,  0,  2,  2,   -1110.0,         0.0],
    [3,  0, -1,  0,    -892.0,      3258.0],
    [2,  1,  1,  0,    -810.0,      2616.0],
    [4, -1, -2,  0,     759.0,     -1897.0],
    [0,  2, -1,  0,    -713.0,     -2117.0],
    [2,  2, -1,  0,    -700.0,      2354.0],
    [2,  1, -2,  0,     691.0,         0.0],
    [2, -1,  0, -2,     596.0,         0.0],
    [4,  0,  1,  0,     549.0,     -1423.0],
    [0,  0,  4,  0,     537.0,     -1117.0],
    [4, -1,  0,  0,     520.0,     -1571.0],
    [1,  0, -2,  0,    -487.0,     -1739.0],
    [2,  1,  0, -2,    -399.0,         0.0],
    [0,  0,  2, -2,    -381.0,     -4421.0],
    [1,  1,  1,  0,     351.0,         0.0],
    [3,  0, -2,  0,    -340.0,         0.0],
    [4,  0, -3,  0,     330.0,         0.0],
    [2, -1,  2,  0,     327.0,         0.0],
    [0,  2,  1,  0,    -323.0,      1165.0],
    [1,  1, -1,  0,     299.0,         0.0],
    [2,  0,  3,  0,     294.0,         0.0],
    [2,  0, -1, -2,       0.0,      8752.0]
]
"""This table contains the periodic terms for the longitude (Sigmal) and
distance (Sigmar) of the Moon. Units are 0.000001 degree for Sigmal, and 0.001
kilometer for Sigmar. In Meeus' book this is Table 47.A and can be found in
pages 339-340."""

PERIODIC_TERMS_B_TABLE = [
    [0,  0,  0,  1, 5128122.0],
    [0,  0,  1,  1,  280602.0],
    [0,  0,  1, -1,  277693.0],
    [2,  0,  0, -1,  173237.0],
    [2,  0, -1,  1,   55413.0],
    [2,  0, -1, -1,   46271.0],
    [2,  0,  0,  1,   32573.0],
    [0,  0,  2,  1,   17198.0],
    [2,  0,  1, -1,    9266.0],
    [0,  0,  2, -1,    8822.0],
    [2, -1,  0, -1,    8216.0],
    [2,  0, -2, -1,    4324.0],
    [2,  0,  1,  1,    4200.0],
    [2,  1,  0, -1,   -3359.0],
    [2, -1, -1,  1,    2463.0],
    [2, -1,  0,  1,    2211.0],
    [2, -1, -1, -1,    2065.0],
    [0,  1, -1, -1,   -1870.0],
    [4,  0, -1, -1,    1828.0],
    [0,  1,  0,  1,   -1794.0],
    [0,  0,  0,  3,   -1749.0],
    [0,  1, -1,  1,   -1565.0],
    [1,  0,  0,  1,   -1491.0],
    [0,  1,  1,  1,   -1475.0],
    [0,  1,  1, -1,   -1410.0],
    [0,  1,  0, -1,   -1344.0],
    [1,  0,  0, -1,   -1335.0],
    [0,  0,  3,  1,    1107.0],
    [4,  0,  0, -1,    1021.0],
    [4,  0, -1,  1,     833.0],
    [0,  0,  1, -3,     777.0],
    [4,  0, -2,  1,     671.0],
    [2,  0,  0, -3,     607.0],
    [2,  0,  2, -1,     596.0],
    [2, -1,  1, -1,     491.0],
    [2,  0, -2,  1,    -451.0],
    [0,  0,  3, -1,     439.0],
    [2,  0,  2,  1,     422.0],
    [2,  0, -3, -1,     421.0],
    [2,  1, -1,  1,    -366.0],
    [2,  1,  0,  1,    -351.0],
    [4,  0,  0,  1,     331.0],
    [2, -1,  1,  1,     315.0],
    [2, -2,  0, -1,     302.0],
    [0,  0,  1,  3,    -283.0],
    [2,  1,  1, -1,    -229.0],
    [1,  1,  0, -1,     223.0],
    [1,  1,  0,  1,     223.0],
    [0,  1, -2, -1,    -220.0],
    [2,  1, -1, -1,    -220.0],
    [1,  0,  1,  1,    -185.0],
    [2, -1, -2, -1,     181.0],
    [0,  1,  2,  1,    -177.0],
    [4,  0, -2, -1,     176.0],
    [4, -1, -1, -1,     166.0],
    [1,  0,  1, -1,    -164.0],
    [4,  0,  1, -1,     132.0],
    [1,  0, -1, -1,    -119.0],
    [4, -1,  0, -1,     115.0],
    [2, -2,  0,  1,     107.0],
]
"""This table contains the periodic terms for the latitude of the Moon Sigmab.
Units are 0.000001 degree. In Meeus' book this is Table 47.B and can be found
in page 341."""


class Moon(object):
    """
    Class Moon models Earth's satellite.
    """

    @staticmethod
    def geocentric_ecliptical_pos(epoch):
        """This method computes the geocentric ecliptical position (longitude,
        latitude) of the Moon for a given instant, referred to the mean equinox
        of the date, as well as the Moon-Earth distance in kilometers and the
        equatorial horizontal parallax.

        :param epoch: Instant to compute the Moon's position, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: Tuple containing:

            * Geocentric longitude of the center of the Moon, as an
              py:class:`Epoch` object.
            * Geocentric latitude of the center of the Moon, as an
              py:class:`Epoch` object.
            * Distance in kilometers between the centers of Earth and Moon, in
              kilometers (float)
            * Equatorial horizontal parallax of the Moon, as an
              py:class:`Epoch` object.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 4, 12.0)
        >>> Lambda, Beta, Delta, ppi = Moon.geocentric_ecliptical_pos(epoch)
        >>> print(round(Lambda, 6))
        133.162655
        >>> print(round(Beta, 6))
        -3.229126
        >>> print(round(Delta, 1))
        368409.7
        >>> print(round(ppi, 5))
        0.99199
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Get the time from J2000.0 in Julian centuries
        t = (epoch - JDE2000) / 36525.0
        # Compute Moon's mean longitude, referred to mean equinox of date
        Lprime = 218.3164477 + (481267.88123421
                                + (-0.0015786
                                   + (1.0/538841.0
                                       - t/65194000.0) * t) * t) * t
        # Mean elongation of the Moon
        D = 297.8501921 + (445267.1114034
                           + (-0.0018819
                              + (1.0/545868.0 - t/113065000.0) * t) * t) * t
        # Sun's mean anomaly
        M = 357.5291092 + (35999.0502909 + (-0.0001536 + t/24490000.0) * t) * t
        # Moon's mean anomaly
        Mprime = 134.9633964 + (477198.8675055
                                + (0.0087414
                                   + (1.0/69699.9
                                      + t/14712000.0) * t) * t) * t
        # Moon's argument of latitude
        F = 93.2720950 + (483202.0175233
                          + (-0.0036539
                             + (-1.0/3526000.0 + t/863310000.0) * t) * t) * t
        # Let's compute some additional arguments
        A1 = 119.75 + 131.849 * t
        A2 = 53.09 + 479264.290 * t
        A3 = 313.45 + 481266.484 * t
        # Eccentricity of Earth's orbit around the Sun
        E = 1.0 + (-0.002516 - 0.0000074 * t) * t
        E2 = E * E
        # Reduce the angles to a [0 360] range
        Lprime = Angle(Angle.reduce_deg(Lprime)).to_positive()
        Lprimer = Lprime.rad()
        D = Angle(Angle.reduce_deg(D)).to_positive()
        Dr = D.rad()
        M = Angle(Angle.reduce_deg(M)).to_positive()
        Mr = M.rad()
        Mprime = Angle(Angle.reduce_deg(Mprime)).to_positive()
        Mprimer = Mprime.rad()
        F = Angle(Angle.reduce_deg(F)).to_positive()
        Fr = F.rad()
        A1 = Angle(Angle.reduce_deg(A1)).to_positive()
        A1r = A1.rad()
        A2 = Angle(Angle.reduce_deg(A2)).to_positive()
        A2r = A2.rad()
        A3 = Angle(Angle.reduce_deg(A3)).to_positive()
        A3r = A3.rad()
        # Let's store this results in a list, in preparation for using tables
        arguments = [Dr, Mr, Mprimer, Fr]
        # Now we use the tables of periodic terms. First for sigmal and sigmar
        sigmal = 0.0
        sigmar = 0.0
        for i, value in enumerate(PERIODIC_TERMS_LR_TABLE):
            argument = 0.0
            for j in range(4):
                if PERIODIC_TERMS_LR_TABLE[i][j]:  # Avoid multiply by zero
                    argument += PERIODIC_TERMS_LR_TABLE[i][j] * arguments[j]
            coeffl = value[4]
            coeffr = value[5]
            if abs(value[1]) == 1:
                coeffl = coeffl * E
                coeffr = coeffr * E
            elif abs(value[1]) == 2:
                coeffl = coeffl * E2
                coeffr = coeffr * E2
            sigmal += coeffl * sin(argument)
            sigmar += coeffr * cos(argument)
        # Add the additive terms to sigmal
        sigmal += (3958.0 * sin(A1r) + 1962.0 * sin(Lprimer - Fr)
                   + 318.0 * sin(A2r))
        # Now use the tabla for sigmab
        sigmab = 0.0
        for i, value in enumerate(PERIODIC_TERMS_B_TABLE):
            argument = 0.0
            for j in range(4):
                if PERIODIC_TERMS_B_TABLE[i][j]:  # Avoid multiply by zero
                    argument += PERIODIC_TERMS_B_TABLE[i][j] * arguments[j]
            coeffb = value[4]
            if abs(value[1]) == 1:
                coeffb = coeffb * E
            elif abs(value[1]) == 2:
                coeffb = coeffb * E2
            sigmab += coeffb * sin(argument)
        # Add the additive terms to sigmab
        sigmab += (-2235.0 * sin(Lprimer) + 382.0 * sin(A3r)
                   + 175.0 * sin(A1r - Fr) + 175.0 * sin(A1r + Fr)
                   + 127.0 * sin(Lprimer - Mprimer)
                   - 115.0 * sin(Lprimer + Mprimer))
        Lambda = Lprime + (sigmal / 1000000.0)
        Beta = Angle(sigmab / 1000000.0)
        Delta = 385000.56 + (sigmar / 1000.0)
        ppii = asin(6378.14 / Delta)
        ppi = Angle(ppii, radians=True)
        return Lambda, Beta, Delta, ppi

    @staticmethod
    def apparent_ecliptical_pos(epoch):
        """This method computes the apparent geocentric ecliptical position
        (longitude, latitude) of the Moon for a given instant, referred to the
        mean equinox of the date, as well as the Moon-Earth distance in
        kilometers and the equatorial horizontal parallax.

        :param epoch: Instant to compute the Moon's position, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: Tuple containing:

            * Apparent geocentric longitude of the center of the Moon, as an
              py:class:`Epoch` object.
            * Apparent geocentric latitude of the center of the Moon, as an
              py:class:`Epoch` object.
            * Distance in kilometers between the centers of Earth and Moon, in
              kilometers (float)
            * Equatorial horizontal parallax of the Moon, as an
              py:class:`Epoch` object.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 4, 12.0)
        >>> Lambda, Beta, Delta, ppi = Moon.apparent_ecliptical_pos(epoch)
        >>> print(round(Lambda, 5))
        133.16726
        >>> print(round(Beta, 6))
        -3.229126
        >>> print(round(Delta, 1))
        368409.7
        >>> print(round(ppi, 5))
        0.99199
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Now, let's call the method geocentric_ecliptical_pos()
        Lambda, Beta, Delta, ppi = Moon.geocentric_ecliptical_pos(epoch)
        # Compute the nutation in longitude (deltaPsi)
        deltaPsi = nutation_longitude(epoch)
        # Correct the longitude to obtain the apparent longitude
        aLambda = Lambda + deltaPsi
        return aLambda, Beta, Delta, ppi

    @staticmethod
    def apparent_equatorial_pos(epoch):
        """This method computes the apparent equatorial position (right
        ascension, declination) of the Moon for a given instant, referred to
        the mean equinox of the date, as well as the Moon-Earth distance in
        kilometers and the equatorial horizontal parallax.

        :param epoch: Instant to compute the Moon's position, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: Tuple containing:

            * Apparent right ascension of the center of the Moon, as an
              py:class:`Epoch` object.
            * Apparent declination of the center of the Moon, as an
              py:class:`Epoch` object.
            * Distance in kilometers between the centers of Earth and Moon, in
              kilometers (float)
            * Equatorial horizontal parallax of the Moon, as an
              py:class:`Epoch` object.
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 4, 12.0)
        >>> ra, dec, Delta, ppi = Moon.apparent_equatorial_pos(epoch)
        >>> print(round(ra, 6))
        134.688469
        >>> print(round(dec, 6))
        13.768367
        >>> print(round(Delta, 1))
        368409.7
        >>> print(round(ppi, 5))
        0.99199
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Let's start calling the method 'apparent_ecliptical_pos()'
        Lambda, Beta, Delta, ppi = Moon.apparent_ecliptical_pos(epoch)
        # Now we need the obliquity of the ecliptic
        epsilon = true_obliquity(epoch)
        # And now let's carry out the transformation ecliptical->equatorial
        ra, dec = ecliptical2equatorial(Lambda, Beta, epsilon)
        return ra, dec, Delta, ppi

    @staticmethod
    def longitude_mean_ascending_node(epoch):
        """This method computes the longitude of the mean ascending node of the
        Moon in degrees, for a given instant, measured from the mean equinox of
        the date.

        :param epoch: Instant to compute the Moon's mean ascending node, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: The longitude of the mean ascending node.
        :rtype: py:class:`Angle`
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1913, 5, 27.0)
        >>> Omega = Moon.longitude_mean_ascending_node(epoch)
        >>> print(round(Omega, 1))
        0.0
        >>> epoch = Epoch(2043, 9, 10.0)
        >>> Omega = Moon.longitude_mean_ascending_node(epoch)
        >>> print(round(Omega, 1))
        0.0
        >>> epoch = Epoch(1959, 12, 7.0)
        >>> Omega = Moon.longitude_mean_ascending_node(epoch)
        >>> print(round(Omega, 1))
        180.0
        >>> epoch = Epoch(2108, 11, 3.0)
        >>> Omega = Moon.longitude_mean_ascending_node(epoch)
        >>> print(round(Omega, 1))
        180.0
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Get the time from J2000.0 in Julian centuries
        t = (epoch - JDE2000) / 36525.0
        # Compute Moon's longitude of the mean ascending node
        Omega = 125.0445479 + (-1934.1362891
                               + (0.0020754
                                  + (1.0/476441.0
                                     - t/60616000.0) * t) * t) * t
        Omega = Angle(Omega).to_positive()
        return Omega

    @staticmethod
    def longitude_true_ascending_node(epoch):
        """This method computes the longitude of the true ascending node of the
        Moon in degrees, for a given instant, measured from the mean equinox of
        the date.

        :param epoch: Instant to compute the Moon's true ascending node, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: The longitude of the true ascending node.
        :rtype: py:class:`Angle`
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1913, 5, 27.0)
        >>> Omega = Moon.longitude_true_ascending_node(epoch)
        >>> print(round(Omega, 4))
        0.8763
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Let's start computing the longitude of the MEAN ascending node
        Omega = Moon.longitude_mean_ascending_node(epoch)
        # Get the time from J2000.0 in Julian centuries
        t = (epoch - JDE2000) / 36525.0
        # Mean elongation of the Moon
        D = 297.8501921 + (445267.1114034
                           + (-0.0018819
                              + (1.0/545868.0 - t/113065000.0) * t) * t) * t
        # Sun's mean anomaly
        M = 357.5291092 + (35999.0502909 + (-0.0001536 + t/24490000.0) * t) * t
        # Moon's mean anomaly
        Mprime = 134.9633964 + (477198.8675055
                                + (0.0087414
                                   + (1.0/69699.9
                                      + t/14712000.0) * t) * t) * t
        # Moon's argument of latitude
        F = 93.2720950 + (483202.0175233
                          + (-0.0036539
                             + (-1.0/3526000.0 + t/863310000.0) * t) * t) * t
        # Reduce the angles to a [0 360] range
        D = Angle(Angle.reduce_deg(D)).to_positive()
        Dr = D.rad()
        M = Angle(Angle.reduce_deg(M)).to_positive()
        Mr = M.rad()
        Mprime = Angle(Angle.reduce_deg(Mprime)).to_positive()
        Mprimer = Mprime.rad()
        F = Angle(Angle.reduce_deg(F)).to_positive()
        Fr = F.rad()
        # Compute the periodic terms
        corr = (-1.4979 * sin(2.0 * (Dr - Fr)) - 0.15 * sin(Mr)
                - 0.1226 * sin(2.0 * Dr) + 0.1176 * sin(2.0 * Fr)
                - 0.0801 * sin(2.0 * (Mprimer - Fr)))
        Omega += Angle(corr)
        return Omega

    @staticmethod
    def longitude_mean_perigee(epoch):
        """This method computes the longitude of the mean perigee of the lunar
        orbitn in degrees, for a given instant, measured from the mean equinox
        of the date.

        :param epoch: Instant to compute the Moon's mean perigee, as an
            py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: The longitude of the mean perigee.
        :rtype: py:class:`Angle`
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(2021, 3, 5.0)
        >>> Pi = Moon.longitude_mean_perigee(epoch)
        >>> print(round(Pi, 5))
        224.89194
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Get the time from J2000.0 in Julian centuries
        t = (epoch - JDE2000) / 36525.0
        # Compute Moon's longitude of the mean perigee
        ppii = 83.3532465 + (4069.0137287
                             + (-0.01032
                                + (-1.0/80053.0 + t/18999000.0) * t) * t) * t
        ppii = Angle(ppii)
        return ppii

    @staticmethod
    def illuminated_fraction_disk(epoch):
        """This method computes the approximate illuminated fraction 'k' of the
        disk of the Moon. The method used has a relatively low accuracy, but it
        is enough to the 2nd decimal place.

        :param epoch: Instant to compute the Moon's illuminated fraction of the
            disk, as a py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: The approximate illuminated fraction of the Moon's disk.
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 4, 12.0)
        >>> k = Moon.illuminated_fraction_disk(epoch)
        >>> print(round(k, 2))
        0.68
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Get the time from J2000.0 in Julian centuries
        t = (epoch - JDE2000) / 36525.0
        # Mean elongation of the Moon
        D = 297.8501921 + (445267.1114034
                           + (-0.0018819
                              + (1.0/545868.0 - t/113065000.0) * t) * t) * t
        # Sun's mean anomaly
        M = 357.5291092 + (35999.0502909 + (-0.0001536 + t/24490000.0) * t) * t
        # Moon's mean anomaly
        Mprime = 134.9633964 + (477198.8675055
                                + (0.0087414
                                   + (1.0/69699.9
                                      + t/14712000.0) * t) * t) * t
        # Reduce the angles to a [0 360] range
        D = Angle(Angle.reduce_deg(D)).to_positive()
        Dr = D.rad()
        M = Angle(Angle.reduce_deg(M)).to_positive()
        Mr = M.rad()
        Mprime = Angle(Angle.reduce_deg(Mprime)).to_positive()
        Mprimer = Mprime.rad()
        # Compute the 'i' angle
        i = Angle(180.0 - D - 6.289 * sin(Mprimer) + 2.1 * sin(Mr)
                  - 1.274 * sin(2.0 * Dr - Mprimer) - 0.658 * sin(2.0 * Dr)
                  - 0.214 * sin(2.0 * Mprimer) - 0.11 * sin(Dr))
        k = (1.0 + cos(i.rad())) / 2.0
        return k

    @staticmethod
    def position_bright_limb(epoch):
        """This method computes the position angle of the Moon's bright limb,
        i.e., the position angle of the midpoint of the illuminated limb,
        reckoned eastward from the North Point of the disk (not from the axis
        of rotation of the lunar globe).

        :param epoch: Instant to compute the position angle of the  Moon's
            bright limb, as a py:class:`Epoch` object.
        :type epoch: :py:class:`Epoch`

        :returns: The position angle of the Moon's bright limb.
        :rtype: :py:class:`Angle`
        :raises: TypeError if input value is of wrong type.

        >>> epoch = Epoch(1992, 4, 12.0)
        >>> xi = Moon.position_bright_limb(epoch)
        >>> print(round(xi, 1))
        285.0
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch)):
            raise TypeError("Invalid input type")
        # Compute the right ascension and declination of the Sun
        a0, d0, r0 = Sun.apparent_rightascension_declination_coarse(epoch)
        # Now compute the right ascension and declination of the Moon
        a, d, r, ppi = Moon.apparent_equatorial_pos(epoch)
        a0r = a0.rad()
        d0r = d0.rad()
        ar = a.rad()
        dr = d.rad()
        # Compute the numerator of the tan(xi) formula
        numerator = cos(d0r) * sin(a0r - ar)
        # Now the denominator
        denominator = sin(d0r) * cos(dr) - cos(d0r) * sin(dr) * cos(a0r - ar)
        # Now let's compute xi
        xi = atan2(numerator, denominator)
        xi = Angle(xi, radians=True).to_positive()
        return xi

    @staticmethod
    def moon_phase(epoch, target="new"):
        """This method computes the phases of the moon closest to the provided
        epoch.

        :param epoch: Epoch we want to compute the moon phase for
        :type year: :py:class:`Epoch`
        :param target: Corresponding phase. It can be "new" (New Moon), "first"
            (First Quarter), "full" (Full Moon) and "last" (Last Quarter).
        :type target: str

        :returns: The instant of time when the provided phase happens
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input values are of wrong type.
        :raises: ValueError if 'target' value is invalid.

        >>> epoch = Epoch(1977, 2, 15.0)
        >>> new_moon = Moon.moon_phase(epoch, target="new")
        >>> y, m, d, h, mi, s = new_moon.get_full_date()
        >>> print("{}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))
        1977/2/18 3:37:42.0
        >>> epoch = Epoch(2044, 1, 15.0)
        >>> new_moon = Moon.moon_phase(epoch, target="last")
        >>> y, m, d, h, mi, s = new_moon.get_full_date()
        >>> print("{}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))
        2044/1/21 23:48:17.0
        """

        # First check that input values are of correct types
        if not (isinstance(epoch, Epoch) and isinstance(target, str)):
            raise TypeError("Invalid input types")
        # Second, check that the target is correct
        if (
            (target != "new")
            and (target != "first")
            and (target != "full")
            and (target != "last")
        ):
            raise ValueError("'target' value is invalid")
        # Let's start computing the year with decimals
        y, m, d = epoch.get_date()
        num_days_year = 365.0
        if Epoch.is_leap(y):
            num_days_year = 366.0
        doy = Epoch.get_doy(y, m, d)
        year = y + doy / num_days_year
        # We start computing the 'k' parameter
        k = iint((year - 2000.0) * 12.3685)
        if target == "first":
            k += 0.25
        elif target == "full":
            k += 0.5
        elif target == "last":
            k += 0.75
        t = k / 1236.85
        # Compute the time of the 'mean' phase of the Moon
        jde = (2451550.09766 + 29.530588861 * k
               + (0.00015437 + (-0.00000015 + 0.00000000073 * t) * t) * t * t)
        E = 1.0 + (-0.002516 - 0.0000074 * t) * t
        # Sun's mean anomaly
        M = 2.5534 + 29.1053567 * k + (-0.0000014 - 0.00000011 * t) * t * t
        # Moon's mean anomaly
        Mprime = (201.5643 + 385.81693528 * k
                  + (0.0107582 + (0.00001238 - 0.000000058 * t) * t) * t * t)
        # Moon's argument of latitude
        F = (160.7108 + 390.67050284 * k
             + (-0.0016118 + (-0.00000227 + 0.000000011 * t) * t) * t * t)
        # Longitude of the ascending node of the lunar orbit
        Omega = (124.7746 - 1.56375588 * k
                 + (0.0020672 + 0.00000215 * t) * t * t)
        M = Angle(Angle.reduce_deg(M)).to_positive()
        Mr = M.rad()
        Mprime = Angle(Angle.reduce_deg(Mprime)).to_positive()
        Mprimer = Mprime.rad()
        F = Angle(Angle.reduce_deg(F)).to_positive()
        Fr = F.rad()
        Omega = Angle(Angle.reduce_deg(Omega)).to_positive()
        Omegar = Omega.rad()
        # Planetary arguments
        a1 = 299.77 + 0.107408 * k - 0.009173 * t * t
        a2 = 251.88 + 0.016321 * k
        a3 = 251.83 + 26.651886 * k
        a4 = 349.42 + 36.412478 * k
        a5 = 84.66 + 18.206239 * k
        a6 = 141.74 + 53.303771 * k
        a7 = 207.14 + 2.453732 * k
        a8 = 154.84 + 7.30686 * k
        a9 = 34.52 + 27.261239 * k
        a10 = 207.19 + 0.121824 * k
        a11 = 291.34 + 1.844379 * k
        a12 = 161.72 + 24.198154 * k
        a13 = 239.56 + 25.513099 * k
        a14 = 331.55 + 3.592518 * k
        a1 = Angle(Angle.reduce_deg(a1)).to_positive()
        a1r = a1.rad()
        a2 = Angle(Angle.reduce_deg(a2)).to_positive()
        a2r = a2.rad()
        a3 = Angle(Angle.reduce_deg(a3)).to_positive()
        a3r = a3.rad()
        a4 = Angle(Angle.reduce_deg(a4)).to_positive()
        a4r = a4.rad()
        a5 = Angle(Angle.reduce_deg(a5)).to_positive()
        a5r = a5.rad()
        a6 = Angle(Angle.reduce_deg(a6)).to_positive()
        a6r = a6.rad()
        a7 = Angle(Angle.reduce_deg(a7)).to_positive()
        a7r = a7.rad()
        a8 = Angle(Angle.reduce_deg(a8)).to_positive()
        a8r = a8.rad()
        a9 = Angle(Angle.reduce_deg(a9)).to_positive()
        a9r = a9.rad()
        a10 = Angle(Angle.reduce_deg(a10)).to_positive()
        a10r = a10.rad()
        a11 = Angle(Angle.reduce_deg(a11)).to_positive()
        a11r = a11.rad()
        a12 = Angle(Angle.reduce_deg(a12)).to_positive()
        a12r = a12.rad()
        a13 = Angle(Angle.reduce_deg(a13)).to_positive()
        a13r = a13.rad()
        a14 = Angle(Angle.reduce_deg(a14)).to_positive()
        a14r = a14.rad()
        # Now let's compute the corrections
        corr = 0.0
        w = 0.0
        if target == "new":
            corr = (-0.4072 * sin(Mprimer) + 0.17241 * E * sin(Mr)
                    + 0.01608 * sin(2.0 * Mprimer) + 0.01039 * sin(2.0 * Fr)
                    + 0.00739 * E * sin(Mprimer - Mr)
                    - 0.00514 * E * sin(Mprimer + Mr)
                    + 0.00208 * E * E * sin(2.0 * Mr)
                    - 0.00111 * sin(Mprimer - 2.0 * Fr)
                    - 0.00057 * sin(Mprimer + 2.0 * Fr)
                    + 0.00056 * E * sin(2.0 * Mprimer + Mr)
                    - 0.00042 * sin(3.0 * Mprimer)
                    + 0.00042 * E * sin(Mr + 2.0 * Fr)
                    + 0.00038 * E * sin(Mr - 2.0 * Fr)
                    - 0.00024 * E * sin(2.0 * Mprimer - Mr)
                    - 0.00017 * sin(Omegar) - 0.00007 * sin(Mprimer + 2.0 * Mr)
                    + 0.00004 * sin(2.0 * (Mprimer - Fr))
                    + 0.00004 * sin(3.0 * Mr)
                    + 0.00003 * sin(Mprimer + Mr - 2.0 * Fr)
                    + 0.00003 * sin(2.0 * (Mprimer + Fr))
                    - 0.00003 * sin(Mprimer + Mr + 2.0 * Fr)
                    + 0.00003 * sin(Mprimer - Mr + 2.0 * Fr)
                    - 0.00002 * sin(Mprimer - Mr - 2.0 * Fr)
                    - 0.00002 * sin(3.0 * Mprimer + Mr)
                    + 0.00002 * sin(4.0 * Mprimer))
        elif target == "full":
            corr = (-0.40614 * sin(Mprimer) + 0.17302 * E * sin(Mr)
                    + 0.01614 * sin(2.0 * Mprimer) + 0.01043 * sin(2.0 * Fr)
                    + 0.00734 * E * sin(Mprimer - Mr)
                    - 0.00515 * E * sin(Mprimer + Mr)
                    + 0.00209 * E * E * sin(2.0 * Mr)
                    - 0.00111 * sin(Mprimer - 2.0 * Fr)
                    - 0.00057 * sin(Mprimer + 2.0 * Fr)
                    + 0.00056 * E * sin(2.0 * Mprimer + Mr)
                    - 0.00042 * sin(3.0 * Mprimer)
                    + 0.00042 * E * sin(Mr + 2.0 * Fr)
                    + 0.00038 * E * sin(Mr - 2.0 * Fr)
                    - 0.00024 * E * sin(2.0 * Mprimer - Mr)
                    - 0.00017 * sin(Omegar) - 0.00007 * sin(Mprimer + 2.0 * Mr)
                    + 0.00004 * sin(2.0 * (Mprimer - Fr))
                    + 0.00004 * sin(3.0 * Mr)
                    + 0.00003 * sin(Mprimer + Mr - 2.0 * Fr)
                    + 0.00003 * sin(2.0 * (Mprimer + Fr))
                    - 0.00003 * sin(Mprimer + Mr + 2.0 * Fr)
                    + 0.00003 * sin(Mprimer - Mr + 2.0 * Fr)
                    - 0.00002 * sin(Mprimer - Mr - 2.0 * Fr)
                    - 0.00002 * sin(3.0 * Mprimer + Mr)
                    + 0.00002 * sin(4.0 * Mprimer))
        elif target == "first" or target == "last":
            corr = (-0.62801 * sin(Mprimer) + 0.17172 * E * sin(Mr)
                    - 0.01183 * E * sin(Mprimer + Mr)
                    + 0.00862 * sin(2.0 * Mprimer) + 0.00804 * sin(2.0 * Fr)
                    + 0.00454 * E * sin(Mprimer - Mr)
                    + 0.00204 * E * E * sin(2.0 * Mr)
                    - 0.0018 * sin(Mprimer - 2.0 * Fr)
                    - 0.0007 * sin(Mprimer + 2.0 * Fr)
                    - 0.0004 * sin(3.0 * Mprimer)
                    - 0.00034 * E * sin(2.0 * Mprimer - Mr)
                    + 0.00032 * E * sin(Mr + 2.0 * Fr)
                    + 0.00032 * E * sin(Mr - 2.0 * Fr)
                    - 0.00028 * E * E * sin(Mprimer + 2.0 * Mr)
                    + 0.00027 * E * sin(2.0 * Mprimer + Mr)
                    - 0.00017 * sin(Omegar)
                    - 0.00005 * sin(Mprimer - Mr - 2.0 * Fr)
                    + 0.00004 * sin(2.0 * (Mprimer + Fr))
                    - 0.00004 * sin(Mprimer + Mr + 2.0 * Fr)
                    + 0.00004 * sin(Mprimer - 2.0 * Mr)
                    + 0.00003 * sin(Mprimer + Mr - 2.0 * Fr)
                    + 0.00003 * sin(3.0 * Mr)
                    + 0.00002 * sin(2.0 * (Mprimer - Fr))
                    + 0.00002 * sin(Mprimer - Mr + 2.0 * Fr)
                    - 0.00002 * sin(3.0 * Mprimer + Mr))
            w = (0.00306 - 0.00038 * E * cos(Mr) + 0.00026 * cos(Mprimer)
                 - 0.00002 * cos(Mprimer - Mr) + 0.00002 * cos(Mprimer + Mr)
                 + 0.00002 * cos(2.0 * Fr))
            if target == "last":
                w = -w
        # Additional corrections for all phases
        corr2 = (0.000325 * sin(a1r) + 0.000165 * sin(a2r)
                 + 0.000164 * sin(a3r) + 0.000126 * sin(a4r)
                 + 0.000110 * sin(a5r) + 0.000062 * sin(a6r)
                 + 0.000060 * sin(a7r) + 0.000056 * sin(a8r)
                 + 0.000047 * sin(a9r) + 0.000042 * sin(a10r)
                 + 0.000040 * sin(a11r) + 0.000037 * sin(a12r)
                 + 0.000035 * sin(a13r) + 0.000023 * sin(a14r))
        jde += corr + corr2 + w
        jde = Epoch(jde)
        return jde


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Saturn class
    print("\n" + 35 * "*")
    print("*** Use of Moon class")
    print(35 * "*" + "\n")

    # Let's compute the Moon geocentric ecliiptical position for a given epoch
    epoch = Epoch(1992, 4, 12.0)
    Lambda, Beta, Delta, ppi = Moon.geocentric_ecliptical_pos(epoch)
    print_me("Longitude (Lambda)", round(Lambda, 6))    # 133.162655
    print_me("Latitude (Beta)", round(Beta, 6))         # -3.229126
    print_me("Distance (Delta)", round(Delta, 1))       # 368409.7
    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))  # 0.991990

    print("")

    # Now let's compute the apparent ecliptical position
    epoch = Epoch(1992, 4, 12.0)
    Lambda, Beta, Delta, ppi = Moon.apparent_ecliptical_pos(epoch)
    print_me("Longitude (Lambda)", round(Lambda, 6))    # 133.167264
    print_me("Latitude (Beta)", round(Beta, 6))         # -3.229126
    print_me("Distance (Delta)", round(Delta, 1))       # 368409.7
    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))  # 0.991990

    print("")

    # Get the apparent equatorial position
    epoch = Epoch(1992, 4, 12.0)
    ra, dec, Delta, ppi = Moon.apparent_equatorial_pos(epoch)
    print_me("Right Ascension (ra)", round(ra, 6))      # 134.688469
    print_me("Declination (dec)", round(dec, 6))        # 13.768367
    print_me("Distance (Delta)", round(Delta, 1))       # 368409.7
    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))  # 0.991990

    print("")

    # Compute the longitude of the Moon's mean ascending node
    epoch = Epoch(1913, 5, 27.0)
    Omega = Moon.longitude_mean_ascending_node(epoch)
    print_me("Longitude of the mean ascending node", round(Omega, 1))   # 0.0
    epoch = Epoch(1959, 12, 7.0)
    Omega = Moon.longitude_mean_ascending_node(epoch)
    print_me("Longitude of the mean ascending node", round(Omega, 1))   # 180.0

    print("")

    # Get the longitude of the Moonś true ascending node
    epoch = Epoch(1913, 5, 27.0)
    Omega = Moon.longitude_true_ascending_node(epoch)
    print_me("Longitude of the true ascending node", round(Omega, 4))  # 0.8763

    print("")

    # Compute the longitude of the Moon's mean perigee
    epoch = Epoch(2021, 3, 5.0)
    Pi = Moon.longitude_mean_perigee(epoch)
    print_me("Longitude of the mean perigee", round(Pi, 5))     # 224.89194

    print("")

    # Compute the approximate illuminated fraction of the Moon's disk
    epoch = Epoch(1992, 4, 12.0)
    k = Moon.illuminated_fraction_disk(epoch)
    print_me("Approximate illuminated fraction of Moon's disk", round(k, 2))
    # 0.68

    print("")

    # Compute the position angle of the bright limb of the Moon
    epoch = Epoch(1992, 4, 12.0)
    xi = Moon.position_bright_limb(epoch)
    print_me("Position angle of the bright limb of the Moon", round(xi, 1))
    # 285.0

    print("")

    # Calculate the instant of a New Moon
    epoch = Epoch(1977, 2, 15.0)
    new_moon = Moon.moon_phase(epoch, target="new")
    y, m, d, h, mi, s = new_moon.get_full_date()
    print("New Moon: {}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))
    # 1977/2/18 3:37:42.0

    # Calculate the time of a Last Quarter
    epoch = Epoch(2044, 1, 15.0)
    new_moon = Moon.moon_phase(epoch, target="last")
    y, m, d, h, mi, s = new_moon.get_full_date()
    print("Last Quarter: {}/{}/{} {}:{}:{}".format(y, m, d, h, mi,
                                                   round(s, 0)))
    # 2044/1/21 23:48:17.0


if __name__ == "__main__":

    main()
