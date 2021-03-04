
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


from math import sin, cos, asin
from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch, JDE2000
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
            py:class:`Epoch` object
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
            py:class:`Epoch` object
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
            py:class:`Epoch` object
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
            py:class:`Epoch` object
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
            py:class:`Epoch` object
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
        orbitn in degrees, for a given instant, measured from the mean equinoxi
        of the date.

        :param epoch: Instant to compute the Moon's mean perigee, as an
            py:class:`Epoch` object
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
        Pi = 83.3532465 + (4069.0137287
                           + (-0.01032
                              + (-1.0/80053.0 + t/18999000.0) * t) * t) * t
        Pi = Angle(Pi)
        return Pi


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


if __name__ == "__main__":

    main()
