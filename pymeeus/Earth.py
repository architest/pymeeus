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


from math import sqrt


"""
.. module:: Earth
   :synopsis: Class to model Earth's globe
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Ellipsoid(object):
    """
    Class Ellipsoid is useful to encapsulate the most important parameters of
    a given reference ellipsoid.
    """

    def __init__(self, a, f, omega):
        """Ellipsoid constructor.

        :param a: Semi-major or equatorial radius, in meters
        :type a: float
        :param f: Flattening
        :type f: float
        :param omega: Angular velocity of the Earth, in rad/s
        :type omega: float
        """

        self._a = a
        self._f = f
        self._omega = omega

    def __str__(self):
        """Method used when trying to print the object.

        :returns: Semi-major equatorial radius, flattening and angular velocity
           as string.
        :rtype: string

        >>> a = Ellipsoid(6378140.0, 0.0033528132, 7.292e-5)
        >>> print(a)
        6378140.0:0.0033528132:7.292e-05
        """

        return "{}:{}:{}".format(self._a, self._f, self._omega)

    def b(self):
        """Method to return the semi-minor radius.

        :returns: Semi-minor radius, in meters
        :rtype: float

        >>> a = Ellipsoid(6378140.0, 0.0033528132, 7.292e-5)
        >>> round(a.b(), 3)
        6356755.288
        """

        return self._a*(1.0 - self._f)

    def e(self):
        """Method to return the eccentricity of the Earth's meridian.

        :returns: Eccentricity of the Earth's meridian
        :rtype: float

        >>> a = Ellipsoid(6378140.0, 0.0033528132, 7.292e-5)
        >>> round(a.e(), 8)
        0.08181922
        """

        f = self._f
        return sqrt(2.0*f - f*f)


IAU76 = Ellipsoid(6378140.0, (1.0/298.257), 7.292114992e-5)
"""Reference ellipsoid defined by the International Astronomic Union in 1976"""


WGS84 = Ellipsoid(6378137.0, (1.0/298.257223563), 7292115e-11)
"""Reference ellipsoid World Geodetic System 1984, a modern ellipsoid used by
the GPS system, and the standard in many applications"""


class Earth(object):
    """
    Class Earth models the figure of the Earth surface and, with the help of a
    configurable reference ellipsoid, provides a set of handy method to compute
    different parameters, like the distance between two points on the surface.

    Please note that here we depart a little bit from Meeus' book because the
    Earth class uses the **World Geodetic System 1984 (WGS84)** as the default
    reference ellipsoid, instead of the International Astronomical Union 1974,
    which Meeus uses. This change is done because WGS84 is regarded as more
    modern.
    """

    def __init__(self, ellipsoid=WGS84):
        """Earth constructor.

        It takes a reference ellipsoid as input. If not provided, the ellipsoid
        used is the WGS84 by default.

        :param ellipsoid: Reference ellipsoid to be used. WGS84 by default.
        :type radians: :class:`Ellipsoid`

        :returns: Earth object.
        :rtype: :py:class:`Earth`
        :raises: TypeError if input value is of wrong type.

        #>>> a = Angle(-13, 30, 0.0)
        #>>> print(a)
        #-13.5
        #>>> b = Angle(a)
        #>>> print(b)
        #-13.5
        """

        # Set an invalid ellipsoid by default
        self._ellip = Ellipsoid(0.0, 0.0, 0.0)
        self.set(ellipsoid)   # Let's use 'set()' method

    def set(self, ellipsoid):
        """Method used to define an Earth object.

        It takes a reference ellipsoid as input. If not provided, the ellipsoid
        used is the WGS84 by default.

        :param ellipsoid: Reference ellipsoid to be used. WGS84 by default.
        :type radians: :class:`Ellipsoid`

        :returns: None
        :rtype: None
        :raises: TypeError if input value is of wrong type.
        """

        if isinstance(ellipsoid, Ellipsoid):
            self._ellip = ellipsoid
        else:
            raise TypeError("Invalid input value")
        return


def main():

    pass
#     # Let's define a small helper function
#     def print_me(msg, val):
#         print("{}: {}".format(msg, val))
#
#     # Let's show some uses of Earth class
#     print('\n' + 35*'*')
#     print("*** Use of Earth class")
#     print(35*'*' + '\n')
#
#     # Create an Angle object, providing degrees, minutes and seconds
#     a = Angle(-23.0, 26.0, 48.999983999)
#
#     # First we print using the __call__ method (note the extra parentheses)
#     print_me("The angle 'a()' is", a())             # -23.44694444
#
#     print("")
#
#     # Use the copy constructor
#     b = Angle(a)
#     print_me("Angle 'b', which is a copy of 'a', is", b)
#
#     print("")


if __name__ == '__main__':

    main()
