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


from math import sqrt, radians, sin, cos, tan, atan
from Angle import Angle


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
           as a string.
        :rtype: string

        >>> a = Ellipsoid(6378140.0, 0.0033528132, 7.292e-5)
        >>> print(a)
        6378140.0:0.0033528132:7.292e-05
        """

        return "{}:{}:{}".format(self._a, self._f, self._omega)

    def __repr__(self):
        """Method providing the 'official' string representation of the object.
        It provides a valid expression that could be used to recreate the
        object.

        :returns: As string with a valid expression to recreate the object
        :rtype: string

        >>> a = Ellipsoid(6378140.0, 0.0033528132, 7.292e-5)
        >>> repr(a)
        'Ellipsoid(6378140.0, 0.0033528132, 7.292e-05)'
        """

        return "{}({}, {}, {})".format(self.__class__.__name__, self._a,
                                       self._f, self._omega)

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

    def __str__(self):
        """Method used when trying to print the Earth object. It essentially
        returns the corresponting '__str__()' method from the reference
        ellipsoid being used.

        :returns: Semi-major equatorial radius, flattening and angular velocity
           of the current reference ellipsoid, as a string.
        :rtype: string

        >>> e = Earth()
        >>> s = str(e)
        >>> v = s.split(':')
        >>> print(v[0] + '|' + str(round(float(v[1]), 14)) + '|' + v[2] )
        6378137.0|0.00335281066475|7.292115e-05
        """

        return str(self._ellip)

    def __repr__(self):
        """Method providing the 'official' string representation of the object.
        It provides a valid expression that could be used to recreate the
        object.

        :returns: As string with a valid expression to recreate the object
        :rtype: string
        """

        return "{}(ellipsoid=Ellipsoid({}, {}, {}))".format(
            self.__class__.__name__, self._ellip._a, self._ellip._f,
            self._ellip._omega)

    def rho(self, latitude):
        """"Method to compute the rho term, which is the observer distance to
        the center of the Earth, when the observer is at sea level. In this
        case, the Earth's equatorial radius is taken as unity.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`

        :returns: Rho: Distance to the center of the Earth from sea level. It
            is a ratio with respect to Earth equatorial radius.
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.rho(0.0), 1)
        1.0
        """

        if not isinstance(latitude, (int, float, Angle)):
            raise TypeError("Invalid input value")
        if isinstance(latitude, (int, float)):
            phi = radians(latitude)         # Convert to radians
        else:
            phi = latitude.rad()            # It is an Angle. Call method rad()
        return 0.9983271 + 0.0016764*cos(2.0*phi) - 0.0000035*cos(4.0*phi)

    def rho_sinphi(self, latitude, height):
        """"Method to compute the rho*sin(phi') term, needed in the calculation
        of diurnal parallaxes, eclipses and occulatitions.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`
        :param height: Height of the observer above the sea level, in meters
        :type height: int, float

        :returns: rho*sin(phi') term
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> lat = Angle(33, 21, 22.0)
        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.rho_sinphi(lat, 1706), 6)
        0.546861
        """

        if not (isinstance(latitude, (int, float, Angle)) and
                isinstance(height, (int, float))):
            raise TypeError("Invalid input value")
        if isinstance(latitude, (int, float)):
            phi = radians(latitude)         # Convert to radians
        else:
            phi = latitude.rad()            # It is an Angle. Call method rad()
        b_a = self._ellip.b()/self._ellip._a
        u = atan(b_a*tan(phi))
        return b_a * sin(u) + height/self._ellip._a * sin(phi)

    def rho_cosphi(self, latitude, height):
        """"Method to compute the rho*cos(phi') term, needed in the calculation
        of diurnal parallaxes, eclipses and occulatitions.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`
        :param height: Height of the observer above the sea level, in meters
        :type height: int, float

        :returns: rho*cos(phi') term
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> lat = Angle(33, 21, 22.0)
        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.rho_cosphi(lat, 1706), 6)
        0.836339
        """

        if not (isinstance(latitude, (int, float, Angle)) and
                isinstance(height, (int, float))):
            raise TypeError("Invalid input value")
        if isinstance(latitude, (int, float)):
            phi = radians(latitude)         # Convert to radians
        else:
            phi = latitude.rad()            # It is an Angle. Call method rad()
        b_a = self._ellip.b()/self._ellip._a
        u = atan(b_a*tan(phi))
        return cos(u) + height/self._ellip._a * cos(phi)

    def rp(self, latitude):
        """"Method to compute the radius of the parallel circle at the given
        latitude.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`

        :returns: Radius of the parallel circle at given latitude, in meters
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.rp(42.0), 1)
        4747001.2
        """

        if not isinstance(latitude, (int, float, Angle)):
            raise TypeError("Invalid input value")
        if isinstance(latitude, (int, float)):
            phi = radians(latitude)         # Convert to radians
        else:
            phi = latitude.rad()            # It is an Angle. Call method rad()
        a = self._ellip._a
        e = self._ellip.e()
        return (a*cos(phi))/sqrt(1.0 - e*e*sin(phi)*sin(phi))

    def linear_velocity(self, latitude):
        """"Method to compute the linear velocity of a point at latitude, due
        to the rotation of the Earth.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`

        :returns: Linear velocity of a point at latitude, in meters per second
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.linear_velocity(42.0), 2)
        346.16
        """

        if not isinstance(latitude, (int, float, Angle)):
            raise TypeError("Invalid input value")
        omega = self._ellip._omega
        return omega*self.rp(latitude)

    def rm(self, latitude):
        """"Method to compute the radius of curvature of the Earth's meridian
        at the given latitude.

        :param latitude: Geodetical or geographical latitude of the observer,
            in degrees
        :type latitude: int, float, :class:`Angle`

        :returns: Radius of curvature of the Earth's meridian at the given
            latitude, in meters
        :rtype: float
        :raises: TypeError if input value is of wrong type.

        >>> e = Earth(ellipsoid=IAU76)
        >>> round(e.rm(42.0), 1)
        6364033.3
        """

        if not isinstance(latitude, (int, float, Angle)):
            raise TypeError("Invalid input value")
        if isinstance(latitude, (int, float)):
            phi = radians(latitude)         # Convert to radians
        else:
            phi = latitude.rad()            # It is an Angle. Call method rad()
        a = self._ellip._a
        e = self._ellip.e()
        return (a*(1.0 - e*e))/(1.0 - e*e*sin(phi)*sin(phi))**1.5

    def distance(self, lon1, lat1, lon2, lat2):
        """"This method computes the distance between two points on the Earth's
        surface using the method from H. Andoyer.

        :param lon1: Longitude of the first point, in degrees
        :type lon1: int, float, :class:`Angle`
        :param lat1: Geodetical or geographical latitude of the first point,
            in degrees
        :type lat1: int, float, :class:`Angle`
        :param lon2: Longitude of the second point, in degrees
        :type lon2: int, float, :class:`Angle`
        :param lat2: Geodetical or geographical latitude of the second point,
            in degrees
        :type lat2: int, float, :class:`Angle`

        :returns: Tuple with distance between the two points along Earth's
            surface, and approximate error, in meters
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> e = Earth(ellipsoid=IAU76)
        >>> lon1 = Angle(-2, 20, 14.0)
        >>> lat1 = Angle(48, 50, 11.0)
        >>> lon2 = Angle(77, 3, 56.0)
        >>> lat2 = Angle(38, 55, 17.0)
        >>> dist, error = e.distance(lon1, lat1, lon2, lat2)
        >>> round(dist, 0)
        6181628.0
        >>> error
        69.0
        >>> lon1 = Angle(-2.09)
        >>> lat1 = Angle(41.3)
        >>> lon2 = Angle(73.99)
        >>> lat2 = Angle(40.75)
        >>> dist, error = e.distance(lon1, lat1, lon2, lat2)
        >>> round(dist, 0)
        6176760.0
        >>> error
        69.0
        """

        if not (isinstance(lon1, (int, float, Angle)) and
                isinstance(lat1, (int, float, Angle)) and
                isinstance(lon2, (int, float, Angle)) and
                isinstance(lat2, (int, float, Angle))):
            raise TypeError("Invalid input value")
        if isinstance(lon1, (int, float)):
            l1 = radians(lon1)              # Convert to radians
        else:
            l1 = lon1.rad()                 # It is an Angle. Call method rad()
        if isinstance(lat1, (int, float)):
            phi1 = radians(lat1)            # Convert to radians
        else:
            phi1 = lat1.rad()               # It is an Angle. Call method rad()
        if isinstance(lon2, (int, float)):
            l2 = radians(lon2)              # Convert to radians
        else:
            l2 = lon2.rad()                 # It is an Angle. Call method rad()
        if isinstance(lat2, (int, float)):
            phi2 = radians(lat2)            # Convert to radians
        else:
            phi2 = lat2.rad()               # It is an Angle. Call method rad()
        f = (phi1 + phi2)/2.0
        g = (phi1 - phi2)/2.0
        lam = (l1 - l2)/2.0
        sin2g = sin(g)**2
        cos2g = cos(g)**2
        cos2f = cos(f)**2
        sin2f = sin(f)**2
        sin2lam = sin(lam)**2
        cos2lam = cos(lam)**2
        s = sin2g*cos2lam + cos2f*sin2lam
        c = cos2g*cos2lam + sin2f*sin2lam
        omega = atan(sqrt(s/c))
        r = sqrt(s*c)/omega
        d = 2.0*omega*self._ellip._a
        h1 = (3.0*r - 1.0)/(2.0*c)
        h2 = (3.0*r + 1.0)/(2.0*s)
        fe = self._ellip._f
        dist = d*(1.0 + fe*(h1*sin2f*cos2g - h2*cos2f*sin2g))
        error = round(dist*fe*fe, 0)
        return dist, error


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Earth class
    print('\n' + 35*'*')
    print("*** Use of Earth class")
    print(35*'*' + '\n')

    # An important concept are the reference ellipsoids, comprising information
    # about the Earth global model we are going to use.

    # A very important reference ellipsoid is WGS84, predefined here
    print_me("WGS84", WGS84)
    # First field is equatorial radius, second field is the flattening, and the
    # third field is the angular rotation velocity, in radians per second

    # Let's print the semi-minor axis (polar radius)
    print_me("Polar radius, b", WGS84.b())

    # And now, let's print the eccentricity of Earth's meridian
    print_me("Eccentricity, e", WGS84.e())

    print("")

    # We create an Earth object with a given reference ellipsoid. By default,
    # it is WGS84, but we can use another
    e = Earth(IAU76)
    print("e = Earth(IAU76)")

    # Print the parameters of reference ellipsoid being used
    print_me("'e' Earth object parameters", e)

    print("")

    # Compute the distance to the center of the Earth from a given point at sea
    # level, and at a certain latitude. It is given as a fraction of equatorial
    # radius
    lat = Angle(65, 45, 30.0)               # We can use an Angle for this
    print_me("Relative distance to Earth's center, from latitude 65d 45' 30''",
             e.rho(lat))

    print("")

    # Parameters rho*sin(lat) and rho*cos(lat) are useful for different
    # astronomical applications
    height = 650.0
    print_me("rho*sin(lat)", round(e.rho_sinphi(lat, height), 6))
    print_me("rho*cos(lat)", round(e.rho_cosphi(lat, height), 6))

    print("")

    # Compute the radius of the parallel circle at given latitude
    print_me("Radius of parallel circle at latitude 65d 45' 30'' (meters)",
             round(e.rp(lat), 1))

    # Compute the radius of curvature of the Earth's meridian at given latitude
    print_me("Radius of Earth's meridian at latitude 65d 45' 30'' (meters)",
             round(e.rm(lat), 1))

    print("")

    # It is easy to compute the linear velocity at different latitudes
    print_me("Linear velocity at the Equator (meters/second)",
             round(e.linear_velocity(0.0), 3))
    print_me("Linear velocity at latitude 65d 45' 30'' (meters/second)",
             round(e.linear_velocity(lat), 3))

    print("")

    # Now let's compute the distance between two points on the Earth:
    # Bangkok:          13d 14' 09'' North, 100d 29' 39'' East
    # Buenos Aires:     34d 36' 12'' South,  58d 22' 54'' West
    # NOTE: We will consider that positions 'East' and 'South' are negative

    # Here we will take advantage of facilities provided by Angle class
    lon_ban = Angle(-100, 29, 39.0)
    lat_ban = Angle(13, 14, 9.0)
    lon_bai = Angle(58, 22, 54.0)
    lat_bai = Angle(-34, 36, 12.0)
    dist, error = e.distance(lon_ban, lat_ban, lon_bai, lat_bai)
    print_me("The distance between Bangkok and Buenos Aires is (km)",
             round(dist/1000.0, 2))
    print_me("The approximate error of the estimation is (meters)",
             round(error, 0))


if __name__ == '__main__':

    main()
