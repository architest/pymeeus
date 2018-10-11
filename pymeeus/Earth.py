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
from Epoch import Epoch
from Coordinates import geometric_vsop_pos, apparent_vsop_pos


"""
.. module:: Earth
   :synopsis: Class to model Earth's globe
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


VSOP87_L = [
    # L0
    [[175347046., 0.0, 0.0], [3341656., 4.6692568, 6283.07585],
     [34894., 4.6261, 12556.1517], [3497., 2.7441, 5753.3849],
     [3418., 2.8289, 3.5231], [3136., 3.6277, 77713.7715],
     [2676., 4.4181, 7860.4194], [2343., 6.1352, 3930.2097],
     [1324., 0.7425, 11506.7698], [1273., 2.0371, 529.691],
     [1199., 1.1096, 1577.3435], [990., 5.233, 5884.927],
     [902., 2.045, 26.298], [857., 3.508, 398.149], [780., 1.179, 5223.694],
     [753., 2.533, 5507.553], [505., 4.583, 18849.228], [492., 4.205, 775.523],
     [357., 2.92, 0.067], [317., 5.849, 11790.629], [284., 1.899, 796.298],
     [271., 0.315, 10977.079], [243., 0.345, 5486.778],
     [206., 4.806, 2544.314], [205., 1.869, 5573.143], [202., 2.458, 6069.777],
     [156., 0.833, 213.299], [132., 3.411, 2942.463], [126., 1.083, 20.775],
     [115., 0.645, 0.98], [103., 0.636, 4694.003], [102., 0.976, 15720.839],
     [102., 4.267, 7.114], [99., 6.21, 2146.17], [98., 0.68, 155.42],
     [86., 5.98, 161000.69], [85., 1.3, 6275.96], [85., 3.67, 71430.7],
     [80., 1.81, 17260.15], [79., 3.04, 12036.46], [75., 1.76, 5088.63],
     [74., 3.5, 3154.69], [74., 4.68, 801.82], [70., 0.83, 9437.76],
     [62., 3.98, 8827.39], [61., 1.82, 7084.9], [57., 2.78, 6286.6],
     [56., 4.39, 14143.5], [56., 3.47, 6279.55], [52., 0.19, 12139.55],
     [52., 1.33, 1748.02], [51., 0.28, 5856.48], [49., 0.49, 1194.45],
     [41., 5.37, 8429.24], [41., 2.4, 19651.05], [39., 6.17, 10447.39],
     [37., 6.04, 10213.29], [37., 2.57, 1059.38], [36., 1.71, 2352.87],
     [36., 1.78, 6812.77], [33., 0.59, 17789.85], [30., 0.44, 83996.85],
     [30., 2.74, 1349.87], [25., 3.16, 4690.48]
     ],
    # L1
    [[628331966747., 0.0, 0.0], [206059., 2.678235, 6283.07585],
     [4303., 2.6351, 12566.1517], [425., 1.59, 3.523], [119., 5.796, 26.298],
     [109., 2.966, 1577.344], [93., 2.59, 18849.23], [72., 1.14, 529.69],
     [68., 1.87, 398.15], [67., 4.41, 5507.55], [59., 2.89, 5223.69],
     [56., 2.17, 155.42], [45., 0.4, 796.3], [36., 0.47, 775.52],
     [29., 2.65, 7.11], [21., 5.34, 0.98], [19., 1.85, 5486.78],
     [19., 4.97, 213.3], [17., 2.99, 6275.96], [16., 0.03, 2544.31],
     [16., 1.43, 2146.17], [15., 1.21, 10977.08], [12., 2.83, 1748.02],
     [12., 3.26, 5088.63], [12., 5.27, 1194.45], [12., 2.08, 4694.0],
     [11., 0.77, 553.57], [10., 1.3, 6286.6], [10., 4.24, 1349.87],
     [9., 2.7, 242.73], [9., 5.64, 951.72], [8., 5.3, 2352.87],
     [6., 2.65, 9437.76], [6., 4.67, 4690.48]
     ],
    # L2
    [[52919., 0.0, 0.0], [8720., 1.0721, 6283.0758], [309., 0.867, 12566.152],
     [27., 0.05, 3.52], [16., 5.19, 26.3], [16., 3.68, 155.42],
     [10., 0.76, 18849.23], [9., 2.06, 77713.77], [7., 0.83, 775.52],
     [5., 4.66, 1577.34], [4., 1.03, 7.11], [4., 3.44, 5573.14],
     [3., 5.14, 796.3], [3., 6.05, 5507.55], [3., 1.19, 242.73],
     [3., 6.12, 529.69], [3., 0.31, 398.15], [3., 2.28, 553.57],
     [2., 4.38, 5223.69], [2., 3.75, 0.98]
     ],
    # L3
    [[289., 5.844, 6283.076], [35., 0.0, 0.0], [17., 5.49, 12566.15],
     [3., 5.2, 155.42], [1., 4.72, 3.52], [1., 5.3, 18849.23],
     [1., 5.97, 242.73]
     ],
    # L4
    [[114., 3.142, 0.0], [8., 4.13, 6283.08], [1., 3.84, 12566.15]
     ],
    # L5
    [[1., 3.14, 0.0]
     ]]
"""This table contains Earth's most important periodic terms from the planetary
theory VSOP87 for the heliocentric longitude. In Meeus' book these values can
be found in pages 418-420."""


VSOP87_B = [
    # B0
    [[280., 3.199, 84334.662], [102., 5.422, 5507.553], [80., 3.88, 5223.69],
     [44., 3.7, 2352.87], [32., 4.0, 1577.34]
     ],
    # B1
    [[9., 3.9, 5507.55], [6., 1.73, 5223.69]
     ]]
"""This table contains Earth's most important periodic terms from the planetary
theory VSOP87 for the heliocentric latitude. In Meeus' book these values can
be found in page 420."""


VSOP87_R = [
    # R0
    [[100013989., 0.0, 0.0], [1670700., 3.0984635, 6283.07585],
     [13956., 3.05525, 12566.1517], [3084., 5.1985, 77713.7715],
     [1628., 1.1739, 5753.3849], [1576., 2.8469, 7860.4194],
     [925., 5.453, 11506.77], [542., 4.564, 3930.21], [472., 3.661, 5884.927],
     [346., 0.964, 5507.553], [329., 5.9, 5223.694], [307., 0.299, 5573.143],
     [243., 4.273, 11790.629], [212., 5.847, 1577.344],
     [186., 5.022, 10977.079], [175., 3.012, 18849.228],
     [110., 5.055, 5486.778], [98.0, 0.89, 6069.78], [86., 5.69, 15720.84],
     [86., 1.27, 161000.69], [65., 0.27, 17260.15], [63., 0.92, 529.69],
     [57., 2.01, 83996.85], [56., 5.24, 71430.7], [49., 3.25, 2544.31],
     [47., 2.58, 775.52], [45., 5.54, 9437.76], [43., 6.01, 6275.96],
     [39., 5.36, 4694.0], [38., 2.39, 8827.39], [37., 0.83, 19651.05],
     [37., 4.9, 12139.55], [36., 1.67, 12036.46], [35., 1.84, 2942.46],
     [33., 0.24, 7084.9], [32., 0.18, 5088.63], [32., 1.78, 398.15],
     [28., 1.21, 6286.6], [28., 1.9, 6279.55], [26., 4.59, 10447.39]
     ],
    # R1
    [[103019., 1.10749, 6283.07585], [1721., 1.0644, 12566.1517],
     [702., 3.142, 0.0], [32.0, 1.02, 18849.23], [31., 2.84, 5507.55],
     [25., 1.32, 5223.69], [18., 1.42, 1577.34], [10., 5.91, 10977.08],
     [9., 1.42, 6275.96], [9., 0.27, 5486.78]
     ],
    # R2
    [[4359., 5.7846, 6283.0758], [124., 5.579, 12566.152], [12., 3.14, 0.0],
     [9., 3.63, 77713.77], [6., 1.87, 5573.14], [3., 5.47, 18849.23]
     ],
    # R3
    [[145., 4.273, 6283.076], [7., 3.92, 12566.15]
     ],
    # R4
    [[4., 2.56, 6283.08]
     ]]
"""This table contains Earth's most important periodic terms from the planetary
theory VSOP87 for the radius vector. In Meeus' book these values can be found
in pages 420-421."""


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

        .. note:: We will consider that positions 'East' and 'South' are
            negative.

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

    def geometric_heliocentric_position(self, epoch, toFK5=True):
        """"This method computes the geometric heliocentric position of the
        Earth for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Earth position, as an Epoch object
        :type epoch: :py:class:`Epoch`
        :param toFK5: Whether or not the small correction to convert to the FK5
            system will be applied or not
        :type toFK5: bool

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> e = Earth()
        >>> epoch = Epoch(1992, 10, 13.0)
        >>> lon, lat, r = e.geometric_heliocentric_position(epoch)
        >>> print(round(lon.to_positive(), 6))
        19.905991
        >>> print(lat.dms_str(n_dec=3))
        -0.621''
        >>> print(round(r, 8))
        0.99760775
        """

        # NOTE: In page 169, Meeus gives a different value for the LONGITUDE
        # (19.907373 degrees) as the one presented above (19.905991 degrees).
        # After many checks and tests, I came to the conclusion that the result
        # above is the right one, and Meeus' result is wrong

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Second, call auxiliary function in charge of computations
        return geometric_vsop_pos(epoch, VSOP87_L, VSOP87_B, VSOP87_R, toFK5)

    def apparent_heliocentric_position(self, epoch):
        """"This method computes the apparent heliocentric position of the
        Earth for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Earth position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> e = Earth()
        >>> epoch = Epoch(1992, 10, 13.0)
        >>> lon, lat, r = e.apparent_heliocentric_position(epoch)
        >>> print(round(lon.to_positive(), 6))
        19.904705
        >>> print(lat.dms_str(n_dec=3))
        -0.621''
        >>> print(round(r, 8))
        0.99760775
        """

        # First check that input values are of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Second, call auxiliary function in charge of computations
        return apparent_vsop_pos(epoch, VSOP87_L, VSOP87_B, VSOP87_R)


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

    print("")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = e.geometric_heliocentric_position(epoch)
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat.dms_str(n_dec=3))
    print_me("Radius vector", r)

    print("")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = e.apparent_heliocentric_position(epoch)
    print_me("Apparent Heliocentric Longitude", lon.to_positive())
    print_me("Apparent Heliocentric Latitude", lat.dms_str(n_dec=3))
    print_me("Radius vector", r)


if __name__ == '__main__':

    main()
