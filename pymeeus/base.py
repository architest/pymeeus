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


from math import pi as pi


"""
.. module:: base
   :synopsis: Basic routines used by the pymeeus module
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


DEG2RAD = pi/180.0
"""Constant to convert from degrees to radians."""

RAD2DEG = 180.0/pi
"""Constant to convert from radians to degrees."""

TOL = 1E-10
"""Internal tolerance being used by default"""


class Angle(object):
    """
    Class Angle deals with angles in either decimal format (d.dd) of in
    sexagesimal format (d m' s'').

    It provides methods to handle an Angle object like it were a simple float,
    but adding the functionality associated with an angle.

    The constructor takes decimals and sexagesimal input. The sexagesimal
    angles can be given as separate degree, minutes, seconds values, or as
    tuples or lists.

    Also, if 'radians=True' is passed to the constructor, then the input value
    is considered as in radians, and converted to degrees.
    """

    def __init__(self, *args, **kwargs):
        """Angle constructor.

        It takes decimals and sexagesimal input. The sexagesimal angles can be
        given as separate degree, minutes, seconds values, or as tuples or
        lists.

        If 'radians=True' is passed, then the input value is converted from
        radians to degrees.

        :param \*args: Input angle, in decimal or sexagesimal format.
        :type \*args: int, float, list, tuple
        :param radians: If True, input angle is in radians. False by default.
        :type radians: bool

        :returns: Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(-13, 30, 0.0)
        >>> print(a)
        -13.5
        """
        self._deg = 0.0         # Angle value is stored here in decimal format
        self._tol = TOL
        self.set(*args, **kwargs)   # Let's use 'set()' method to set angle

    @staticmethod
    def reduce_deg(deg):
        """Takes a degree value in decimal format and converts it to a float
        value in the +/-[0:360) range.

        :param deg: Input degree angle in decimal format.
        :type deg: int, float, Angle

        :returns: Float value of the angle in the +/-[0:360) range.
        :rtype: float

        >>> a = 386.3
        >>> b = Angle.reduce_deg(a)
        >>> print(round(b, 1))
        26.3
        """
        if abs(deg) >= 360.0:
            # Extract the sign
            sign = 1.0 if deg >= 0 else -1.0
            frac = abs(deg) % 1         # Separate the fractional part
            deg = int(abs(deg)) % 360   # Reduce to [0:360) range
            deg = sign * (deg + frac)   # Rebuild the value
        return float(deg)

    @staticmethod
    def reduce_dms(degrees, minutes, seconds=0.0):
        """Takes a degree value in sexagesimal format and converts it to a
        value in the +/-[0:360) range (degrees) and [0:60) range (minutes and
        seconds). It also takes care of fractional degrees and minutes.

        :param degrees: Degrees.
        :type degrees: int, float
        :param minutes: Minutes.
        :type minutes: int, float
        :param seconds: Seconds. 0.0 by default.
        :type seconds: int, float

        :returns: Angle in sexagesimal format, with ranges properly adjusted.
        :rtype: tuple

        >>> print(Angle.reduce_dms(-743.0, 26.0, 49.6))
        (23, 26, 49.6, -1.0)
        """
        # If any of the input values is negative, the sign is negative
        sign = -1.0 if (degrees < 0) or (minutes < 0) or (seconds < 0) else 1.0
        degrees = abs(degrees)
        minutes = abs(minutes)
        seconds = abs(seconds)
        # We need to work first from degrees to seconds
        if (degrees % 1 > 0.0):
            # The degrees value has decimals, push them to minutes
            minutes += (degrees % 1) * 60.0
        degrees = int(degrees)                  # Keep the integer part
        if (minutes % 1 > 0.0):
            # The minutes value has decimals, push them to seconds
            seconds += (minutes % 1) * 60.0
        minutes = int(minutes)
        # We now need to work from seconds to degrees, because of overflow
        if (seconds >= 60.0):
            minutes += int(seconds / 60.0)      # Push the excess to minutes
            seconds = seconds % 60              # Keep the rest
        if (minutes >= 60.0):
            degrees += int(minutes / 60.0)      # Push the excess to degrees
            minutes = minutes % 60              # Keep the rest
        degrees = degrees % 360                 # Keep degrees in [0:360) range
        return (degrees, minutes, seconds, sign)

    @staticmethod
    def deg2dms(deg):
        """Converts input from decimal to sexagesimal angle format.

        :param deg: Degrees decimal format.
        :type deg: int, float

        :returns: Angle in sexagesimal format, with ranges adjusted.
        :rtype: tuple

        :note: The output format is (Degrees, Minutes, Seconds, sign)

        >>> print(Angle.deg2dms(23.44694444))
        (23, 26, 48.999983999997596, 1.0)
        """
        deg = Angle.reduce_deg(deg)  # Reduce the degrees to the [0:360) range
        # Extract the sign
        sign = 1.0 if deg >= 0 else -1.0
        # We have the sign, now let's work with positive numbers
        deg = abs(deg)
        mi = (deg % 1)*60.0         # Get the minutes, with decimals
        de = int(deg)               # Get the integer part of the degrees
        se = (mi % 1)*60.0          # Get the seconds
        mi = int(mi)
        return (de, mi, se, sign)

    @staticmethod
    def dms2deg(degrees, minutes, seconds=0.0):
        """Converts an angle from sexagesimal to decimal format.

        :param degrees: Degrees.
        :type degrees: int, float
        :param minutes: Minutes.
        :type minutes: int, float
        :param seconds: Seconds. 0.0 by default.
        :type seconds: int, float

        :returns: Angle in decimal format, within +/-[0:360) range.
        :rtype: float

        >>> print(Angle.dms2deg(-23, 26, 48.999983999997596))
        -23.44694444
        """
        (de, mi, se, sign) = Angle.reduce_dms(degrees, minutes, seconds)
        deg = sign * (de + mi/60.0 + se/3600.0)
        return float(deg)

    def get_tolerance(self):
        """Gets the internal tolerance value used to compare Angles.

        :note: The default tolerance value is TOL.

        :returns: Internal tolerance.
        :rtype: float
        """
        return self._tol

    def set_tolerance(self, tol):
        """Changes the internal tolerance value used to compare Angles.

        :param tol: New tolerance value.
        :type tol: int, float

        :returns: None
        :rtype: None
        """
        self._tol = tol
        return

    def __call__(self):
        """Method used when object is called only with parenthesis.

        :returns: The internal value of the Angle object.
        :rtype: int, float

        >>> a = Angle(54.6)
        >>> print(a())
        54.6
        """
        return self._deg

    def __str__(self):
        """Method used when trying to print the object.

        :returns: Angle as string.
        :rtype: string

        >>> a = Angle(12.5)
        >>> print(a)
        12.5
        """
        return str(self._deg)

    def __repr__(self):
        return "Angle({})".format(self._deg)

    def set(self, *args, **kwargs):
        """Method used to define the value of the Angle object.

        It takes decimals and sexagesimal input. The sexagesimal angles can be
        given as separate degree, minutes, seconds values, or as tuples or
        lists.

        If 'radians=True' is passed, then the input value is converted from
        radians to degrees

        :param \*args: Input angle, in decimal or sexagesimal format.
        :type \*args: int, float, list, tuple
        :param radians: If True, input angle is in radians. False by default.
        :type radians: bool

        :returns: None.
        :rtype: None
        :raises: TypeError if input values are of wrong type.
        """
        # If no arguments are given, internal angle is set to zero
        if len(args) == 0:
            self._deg = 0.0
        # If we have only one argument, it can be a single value or tuple/list
        elif len(args) == 1:
            deg = args[0]
            if isinstance(deg, (int, float, Angle)):
                if 'radians' in kwargs:
                    if kwargs['radians']:
                        # Input value is in radians. Convert to degrees
                        deg = deg*RAD2DEG
                # This works for ints, floats and Angles
                self._deg = Angle.reduce_deg(deg)
            elif isinstance(deg, (list, tuple)):
                if len(deg) == 0:
                    raise TypeError("Invalid input value")
                elif len(deg) == 1:
                    # This is a single value
                    if 'radians' in kwargs:
                        if kwargs['radians']:
                            # Input value is in radians. Convert to degrees
                            deg[0] = deg[0]*RAD2DEG
                    self._deg = Angle.reduce_deg(deg[0])
                elif len(deg) == 2:
                    # Seconds value is set to zero
                    self._deg = Angle.dms2deg(deg[0], deg[1])
                elif len(deg) == 3:
                    # The first three values are taken into account
                    self._deg = Angle.dms2deg(deg[0], deg[1], deg[2])
                else:
                    # Only the first four values are taken into account
                    sign = -1.0 if deg[0] < 0 or deg[1] < 0 or deg[2] < 0 \
                        or deg[3] < 0 else 1.0
                    # If sign < 0, make sure deg[0] is negative
                    degrees = sign * abs(deg[0])
                    self._deg = Angle.dms2deg(degrees, deg[1], deg[2])
            else:
                raise TypeError("Invalid input value")
        elif len(args) == 2:
            # Seconds value is set to zero
            self._deg = Angle.dms2deg(args[0], args[1])
        elif len(args) == 3:
            # The first three values are taken into account
            self._deg = Angle.dms2deg(args[0], args[1], args[2])
        elif len(args) >= 4:
            # Only the first four values are taken into account
            sign = -1.0 if args[0] < 0 or args[1] < 0 or args[2] < 0 \
                or args[3] < 0 else 1.0
            # If sign < 0, make sure args[0] is negative
            degrees = sign * abs(args[0])
            self._deg = Angle.dms2deg(degrees, args[1], args[2])

    def set_radians(self, rads):
        """Method to define the value of the Angle objecti from radians.

        :param rads: Input angle, in radians.
        :type rads: int, float

        :returns: None.
        :rtype: None
        :raises: TypeError if input value is of wrong type.

        >>> a = Angle()
        >>> a.set_radians(pi)
        >>> print(a)
        180.0
        """
        self.set(rads, radians=True)

    def set_ra(self, *args):
        """Define the value of the Angle object from a Right Ascension.

        It takes decimals and sexagesimal input. The sexagesimal Right
        Ascensions can be given as separate hours, minutes, seconds values, or
        as tuples or lists.

        :param \*args: Input Right Ascension, in decimal or sexagesimal format.
        :type \*args: int, float, list, tuple

        :returns: None.
        :rtype: None
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle()
        >>> a.set_ra(9, 14, 55.8)
        >>> print(a)
        138.7325
        """
        self.set(*args)
        self._deg *= 15.0   # Multipy Right Ascension by 15.0 to get degrees

    def dms_str(self, fancy=True):
        """Returns the Angle value as a sexagesimal string.

        The parameter 'fancy' allows to print in "Dd M' S''" format if true,
        and in "D:M:S" (easier to parse) if False.

        :param fancy: Format of output string. True by default.
        :type fancy: bool

        :returns: Angle value as string in sexagesimal format.
        :rtype: string

        >>> a = Angle(42.75)
        >>> print(a.dms_str())
        42d 45' 0.0''
        >>> print(a.dms_str(fancy=False))
        42:45:0.0
        """
        d, m, s, sign = Angle.deg2dms(self._deg)
        if fancy:
            if d != 0:
                return "{}d {}' {}''".format(int(sign*d), m, s)
            elif m != 0:
                return "{}' {}''".format(int(sign*m), s)
            elif s != 0.0:
                return "{}''".format(sign*s)
            else:
                return "0d 0' 0.0''"
        else:
            if d != 0:
                return "{}:{}:{}".format(int(sign*d), m, s)
            elif m != 0:
                return "0:{}:{}".format(int(sign*m), s)
            elif s != 0.0:
                return "0:0:{}".format(sign*s)
            else:
                return "0:0:0.0"

    def get_ra(self):
        """Returns the Angle value as a Right Ascension in float format

        :returns: The internal value of the Angle object as Right Ascension.
        :rtype: int, float

        >>> a = Angle(138.75)
        >>> print(a.get_ra())
        9.25
        """
        return self._deg/15.0

    def ra_str(self, fancy=True):
        """Returns the Angle value as a sexagesimal string in Right Ascension.

        The parameter 'fancy' allows to print in "Hh M' S''" format if true,
        and in "H:M:S" (easier to parse) if False.

        :param fancy: Format of output string. True by default.
        :type fancy: bool

        :returns: Angle value as Right Ascension in sexagesimal format.
        :rtype: string

        >>> a = Angle(138.75)
        >>> print(a.ra_str())
        9h 15' 0.0''
        >>> print(a.ra_str(fancy=False))
        9:15:0.0
        """
        a = Angle(self())/15.0
        s = a.dms_str(fancy)
        if fancy:
            s = s.replace('d', 'h')
        return s

    def rad(self):
        """Returns the Angle value in radians.

        :returns: Angle value in radians.
        :rtype: float

        >>> a = Angle(47.762)
        >>> print(round(a.rad(), 8))
        0.83360416
        """
        return self._deg * DEG2RAD

    def radians(self):
        """Returns the angle value in radians.

        :returns: Angle value in radians.
        :rtype: float

        >>> a = Angle(47.762)
        >>> print(round(a.radians(), 8))
        0.83360416
        """
        return self.rad()

    def to_positive(self):
        """Converts the internal angle value from negative to positive.

        :returns: This angle object.
        :rtype: Angle

        >>> a = Angle(-87.32)
        >>> print(a.to_positive())
        272.68
        """
        if self._deg < 0:
            self._deg = 360.0 - abs(self._deg)
        return self

    def __eq__(self, b):
        """This method defines the 'is equal' operator between Angles.

        :note: For the comparison, the internal tolerance value is used.

        :returns: A boolean.
        :rtype: bool
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(172.01)
        >>> b = Angle(172.009)
        >>> a == b
        False
        """
        if isinstance(b, (int, float)):
            return abs(self._deg - float(b)) < self._tol
        elif isinstance(b, Angle):
            return abs(self._deg - b._deg) < self._tol
        else:
            raise TypeError("Wrong operand type")

    def __ne__(self, b):
        """This method defines the 'is not equal' operator between Angles.

        :note: For the comparison, the internal tolerance value is used.

        :returns: A boolean.
        :rtype: bool

        >>> a = Angle(11.200001)
        >>> b = Angle(11.200000)
        >>> a != b
        True
        """
        return not self.__eq__(b)           # '!=' == 'not(==)'

    def __lt__(self, b):
        """This method defines the 'is less than' operator between Angles.

        :returns: A boolean.
        :rtype: bool
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(72.0)
        >>> b = Angle(72.0)
        >>> a < b
        False
        """
        if isinstance(b, (int, float)):
            return self._deg < float(b)
        elif isinstance(b, Angle):
            return self._deg < b._deg
        else:
            raise TypeError("Wrong operand type")

    def __ge__(self, b):
        """This method defines 'is equal or greater' operator between Angles.

        :returns: A boolean.
        :rtype: bool

        >>> a = Angle(172.01)
        >>> b = Angle(172.009)
        >>> a >= b
        True
        """
        return not self.__lt__(b)           # '>=' == 'not(<)'

    def __gt__(self, b):
        """This method defines the 'is greater than' operator between Angles.

        :returns: A boolean.
        :rtype: bool
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(172.01)
        >>> b = Angle(172.009)
        >>> a > b
        True
        """
        if isinstance(b, (int, float)):
            return self._deg > float(b)
        elif isinstance(b, Angle):
            return self._deg > b._deg
        else:
            raise TypeError("Wrong operand type")

    def __le__(self, b):
        """This method defines 'is equal or less' operator between Angles.

        :returns: A boolean.
        :rtype: bool

        >>> a = Angle(72.0)
        >>> b = Angle(72.0)
        >>> a <= b
        True
        """
        return not self.__gt__(b)           # '<=' == 'not(>)'

    def __neg__(self):
        """This method is used to obtain the negative version of this Angle.

        :returns: A new Angle object.
        :rtype: Angle

        >>> a = Angle(-11.2)
        >>> print(-a)
        11.2
        """
        return Angle(-self._deg)

    def __abs__(self):
        """This method is used to obtain the absolute value of this Angle.

        :returns: A new Angle object.
        :rtype: Angle

        >>> a = Angle(-303.67)
        >>> print(abs(a))
        303.67
        """
        return Angle(abs(self._deg))

    def __mod__(self, b):
        """This method is used to obtain the module b of this Angle.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(333.0)
        >>> b = Angle(72.0)
        >>> print(a % b)
        45.0
        """
        # Negative values will be treated as if they were positive
        sign = 1.0 if self._deg >= 0.0 else -1.0
        if isinstance(b, (int, float)):
            return Angle(sign*(abs(self._deg) % b))
        elif isinstance(b, Angle):
            return Angle(sign*(abs(self._deg) % b._deg))
        else:
            raise TypeError("Wrong operand type")

    def __add__(self, b):
        """This method defines the addition between Angles.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(83.1)
        >>> b = Angle(18.4)
        >>> print(a + b)
        101.5
        """
        if isinstance(b, (int, float)):
            return Angle(self._deg + float(b))
        elif isinstance(b, Angle):
            return Angle(self._deg + b._deg)
        else:
            raise TypeError("Wrong operand type")

    def __sub__(self, b):
        """This method defines the subtraction between Angles.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(25.4)
        >>> b = Angle(10.2)
        >>> print(a - b)
        15.2
        """
        return self.__add__(-b)

    def __mul__(self, b):
        """This method defines the multiplication between Angles.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(33.0)
        >>> b = Angle(72.0)
        >>> print(a * b)
        216.0
        """
        if isinstance(b, (int, float)):
            return Angle(self._deg * float(b))
        elif isinstance(b, Angle):
            return Angle(self._deg * b._deg)
        else:
            raise TypeError("Wrong operand type")

    def __div__(self, b):
        """This method defines the division between Angles.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(172.0)
        >>> b = Angle(86.0)
        >>> print(a/b)
        2.0
        """
        if b == 0.0:
            raise ZeroDivisionError("Division by zero is not allowed")
        if isinstance(b, (int, float)):
            return Angle(self._deg / float(b))
        elif isinstance(b, Angle):
            return Angle(self._deg / b._deg)
        else:
            raise TypeError("Wrong operand type")

    def __truediv__(self, b):
        """This method defines the division between Angles (Python 3).

        :returns: A new Angle object.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.
        :see: __div__
        """
        return self.__div__(b)

    def __pow__(self, b):
        """This method defines the power operation for Angles.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(12.5)
        >>> b = Angle(4.0)
        >>> print(a ** b)
        294.0625
        """
        if isinstance(b, (int, float)):
            return Angle(self._deg**b)
        elif isinstance(b, Angle):
            return Angle(self._deg**b._deg)
        else:
            raise TypeError("Wrong operand type")

    def __imod__(self, b):
        """This method defines the accumulative module b of this Angle.

        :returns: This Angle.
        :rtype: Angle

        >>> a = Angle(330.0)
        >>> b = Angle(45.0)
        >>> a %= b
        >>> print(a)
        15.0
        """
        # Negative values will be treated as if they were positive
        self = self % b
        return self

    def __iadd__(self, b):
        """This method defines the accumulative addition to this Angle.

        :returns: This Angle.
        :rtype: Angle

        >>> a = Angle(172.1)
        >>> b = Angle(54.6)
        >>> a += b
        >>> print(a)
        226.7
        """
        self = self + b
        return self

    def __isub__(self, b):
        """This method defines the accumulative subtraction to this Angle.

        :returns: This Angle.
        :rtype: Angle

        >>> a = Angle(97.0)
        >>> b = Angle(39.0)
        >>> a -= b
        >>> print(a)
        58.0
        """
        self = self - b
        return self

    def __imul__(self, b):
        """This method defines the accumulative multiplication to this Angle.

        :returns: This Angle.
        :rtype: Angle

        >>> a = Angle(30.0)
        >>> b = Angle(55.0)
        >>> a *= b
        >>> print(a)
        210.0
        """
        self = self * b
        return self

    def __idiv__(self, b):
        """This method defines the accumulative division to this Angle.

        :returns: This Angle.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(330.0)
        >>> b = Angle(30.0)
        >>> a /= b
        >>> print(a)
        11.0
        """
        if b == 0.0:
            raise ZeroDivisionError("Division by zero is not allowed")
        if not isinstance(b, (int, float, Angle)):
            raise TypeError("Wrong operand type")
        self = self / b
        return self

    def __itruediv__(self, b):
        """This method defines accumulative division to this Angle (Python3).

        :returns: This Angle.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.
        :see: __idiv__
        """
        return self.__idiv__(b)

    def __ipow__(self, b):
        """This method defines the accumulative power to this Angle.

        :returns: This Angle.
        :rtype: Angle

        >>> a = Angle(37.0)
        >>> b = Angle(3.0)
        >>> a **= b
        >>> print(a)
        253.0
        """
        self = self ** b
        return self

    def __rmod__(self, b):
        """This method defines module operation between Angles by the right.

        :returns: A new Angle object.
        :rtype: Angle

        >>> a = Angle(80.0)
        >>> print(350 % a)
        30.0
        """
        if isinstance(b, (int, float)):
            b = Angle(b)
        # Negative values will be treated as if they were positive
        sign = 1.0 if b._deg >= 0.0 else -1.0
        return Angle(sign*(abs(b._deg) % self._deg))

    def __radd__(self, b):
        """This method defines the addition between Angles by the right

        :returns: A new Angle object.
        :rtype: Angle

        >>> a = Angle(83.1)
        >>> print(8.5 + a)
        91.6
        """
        return self.__add__(b)  # In this case, it is the same as by the left

    def __rsub__(self, b):
        """This method defines the subtraction between Angles by the right.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(25.0)
        >>> print(24.0 - a)
        -1.0
        """
        return -self.__sub__(b)    # b - a = -(a - b)

    def __rmul__(self, b):
        """This method defines multiplication between Angles by the right.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(11.0)
        >>> print(250.0 * a)
        230.0
        """
        return self.__mul__(b)  # In this case, it is the same as by the left

    def __rdiv__(self, b):
        """This method defines division between Angles by the right.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(80.0)
        >>> print(350 / a)
        4.375
        """
        if self == 0.0:
            raise ZeroDivisionError("Division by zero is not allowed")
        if isinstance(b, (int, float)):
            return Angle(float(b) / self._deg)
        elif isinstance(b, Angle):
            return Angle(b._deg / self._deg)
        else:
            raise TypeError("Wrong operand type")

    def __rtruediv__(self, b):
        """This method defines division between Angle by the right (Python3).

        :returns: A new Angle object.
        :rtype: Angle
        :raises: ZeroDivisionError if divisor is zero.
        :raises: TypeError if input values are of wrong type.
        :see: __rdiv__
        """
        return self.__rdiv__(b)

    def __rpow__(self, b):
        """This method defines the power operation for Angles by the right.

        :returns: A new Angle object.
        :rtype: Angle
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(5.0)
        >>> print(24.0 ** a)
        144.0
        """
        if isinstance(b, (int, float)):
            return Angle(b**self._deg)
        elif isinstance(b, Angle):
            return Angle(b._deg**self._deg)
        else:
            raise TypeError("Wrong operand type")

    def __float__(self):
        """This method returns Angle value as a float.

        :returns: Internal angle value as a float.
        :rtype: float

        >>> a = Angle(213.8)
        >>> float(a)
        213.8
        """
        return float(self._deg)

    def __int__(self):
        """This method returns Angle value as an int.

        :returns: Internal angle value as an int.
        :rtype: int

        >>> a = Angle(213.8)
        >>> int(a)
        213
        """
        return int(self._deg)

    def __round__(self, n=0):
        """This method returns an Angle with content rounded to 'n' decimal.

        :returns: A new Angle object.
        :rtype: Angle

        >>> a = Angle(11.4361)
        >>> print(round(a, 2))
        11.44
        """
        # NOTE: This method is only called in Python 3
        return Angle(round(self._deg, n))


class Interpolation(object):
    """
    Class Interpolation deals with finding intermediate values to those given
    in a table.

    Besides the basic interpolation, this class also provides methods to find
    the extrema (maximum or minimum value) in a given interval (if present),
    and it also has methods to find roots (the values of the argument 'x' for
    which the function 'y' becomes zero).

    Please note that it seems that Meeus uses the Bessel interpolation method
    (Chapter 3). However, this class uses the Newton interpolation method
    because it is (arguably) more flexible regarding the number of entries
    provided in the interpolation table.

    The constructor takes pairs of (x, y) values from the table of interest.
    These pairs of values can be given as a sequence of int/floats, tuples,
    lists or Angles.

    :note: When using Angles, be careful with the 360-to-0 discontinuity.

    If a sequence of int/floats is given, the values in the odd positions are
    considered to belong to the 'x' set, while the values in the even positions
    belong to the 'y' set. If only one tuple or list is provided, it is assumed
    that it is the 'y' set, and the 'x' set is build from 0 onwards with steps
    of length 1.

    Please keep in mind that a minimum of two data pairs are needed in order to
    carry out any interpolation. If only one value is provided, a ValueError
    exception will be raised.
    """

    def __init__(self, *args):
        """Interpolation constructor.

        This takes pairs of (x, y) values from the table of interest. These
        pairs of values can be given as a sequence of int/floats, tuples, or
        lists.

        If a sequence of int, floats or Angles is given, the values in the odd
        positions are considered to belong to the 'x' set, while the values in
        the even positions belong to the 'y' set. If only one tuple or list is
        provided, it is assumed that it is the 'y' set, and the 'x' set is
        build from 0 onwards with steps of length 1.

        Please keep in mind that a minimum of two data pairs are needed in
        order to carry out any interpolation. If only one value is provided, a
        ValueError exception will be raised.

        :param \*args: Input tabular values.
        :type \*args: int, float, list, tuple, Angle

        :returns: Interpolation object.
        :rtype: Interpolation
        :raises: ValueError if not enough input data pairs are provided.
        :raises: TypeError if input values are of wrong type.

        >>> i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
        >>> print(i._x)
        [1, 2, 3, 4, 5, 6]
        >>> print(i._y)
        [2, 4, 6, 8, 10, 12]
        >>> j = Interpolation([3, -8, 1, 12, 2, 5, 8])
        >>> print(j._x)
        [0, 1, 2, 3, 4, 5, 6]
        >>> print(j._y)
        [3, -8, 1, 12, 2, 5, 8]
        >>> k = Interpolation(3, -8, 1, 12, 2, 5, 8)
        >>> print(k._x)
        [1, 2, 3]
        >>> print(k._y)
        [12, 5, -8]
        """
        self._x = []
        self._y = []
        self._table = []
        self.set(*args)         # Let's use 'set()' method to handle the setup

    def _order_points(self):
        """Method to put in ascending order (w.r.t. 'x') the data points."""

        # Let's work with copies of the original lists
        x = list(self._x)
        y = list(self._y)

        xnew = []
        ynew = []
        xmax = max(x) + 1.0
        for _ in range(len(x)):
            # Get the index of the minimum value
            imin = x.index(min(x))
            # Append the current minimum value to the new 'x' list
            xnew.append(x[imin])
            # Store the *position* of the current minimum value to new 'y' list
            ynew.append(imin)
            # The current minimum value will no longer be the minimum
            x[imin] = xmax

        # In the new 'y' list, substitute the positions by the real values
        for i in range(len(x)):
            ynew[i] = y[ynew[i]]

        # Store the results in the corresponding fields
        self._x = xnew
        self._y = ynew

    def set(self, *args):
        """Method used to define the value pairs of Interpolation object.

        This takes pairs of (x, y) values from the table of interest. These
        pairs of values can be given as a sequence of int/floats, tuples,
        lists, or Angles.

        :note: When using Angles, be careful with the 360-to-0 discontinuity.

        If a sequence of int, floats or Angles is given, the values in the odd
        positions are considered to belong to the 'x' set, while the values in
        the even positions belong to the 'y' set. If only one tuple or list is
        provided, it is assumed that it is the 'y' set, and the 'x' set is
        build from 0 onwards with steps of length 1.

        Please keep in mind that a minimum of two data pairs are needed in
        order to carry out any interpolation. If only one value is provided, a
        ValueError exception will be raised.

        :param \*args: Input tabular values.
        :type \*args: int, float, list, tuple, Angle

        :returns: None.
        :rtype: None
        :raises: ValueError if not enough input data pairs are provided.
        :raises: TypeError if input values are of wrong type.

        >>> i = Interpolation()
        >>> i.set([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
        >>> print(i._x)
        [1, 2, 3, 4, 5, 6]
        >>> print(i._y)
        [2, 4, 6, 8, 10, 12]
        >>> j = Interpolation()
        >>> j.set([3, -8, 1, 12, 2, 5, 8])
        >>> print(j._x)
        [0, 1, 2, 3, 4, 5, 6]
        >>> print(j._y)
        [3, -8, 1, 12, 2, 5, 8]
        >>> k = Interpolation(3, -8, 1, 12, 2, 5, 8)
        >>> print(k._x)
        [1, 2, 3]
        >>> print(k._y)
        [12, 5, -8]
        """
        # Clean up the internal data tables
        self._x = []
        self._y = []
        self._table = []
        # If no arguments are given, return. Internal data tables are empty
        if len(args) == 0:
            return
        # If we have only one argument, it can be a single value or tuple/list
        elif len(args) == 1:
            if isinstance(args[0], (int, float, Angle)):
                # Insuficient data to interpolate. Raise ValueError exception
                raise ValueError("Invalid number of input values")
            elif isinstance(args[0], (list, tuple)):
                seq = args[0]
                if len(seq) < 2:
                    raise ValueError("Invalid number of input values")
                else:
                    # Read input values into 'y', and create 'x'
                    i = 0
                    for value in seq:
                        self._x.append(i)
                        self._y.append(value)
                        i += 1
            else:
                raise TypeError("Invalid input value")
        elif len(args) == 2:
            if isinstance(args[0], (int, float, Angle)) or \
                    isinstance(args[1], (int, float, Angle)):
                # Insuficient data to interpolate. Raise ValueError exception
                raise ValueError("Invalid number of input values")
            elif isinstance(args[0], (list, tuple)) and \
                    isinstance(args[1], (list, tuple)):
                x = args[0]
                y = args[1]
                # Check if they have the same length. If not, make them equal
                length_min = min(len(x), len(y))
                x = x[:length_min]
                y = y[:length_min]
                if len(x) < 2 or len(y) < 2:
                    raise ValueError("Invalid number of input values")
                else:
                    # Read input values into 'x' and 'y'
                    for xval, yval in zip(x, y):
                        self._x.append(xval)
                        self._y.append(yval)
            else:
                raise TypeError("Invalid input value")
        elif len(args) == 3:
            # In this case, no combination of input values is valid
            raise ValueError("Invalid number of input values")
        else:
            # If there is an odd number of arguments, drop the last one
            if len(args) % 2 != 0:
                args = args[:-1]
            # Check that all the arguments are ints, floats or Angles
            all_numbers = True
            for arg in args:
                all_numbers = all_numbers and \
                            isinstance(arg, (int, float, Angle))
            # If any of the values failed the test, raise an exception
            if not all_numbers:
                raise TypeError("Invalid input value")
            # Now, extract the data: Odds are x's, evens are y's
            for i in range(int(len(args)/2.0)):
                self._x.append(args[2 * i])
                self._y.append(args[2 * i + 1])
        # After self._x is found, confirm that x's are different to each other
        for i in range(len(self._x) - 1):
            for k in range(i + 1, len(self._x)):
                if abs(self._x[i] - self._x[k]) < TOL:
                    raise ValueError("Invalid input: Values in 'x' are equal")
        # Order the data points if needed
        self._order_points()
        # Create table containing Newton coefficientes, only if values given
        if len(self._x) > 0:
            self._compute_table()

    def __str__(self):
        """Method used when trying to print the object.

        :returns: Internal tabular values as strings.
        :rtype: string

        >>> i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
        >>> print(i)
        X: [1, 2, 3, 4, 5, 6]
        Y: [2, 4, 6, 8, 10, 12]
        """
        xstr = "X: " + str(self._x) + "\n"
        ystr = "Y: " + str(self._y)
        return xstr + ystr

    def _compute_table(self):
        """Method to compute coefficients of Newton interpolation method."""
        for i in range(len(self._x)):
            self._table.append(self._newton_diff(0, i))

    def _newton_diff(self, start, end):
        """Auxiliary method to compute the elements of the Newton table.

        :param start: Starting index
        :type start: int
        :param end: Ending index
        :type end: int

        :returns: Resulting value of the element of the Newton table.
        :rtype: float
        """
        if abs(end - start) < TOL:
            val = self._y[start]
        else:
            x = list(self._x)       # Let's make a copy, just in case
            val = (self._newton_diff(start, end - 1) -
                   self._newton_diff(start + 1, end)) / (x[start] - x[end])

        return val

    def __call__(self, x):
        """Method to interpolate the function at a given 'x'.

        :param x: Point where the interpolation will be carried out.
        :type x: int, float, Angle

        :returns: Resulting value of the interpolation.
        :rtype: float
        :raises: ValueError if input value is outside of interpolation range.
        :raises: TypeError if input value is of wrong type.

        >>> i = Interpolation([7, 8, 9], [0.884226, 0.877366, 0.870531])
        >>> y = round(i(8.18125), 6)
        >>> print(y)
        0.876125
        """
        # Check if input value is of correct type
        if isinstance(x, (int, float, Angle)):
            # Check if 'x' already belongs to the data table
            for i in range(len(self._x)):
                if abs(x - self._x[i]) < TOL:
                    return self._y[i]           # We don't need to look further
            # Check if Newton coefficients table is not empty
            if len(self._table) == 0:
                raise RuntimeError("Internal table is empty. Use set().")
            # Check that x is within interpolation table values
            if x < self._x[0] or x > self._x[-1]:
                raise ValueError("Input value outside of interpolation range.")
            # Horner's method is used to efficiently compute the result
            val = self._table[-1]
            for i in range(len(self._table) - 1, 0, -1):
                val = self._table[i - 1] + (x - self._x[i - 1]) * val

            return val
        else:
            raise TypeError("Invalid input value")

    def derivative(self, x):
        """Method to compute the derivative from interpolation polynomial.

        :param x: Point where the interpolation derivative will be carried out.
        :type x: int, float, Angle

        :returns: Resulting value of the interpolation derivative.
        :rtype: float
        :raises: ValueError if input value is outside of interpolation range.
        :raises: TypeError if input value is of wrong type.

        >>> m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])
        >>> m.derivative(-1.0)
        8.0
        >>> m.derivative(0.5)
        -1.0
        """
        # Check if input value is of correct type
        if isinstance(x, (int, float, Angle)):
            # Check that x is within interpolation table values
            if x < self._x[0] or x > self._x[-1]:
                raise ValueError("Input value outside of interpolation range.")
            # If we only have two interpolation points, derivative is simple
            if len(self._x) == 2:
                return (self._y[1] - self._y[0])/(self._x[1] - self._x[0])
            else:
                res = self._table[1]
                for k in range(len(self._table) - 1, 1, -1):
                    val = 0.0
                    for j in range(k):
                        s = 1.0
                        for i in range(k):
                            if i != j:
                                s *= (x - self._x[i])
                        val += s
                    res += val*self._table[k]
                return res
        else:
            raise TypeError("Invalid input value")

    def root(self, xl=0, xh=0, max_iter=1000):
        """Method to find the root inside the [xl, xh] range.

        This method applies, in principle, the Newton method to find the root;
        however, if conditions are such that Newton method may not bei properly
        behaving or converging, then it switches to the linear Interpolation
        method.

        If values xl, xh are not given, the limits of the interpolation table
        values will be used.

        :note: This method returns a ValueError exception if the corresponding
        yl = f(xl) and yh = f(xh) values have the same sign. In that case, the
        method assumes there is no root in the [xl, xh] interval.

        :note: If any of the xl, xh values is beyond the limits given by the
        interpolation values, its value will be set to the corresponding limit.

        :note: If xl == xh (and not zero), a ValueError exception is raised.

        :note: If the method doesn't converge within max_iter ierations, then a
        ValueError exception is raised.

        :param xl: Lower limit of interval where the root will be looked for.
        :type xl: int, float, Angle
        :param xh: Higher limit of interval where the root will be looked for.
        :type xh: int, float, Angle
        :param max_iter: Maximum number of iterations allowed.
        :type max_iter: int

        :returns: Root of the interpolated function within [xl, xh] interval.
        :rtype: int, float, Angle
        :raises: ValueError if yl = f(xl), yh = f(xh) have same sign.
        :raises: ValueError if xl == xh.
        :raises: ValueError if maximum number of iterations is exceeded.
        :raises: TypeError if input value is of wrong type.

        >>> m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])
        >>> round(m.root(), 8)
        -0.72075922
        """
        # Get the limits of the interpolation table
        xmin = self._x[0]
        xmax = self._x[-1]
        # Check if input value is of correct type
        if isinstance(xl, (int, float, Angle)) and \
           isinstance(xh, (int, float, Angle)):
            # Check if BOTH values are zero
            if xl == 0 and xh == 0:
                xl = xmin
                xh = xmax
            # Check if limits are equal
            if abs(xl - xh) < TOL:
                raise ValueError("Invalid limits: xl and xh values are equal")
            # Put limits in order. Swap them if necessary
            if xl > xh:
                xl, xh = xh, xl
            # Check limits against interpolation table. Reset if necessary
            if xl < self._x[0]:
                xl = xmin
            if xh < self._x[-1]:
                xh = xmax
            yl = self.__call__(xl)
            yh = self.__call__(xh)
            # Check for a couple special cases
            if abs(yl) < TOL:
                return xl               # xl is a root
            if abs(yh) < TOL:
                return xh               # xh is a root
            # Check if signs of ordinates are the same
            if (yl * yh) > 0.0:
                raise ValueError("Invalid interval: Probably no root exists")
            # We are good to go. First option: Newton's root-finding method
            x = (xl + xh)/2.0           # Start in the middle of interval
            y = self.__call__(x)
            num_iter = 0                # Count the number of iterations
            while abs(y) > TOL:
                if num_iter >= max_iter:
                    raise ValueError("Too many iterations: Probably no root\
                                     exists")
                num_iter += 1
                yp = self.derivative(x)
                # If derivative is too small, switch to linear interpolation
                if abs(yp) < 1e-3:
                    x = (xl*yh - xh*yl)/(yh - yl)
                    y = self.__call__(x)
                else:
                    x = x - y/yp
                    # Check if x is within limits
                    if x < xmin or x > xmax:
                        # Switch to linear interpolation
                        x = (xl*yh - xh*yl)/(yh - yl)
                        y = self.__call__(x)
                    else:
                        y = self.__call__(x)
                if (y * yl) >= 0.0:
                    xl = x
                    yl = y
                else:
                    xh = x
                    yh = y
            return x
        else:
            raise TypeError("Invalid input value")

    def minmax(self, xl=0, xh=0, max_iter=1000):
        """Method to find the minimum or maximum inside the [xl, xh] range.

        Finding the minimum or maximum (extremum) of a function within a given
        interval is akin to find the root of its derivative. Therefore, this
        method creates an interpolation object for the derivative function, and
        calls the root method of that object. See root() method for more
        details.

        If values xl, xh are not given, the limits of the interpolation table
        values will be used.

        :note: This method returns a ValueError exception if the corresponding
        derivatives yl' = f'(xl) and yh' = f'(xh) values have the same sign.
        In that case, the method assumes there is no extremum in the [xl, xh]
        interval.

        :note: If any of the xl, xh values is beyond the limits given by the
        interpolation values, its value will be set to the corresponding limit.

        :note: If xl == xh (and not zero), a ValueError exception is raised.

        :note: If the method doesn't converge within max_iter ierations, then
        a ValueError exception is raised.

        :param xl: Lower limit of interval where a extremum will be looked for.
        :type xl: int, float, Angle
        :param xh: Higher limit of interval where extremum will be looked for.
        :type xh: int, float, Angle
        :param max_iter: Maximum number of iterations allowed.
        :type max_iter: int

        :returns: Extremum of interpolated function within [xl, xh] interval.
        :rtype: int, float, Angle
        :raises: ValueError if yl = f(xl), yh = f(xh) have same sign.
        :raises: ValueError if xl == xh.
        :raises: ValueError if maximum number of iterations is exceeded.
        :raises: TypeError if input value is of wrong type.

        >>> m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])
        >>> round(m.minmax(), 8)
        0.33333333
        """
        # Compute the derivatives for the current data
        x = list(self._x)
        y = []
        for xi in x:
            y.append(self.derivative(xi))
        # Create a new Interpolation object
        prime = Interpolation(x, y)
        # Find the root within that object, and return it
        return prime.root(xl, xh, max_iter)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Angle class
    print('\n' + 35*'*')
    print("*** Use of Angle class")
    print(35*'*' + '\n')

    # Create an Angle object, providing degrees, minutes and seconds
    a = Angle(-23.0, 26.0, 48.999983999)

    # First we print using the __call__ method (note the extra parentheses)
    print_me("The angle 'a()' is", a())             # -23.44694444

    # Second we print using the __str__ method (no extra parentheses needed)
    print_me("The angle 'a' is", a)                 # -23.44694444

    print("")

    # Use the static 'deg2dms()' method to carry out conversions
    d, m, s, sign = Angle.deg2dms(23.44694444)
    val = "{}d {}' {}''".format(int(sign*d), m, s)
    print_me("{Deg}d {Min}' {Sec}''", val)          # 23d 26' 48.999984''

    # We can print Angle 'a' directly in sexagesimal format
    # In 'fancy' format:
    print_me("{Deg}d {Min}' {Sec}''", a.dms_str())  # 23d 26' 48.999984''
    # In plain format:
    print_me("{Deg}:{Min}:{Sec}", a.dms_str(False))  # -23:26:48.999983999

    print("")

    # Redefine Angle 'a' several times
    a.set(-0.44694444)
    print("a.set(-0.44694444)")
    print_me("   a.dms_str()", a.dms_str())             # -26' 48.999984''
    a.set(0, 0, -46.31)
    print("a.set(0, 0, -46.31)")
    print_me("   a.dms_str(False)", a.dms_str(False))   # 0:0:-46.31

    print("")

    # We can use decimals in degrees/minutes. They are converted automatically
    a.set(0, -46.25, 0.0)
    print("a.set(0, -46.25, 0.0)")
    print_me("   a.dms_str()", a.dms_str())             # -46' 15.0''
    a.set(0, 0, 0.0)
    print("a.set(0, 0, 0.0)")
    print_me("   a.dms_str()", a.dms_str())             # 0d 0' 0.0''

    print("")

    # We can define the angle as in radians. It will be converted to degrees
    b = Angle(pi, radians=True)
    print_me("b = Angle(pi, radians=True); print(b)", b)    # 180.0

    # And we can easily carry out the 'degrees to radians' conversion
    print_me("print(b.radians())", b.radians())             # 3.14159265359

    print("")

    # We can also specify the angle as a Right Ascension
    print("Angle can be given as a Right Ascension: Hours, Minutes, Seconds")
    a.set_ra(9, 14, 55.8)
    print("a.set_ra(9, 14, 55.8)")
    print_me("   print(a)", a)

    print("")

    # We can print the Angle as Right Ascension, as a float and as string
    a = Angle(138.75)
    print("a = Angle(138.75)")
    print_me("   print(a.get_ra())", a.get_ra())
    print_me("   print(a.ra_str())", a.ra_str())
    print_me("   print(a.ra_str(False))", a.ra_str(False))

    print("")

    # Use the 'to_positive()' method to get the positive version of an angle
    a = Angle(-87.32)                                       # 272.68
    print("a = Angle(-87.32)")
    print_me("   print(a.to_positive())", a.to_positive())

    print("")

    # Call the __repr__() method to get a string defining the current object
    # This string can then be fed to 'eval()' function to generate the object
    print_me("print(b.__repr__())", b.__repr__())           # Angle(180.0)
    c = eval(repr(b))
    print_me("c = eval(repr(b)); print(c)", c)              # 180.0

    print("")

    print_me("c", c)                                        # 180.0

    # Negate an angle
    d = Angle(13, 30)
    print_me("d", d)                                        # 13.5
    e = -d
    print_me("   e = -d", e)                                # -13.5

    # Get the absolute value of an angle
    e = abs(e)
    print_me("   e = abs(e)", e)                            # 13.5

    # Module operation on an angle
    d = Angle(17.0)
    print_me("d", d)                                        # 17.0
    e = c % d
    print_me("   e = c % d", e)                             # 10.0

    print("")

    # Convert the angle to an integer
    d = Angle(13.95)
    print_me("d", d)                                        # 13.95
    print_me("   int(d)", int(d))                           # 13.0
    d = Angle(-4.95)
    print_me("d", d)                                        # -4.95
    print_me("   int(d)", int(d))                           # -4.0

    # Convert the angle to a float
    print_me("   float(d)", float(d))                       # -4.95

    # Round the angle to a float
    e = Angle(-4.951648)
    print_me("e", e)                                        # -4.951648
    print_me("   round(e)", round(e))                       # -5.0
    print_me("   round(e, 2)", round(e, 2))                 # -4.95
    print_me("   round(e, 3)", round(e, 3))                 # -4.952
    print_me("   round(e, 4)", round(e, 4))                 # -4.9516

    print("")

    # Comparison operators
    print_me("   d == e", d == e)                           # False
    print_me("   d != e", d != e)                           # True
    print_me("   d > e", d > e)                             # True
    print_me("   c >= 180.0", c >= 180.0)                   # True
    print_me("   c < 180.0", c < 180.0)                     # False
    print_me("   c <= 180.0", c <= 180.0)                   # True

    print("")

    # It is very easy to add Angles to obtain a new Angle
    e = c + d
    print_me("   c + d", e)                                 # 193.5

    # We can also directly add a decimal angle
    e = c + 11.5
    print_me("   c + 11.5", e)                              # 193.5

    print("")

    # Types allowed are int, float and Angle
    print('e = c + "32.5"')
    try:
        e = c + "32.5"
    except TypeError:
        print("TypeError!: Valid types are int, float, and Angle, not string!")

    print("")

    # Subtraction
    e = c - d
    print_me("   c - d", e)                                 # 184.95

    # Multiplication
    c.set(150.0)
    d.set(5.0)
    print_me("c", c)                                        # 150.0
    print_me("d", d)                                        # 5.0
    e = c * d
    print_me("   c * d", e)                                 # 30.0

    # Division
    c.set(150.0)
    d.set(6.0)
    print_me("d", d)                                        # 6.0
    e = c / d
    print_me("   c / d", e)                                 # 25.0

    print("")

    # Division by zero is not allowed
    d.set(0.0)
    print_me("d", d)                                        # 6.0
    print('e = c / d')
    try:
        e = c / d
    except ZeroDivisionError:
        print("ZeroDivisionError!: Division by zero is not allowed!")

    print("")

    # Power
    d.set(2.2)
    print_me("d", d)                                        # 2.2
    e = c ** d
    print_me("   c ** d", e)                                # 91.57336709992524

    print("")

    # Accumulative module operation
    d.set(17.0)
    print_me("d", d)                                        # 17.0
    e %= d
    print_me("   e %= d", e)                                # 6.573367099925235

    # Accumulative addition
    c += d
    print_me("   c += d", c)                                # 156.0

    # Accumulative subtraction
    print_me("b", b)                                        # 180.0
    c -= b
    print_me("   c -= b", c)                                # -24.0

    # Accumulative multiplication
    print_me("b", b)                                        # 180.0
    c *= b
    print_me("   c *= b", c)                                # -180.0

    # Accumulative division
    print_me("b", b)                                        # 180.0
    d.set(6.0)
    print_me("d", d)                                        # 6.0
    b /= d
    print_me("   b /= d", b)                                # 30.0

    # Accumulative power
    d.set(2.2)
    print_me("d", d)                                        # 2.2
    c = abs(c)
    print_me("   c = abs(c)", c)                            # 180.0
    c **= d
    print_me("   c **= d", c)                               # 97.5978160305

    print("")

    # The same operation, but by the right side
    e = 3.5 + b
    print_me("   e = 3.5 + b", e)                           # 33.5
    e = 3.5 - b
    print_me("   e = 3.5 - b", e)                           # -26.5
    e = 3.5 * b
    print_me("   e = 3.5 * b", e)                           # 105.0
    e = 3.5 / b
    print_me("   e = 3.5 / b", e)                           # 0.116666666667
    e = 3.5 ** b
    print_me("   e = 3.5 ** b", e)                          # 220.0

    # Let's now work with the Interpolation class
    print('\n' + 35*'*')
    print("*** Use of Interpolation class")
    print(35*'*' + '\n')

    i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
    print("i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])")
    print(i)
    print("NOTE:")
    print("   a. They are ordered in 'x'")
    print("   b. The extra value in 'x' was dropped")
    print("")

    j = Interpolation([0.0, 1.0, 3.0], [-1.0, -2.0, 2.0])
    print("j = Interpolation([0.0, 1.0, 3.0], [-1.0, -2.0, 2.0])")
    print(j)
    print_me("j(2)", j(2))
    print_me("j(0.5)", j(0.5))
    # Test with a value already in the data table
    print_me("j(1)", j(1))
    print("")

    # We can interpolate Angles too
    k = Interpolation([27.0, 27.5, 28.0, 28.5, 29.0],
                      [Angle(0, 54, 36.125), Angle(0, 54, 24.606),
                       Angle(0, 54, 15.486), Angle(0, 54, 8.694),
                       Angle(0, 54, 4.133)])
    print("k = Interpolation([27.0, 27.5, 28.0, 28.5, 29.0],\n\
                      [Angle(0, 54, 36.125), Angle(0, 54, 24.606),\n\
                       Angle(0, 54, 15.486), Angle(0, 54, 8.694),\n\
                       Angle(0, 54, 4.133)])")

    print_me("k(28.27777778)", Angle(k(28.1388888889)).dms_str())
    print("")

    m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])
    print("m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])")
    print(m)
    # Get interpolated values
    print_me("m(-0.5)", m(-0.5))
    print_me("m(0.5)", m(0.5))
    # Get derivatives
    print_me("m'(-1.0)", m.derivative(-1.0))
    print_me("m'(-0.5)", m.derivative(-0.5))
    print_me("m'(0.0)", m.derivative(0.0))
    print_me("m'(0.5)", m.derivative(0.5))
    print_me("m'(1.0)", m.derivative(1.0))
    # Get the root within the interval
    print_me("m.root()", m.root())
    # Get the extremum within the interval
    print_me("m.minmax()", m.minmax())


if __name__ == '__main__':

    main()
