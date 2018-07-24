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


from math import pi, degrees, radians


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

        :param a: Semi-major or equatorial radius
        :type a: float
        :param f: Flattening
        :type f: float
        :param omega: Angular velocity of the Earth
        :type omega: float
        """
        self._a = a
        self._f = f
        self._omega = omega


class Earth(object):
    """
    Class Angle deals with angles in either decimal format (d.dd) of in
    sexagesimal format (d m' s'').

    It provides methods to handle an Angle object like it were a simple float,
    but adding the functionality associated with an angle.

    The constructor takes decimals and sexagesimal input. The sexagesimal
    angles can be given as separate degree, minutes, seconds values, or as
    tuples or lists. It is also possible to provide another Angle object as
    input.

    Also, if **radians=True** is passed to the constructor, then the input
    value is considered as in radians, and converted to degrees.
    """

    def __init__(self, *args, **kwargs):
        """Angle constructor.

        It takes decimals and sexagesimal input. The sexagesimal angles can be
        given as separate degree, minutes, seconds values, or as tuples or
        lists. It is also possible to provide another Angle object as input.

        If **radians=True** is passed, then the input value is converted from
        radians to degrees.

        If **ra=True** is passed, then the input value is converted from Right
        Ascension to degrees

        :param \*args: Input angle, in decimal or sexagesimal format, or Angle
        :type \*args: int, float, list, tuple, :py:class:`Angle`
        :param radians: If True, input angle is in radians. False by default.
        :type radians: bool
        :param ra: If True, input angle is in Right Ascension. False by default
        :type ra: bool

        :returns: Angle object.
        :rtype: :py:class:`Angle`
        :raises: TypeError if input values are of wrong type.

        >>> a = Angle(-13, 30, 0.0)
        >>> print(a)
        -13.5
        >>> b = Angle(a)
        >>> print(b)
        -13.5
        """
        self._deg = 0.0         # Angle value is stored here in decimal format
        self._tol = TOL
        self.set(*args, **kwargs)   # Let's use 'set()' method to set angle

    def set(self, *args, **kwargs):
        """Method used to define the value of the Angle object.

        It takes decimals and sexagesimal input. The sexagesimal angles can be
        given as separate degree, minutes, seconds values, or as tuples or
        lists. It is also possible to provide another Angle object as input.

        If **radians=True** is passed, then the input value is converted from
        radians to degrees

        If **ra=True** is passed, then the input value is converted from Right
        Ascension to degrees

        :param \*args: Input angle, in decimal or sexagesimal format, or Angle
        :type \*args: int, float, list, tuple, :py:class:`Angle`
        :param radians: If True, input angle is in radians. False by default.
        :type radians: bool
        :param ra: If True, input angle is in Right Ascension. False by default
        :type ra: bool

        :returns: None.
        :rtype: None
        :raises: TypeError if input values are of wrong type.
        """
        if 'ra' in kwargs:
            if kwargs['ra']:
                # Input values are a Right Ascension
                self.set_ra(*args)
                return
        # If no arguments are given, internal angle is set to zero
        if len(args) == 0:
            self._deg = 0.0
            return
        # If we have only one argument, it can be a single value, a tuple/list
        # or an Angle
        elif len(args) == 1:
            deg = args[0]
            if isinstance(deg, Angle):                  # Copy constructor
                self._deg = deg._deg
                self._tol = deg._tol
                return
            if isinstance(deg, (int, float)):
                if 'radians' in kwargs:
                    if kwargs['radians']:
                        # Input value is in radians. Convert to degrees
                        deg = degrees(deg)
                # This works for ints, floats and Angles
                self._deg = Angle.reduce_deg(deg)
                return
            elif isinstance(deg, (list, tuple)):
                if len(deg) == 0:
                    raise TypeError("Invalid input value")
                elif len(deg) == 1:
                    # This is a single value
                    if 'radians' in kwargs:
                        if kwargs['radians']:
                            # Input value is in radians. Convert to degrees
                            deg[0] = degrees(deg[0])
                    self._deg = Angle.reduce_deg(deg[0])
                    return
                elif len(deg) == 2:
                    # Seconds value is set to zero
                    self._deg = Angle.dms2deg(deg[0], deg[1])
                    return
                elif len(deg) == 3:
                    # The first three values are taken into account
                    self._deg = Angle.dms2deg(deg[0], deg[1], deg[2])
                    return
                else:
                    # Only the first four values are taken into account
                    sign = -1.0 if deg[0] < 0 or deg[1] < 0 or deg[2] < 0 \
                        or deg[3] < 0 else 1.0
                    # If sign < 0, make all values negative, to be sure
                    deg0 = sign * abs(deg[0])
                    deg1 = sign * abs(deg[1])
                    deg2 = sign * abs(deg[2])
                    self._deg = Angle.dms2deg(deg0, deg1, deg2)
                    return
            else:
                raise TypeError("Invalid input value")
        elif len(args) == 2:
            # Seconds value is set to zero
            self._deg = Angle.dms2deg(args[0], args[1])
            return
        elif len(args) == 3:
            # The first three values are taken into account
            self._deg = Angle.dms2deg(args[0], args[1], args[2])
            return
        else:
            # Only the first four values are taken into account
            sign = -1.0 if args[0] < 0 or args[1] < 0 or args[2] < 0 \
                or args[3] < 0 else 1.0
            # If sign < 0, make all values negative, to be sure
            args0 = sign * abs(args[0])
            args1 = sign * abs(args[1])
            args2 = sign * abs(args[2])
            self._deg = Angle.dms2deg(args0, args1, args2)
            return


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Earth class
    print('\n' + 35*'*')
    print("*** Use of Earth class")
    print(35*'*' + '\n')

    # Create an Angle object, providing degrees, minutes and seconds
    a = Angle(-23.0, 26.0, 48.999983999)

    # First we print using the __call__ method (note the extra parentheses)
    print_me("The angle 'a()' is", a())             # -23.44694444

    print("")

    # Use the copy constructor
    b = Angle(a)
    print_me("Angle 'b', which is a copy of 'a', is", b)

    print("")


if __name__ == '__main__':

    main()
