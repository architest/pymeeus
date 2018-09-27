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


from math import sin, cos

# from base import TOL
from Angle import Angle
from Epoch import Epoch, JDE2000


"""
.. module:: Sun
   :synopsis: Module including functions regarding Sun position
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


def sun_true_longitude_coarse(epoch):
    """This function provides the Sun's true longitude with a relatively low
    accuracy of 0.01 degree.

    :param epoch: Epoch to compute the position of the Sun
    :type epoch: :py:class:`Epoch`

    :returns: A tuple containing the true (ecliptical) longitude and the radius
        vector in astronomical units.
    :rtype: tuple
    :raises: TypeError if input value is of wrong type.

    >>> epoch = Epoch(1992, 10, 13)
    >>> true_lon, r = sun_true_longitude_coarse(epoch)
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


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Sun functions
    print('\n' + 35*'*')
    print("*** Use of Sun functions")
    print(35*'*' + '\n')

    # Here follows a series of important parameters related to the angle
    # between Earth's rotation axis and the ecliptic
#    e0 = mean_obliquity(1987, 4, 10)
#    print("The mean angle between Earth rotation axis and ecliptic axis for "+
#          "1987/4/10 is:")
#    print_me("Mean obliquity", e0.dms_str(n_dec=3))         # 23d 26' 27.407''
#    epsilon = true_obliquity(1987, 4, 10)
#    print("'True' (instantaneous) angle between those axes for 1987/4/10 is:")
#    print_me("True obliquity", epsilon.dms_str(n_dec=3))    # 23d 26' 36.849''
#    epsilon = true_obliquity(2018, 7, 29)
#    print("'True' (instantaneous) angle between those axes for 2018/7/29 is:")
#    print_me("True obliquity", epsilon.dms_str(True, 4))    # 23d 26' 7.2157''
#
#    # The nutation effect is separated in two components: One parallel to the
#    # ecliptic (nutation in longitude) and other perpendicular to the ecliptic
#    # (nutation in obliquity)
#    print("Nutation correction in longitude for 1987/4/10:")
#    dpsi = nutation_longitude(1987, 4, 10)
#    print_me("Nutation in longitude", dpsi.dms_str(n_dec=3))  # 0d 0' -3.788''
#    print("Nutation correction in obliquity for 1987/4/10:")
#    depsilon = nutation_obliquity(1987, 4, 10)            # 0d 0' 9.443''
#    print_me("Nutation in obliquity", depsilon.dms_str(n_dec=3))

    print("")


if __name__ == '__main__':

    main()
