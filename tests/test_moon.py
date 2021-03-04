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


from pymeeus.base import TOL
from pymeeus.Moon import Moon
from pymeeus.Epoch import Epoch


# Moon class

def test_moon_geocentric_ecliptical_pos():
    """Tests the method 'geocentric_ecliptical_pos()' of Moon class"""

    epoch = Epoch(1992, 4, 12.0)
    Lambda, Beta, Delta, ppi = Moon.geocentric_ecliptical_pos(epoch)
    Lambda = round(Lambda, 6)
    Beta = round(Beta, 6)
    Delta = round(Delta, 1)
    ppi = round(ppi, 6)

    assert abs(Lambda - 133.162655) < TOL, \
        "ERROR: 1st 'geocentric_ecliptical_pos()' test, 'Lambda' value doesn't\
            match"

    assert abs(Beta - (-3.229126)) < TOL, \
        "ERROR: 2nd 'geocentric_ecliptical_pos()' test, 'Beta' value doesn't\
            match"

    assert abs(Delta - 368409.7) < TOL, \
        "ERROR: 3rd 'geocentric_ecliptical_pos()' test, 'Delta' value doesn't\
            match"

    assert abs(ppi - 0.991990) < TOL, \
        "ERROR: 4th 'geocentric_ecliptical_pos()' test, 'ppi' value doesn't\
            match"


def test_moon_apparent_ecliptical_pos():
    """Tests the method 'apparent_ecliptical_pos()' of Moon class"""

    epoch = Epoch(1992, 4, 12.0)
    Lambda, Beta, Delta, ppi = Moon.apparent_ecliptical_pos(epoch)
    Lambda = round(Lambda, 5)
    Beta = round(Beta, 6)
    Delta = round(Delta, 1)
    ppi = round(ppi, 6)

    assert abs(Lambda - 133.16726) < TOL, \
        "ERROR: 1st 'apparent_ecliptical_pos()' test, 'Lambda' value doesn't\
            match"

    assert abs(Beta - (-3.229126)) < TOL, \
        "ERROR: 2nd 'apparent_ecliptical_pos()' test, 'Beta' value doesn't\
            match"

    assert abs(Delta - 368409.7) < TOL, \
        "ERROR: 3rd 'apparent_ecliptical_pos()' test, 'Delta' value doesn't\
            match"

    assert abs(ppi - 0.991990) < TOL, \
        "ERROR: 4th 'apparent_ecliptical_pos()' test, 'ppi' value doesn't\
            match"


def test_moon_apparent_equatorial_pos():
    """Tests the method 'apparent_equatorial_pos()' of Moon class"""

    epoch = Epoch(1992, 4, 12.0)
    ra, dec, Delta, ppi = Moon.apparent_equatorial_pos(epoch)
    ra = round(ra, 6)
    dec = round(dec, 6)
    Delta = round(Delta, 1)
    ppi = round(ppi, 6)

    assert abs(ra - 134.688469) < TOL, \
        "ERROR: 1st 'apparent_equatorial_pos()' test, 'ra' value doesn't\
            match"

    assert abs(dec - 13.768367) < TOL, \
        "ERROR: 2nd 'apparent_equatorial_pos()' test, 'dec' value doesn't\
            match"

    assert abs(Delta - 368409.7) < TOL, \
        "ERROR: 3rd 'apparent_equatorial_pos()' test, 'Delta' value doesn't\
            match"

    assert abs(ppi - 0.991990) < TOL, \
        "ERROR: 4th 'apparent_equatorial_pos()' test, 'ppi' value doesn't\
            match"


def test_moon_longitude_mean_ascending_node():
    """Tests the method 'longitude_mean_ascending_node()' of Moon class"""

    epoch = Epoch(1913, 5, 27.0)
    Omega1 = Moon.longitude_mean_ascending_node(epoch)
    Omega1 = round(Omega1, 1)
    epoch = Epoch(2043, 9, 10.0)
    Omega2 = Moon.longitude_mean_ascending_node(epoch)
    Omega2 = round(Omega2, 1)
    epoch = Epoch(1959, 12, 7.0)
    Omega3 = Moon.longitude_mean_ascending_node(epoch)
    Omega3 = round(Omega3, 1)
    epoch = Epoch(2108, 11, 3.0)
    Omega4 = Moon.longitude_mean_ascending_node(epoch)
    Omega4 = round(Omega4, 1)

    assert abs(Omega1 - 0.0) < TOL, \
        "ERROR: 1st 'longitude_mean_ascending_node()' test, 'Omega1' value\
            doesn't match"

    assert abs(Omega2 - 0.0) < TOL, \
        "ERROR: 2nd 'longitude_mean_ascending_node()' test, 'Omega2' value\
            doesn't match"

    assert abs(Omega3 - 180.0) < TOL, \
        "ERROR: 3rd 'longitude_mean_ascending_node()' test, 'Omega3' value\
            doesn't match"

    assert abs(Omega4 - 180.0) < TOL, \
        "ERROR: 4th 'longitude_mean_ascending_node()' test, 'Omega4' value\
            doesn't match"


def test_moon_longitude_true_ascending_node():
    """Tests the method 'longitude_true_ascending_node()' of Moon class"""

    epoch = Epoch(1913, 5, 27.0)
    Omega = Moon.longitude_true_ascending_node(epoch)
    Omega = round(Omega, 4)

    assert abs(Omega - 0.8763) < TOL, \
        "ERROR: 1st 'longitude_true_ascending_node()' test, 'Omega' value\
            doesn't match"


def test_moon_longitude_mean_perigee_node():
    """Tests the method 'longitude_mean_perigee()' of Moon class"""

    epoch = Epoch(2021, 3, 5.0)
    Pi = Moon.longitude_mean_perigee(epoch)
    Pi = round(Pi, 5)

    assert abs(Pi - 224.89194) < TOL, \
        "ERROR: 1st 'longitude_mean_perigee()' test, 'Pi' value doesn't match"
