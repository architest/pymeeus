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
