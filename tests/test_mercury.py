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
from pymeeus.Mercury import Mercury
from pymeeus.Epoch import Epoch


# Mercury class

def test_mercury_geometric_heliocentric_position():
    """Tests the geometric_heliocentric_position() method of Mercury class"""

    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Mercury.geometric_heliocentric_position(epoch)

    assert abs(round(lon.to_positive(), 4) - 287.4887) < TOL, \
        "ERROR: 1st geometric_heliocentric_position() test doesn't match"

    assert abs(round(lat, 4) - (-6.0086)) < TOL, \
        "ERROR: 2nd geometric_heliocentric_position() test doesn't match"

    assert abs(round(r, 5) - 0.45113) < TOL, \
        "ERROR: 3rd geometric_heliocentric_position() test doesn't match"


def test_mercury_orbital_elements_mean_equinox():
    """Tests the orbital_elements_mean_equinox() method of Mercury class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mercury.orbital_elements_mean_equinox(epoch)

    assert abs(round(l, 6) - 203.494701) < TOL, \
        "ERROR: 1st orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(a, 8) - 0.38709831) < TOL, \
        "ERROR: 2nd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(e, 7) - 0.2056451) < TOL, \
        "ERROR: 3rd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(i, 6) - 7.006171) < TOL, \
        "ERROR: 4th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(ome, 5) - 49.10765) < TOL, \
        "ERROR: 5th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(arg, 6) - 29.367732) < TOL, \
        "ERROR: 6th orbital_elements_mean_equinox() test doesn't match"


def test_mercury_orbital_elements_j2000():
    """Tests the orbital_elements_j2000() method of Mercury class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mercury.orbital_elements_j2000(epoch)

    assert abs(round(l, 6) - 202.579453) < TOL, \
        "ERROR: 1st orbital_elements_j2000() test doesn't match"

    assert abs(round(a, 8) - 0.38709831) < TOL, \
        "ERROR: 2nd orbital_elements_j2000() test doesn't match"

    assert abs(round(e, 7) - 0.2056451) < TOL, \
        "ERROR: 3rd orbital_elements_j2000() test doesn't match"

    assert abs(round(i, 6) - 7.001089) < TOL, \
        "ERROR: 4th orbital_elements_j2000() test doesn't match"

    assert abs(round(ome, 5) - 48.24873) < TOL, \
        "ERROR: 5th orbital_elements_j2000() test doesn't match"

    assert abs(round(arg, 6) - 29.311401) < TOL, \
        "ERROR: 6th orbital_elements_j2000() test doesn't match"
