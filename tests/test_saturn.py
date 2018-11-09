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
from pymeeus.Saturn import Saturn
from pymeeus.Epoch import Epoch


# Saturn class

def test_saturn_geometric_heliocentric_position():
    """Tests the geometric_heliocentric_position() method of Saturn class"""

    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Saturn.geometric_heliocentric_position(epoch)

    assert abs(round(lon.to_positive(), 4) - 279.5108) < TOL, \
        "ERROR: 1st geometric_heliocentric_position() test doesn't match"

    assert abs(round(lat, 4) - 0.6141) < TOL, \
        "ERROR: 2nd geometric_heliocentric_position() test doesn't match"

    assert abs(round(r, 5) - 10.06266) < TOL, \
        "ERROR: 3rd geometric_heliocentric_position() test doesn't match"


def test_saturn_orbital_elements_mean_equinox():
    """Tests the orbital_elements_mean_equinox() method of Saturn class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Saturn.orbital_elements_mean_equinox(epoch)

    assert abs(round(l, 6) - 131.196871) < TOL, \
        "ERROR: 1st orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(a, 8) - 9.55490779) < TOL, \
        "ERROR: 2nd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(e, 7) - 0.0553209) < TOL, \
        "ERROR: 3rd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(i, 6) - 2.486426) < TOL, \
        "ERROR: 4th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(ome, 5) - 114.23974) < TOL, \
        "ERROR: 5th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(arg, 6) - (-19.896331)) < TOL, \
        "ERROR: 6th orbital_elements_mean_equinox() test doesn't match"


def test_saturn_orbital_elements_j2000():
    """Tests the orbital_elements_j2000() method of Saturn class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Saturn.orbital_elements_j2000(epoch)

    assert abs(round(l, 6) - 130.28188) < TOL, \
        "ERROR: 1st orbital_elements_j2000() test doesn't match"

    assert abs(round(a, 8) - 9.55490779) < TOL, \
        "ERROR: 2nd orbital_elements_j2000() test doesn't match"

    assert abs(round(e, 7) - 0.0553209) < TOL, \
        "ERROR: 3rd orbital_elements_j2000() test doesn't match"

    assert abs(round(i, 6) - 2.490529) < TOL, \
        "ERROR: 4th orbital_elements_j2000() test doesn't match"

    assert abs(round(ome, 5) - 113.49736) < TOL, \
        "ERROR: 5th orbital_elements_j2000() test doesn't match"

    assert abs(round(arg, 6) - (-20.068943)) < TOL, \
        "ERROR: 6th orbital_elements_j2000() test doesn't match"
