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
from pymeeus.Uranus import Uranus
from pymeeus.Epoch import Epoch


# Uranus class

def test_uranus_geometric_heliocentric_position():
    """Tests the geometric_heliocentric_position() method of Uranus class"""

    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Uranus.geometric_heliocentric_position(epoch)

    assert abs(round(lon.to_positive(), 4) - 30.5888) < TOL, \
        "ERROR: 1st geometric_heliocentric_position() test doesn't match"

    assert abs(round(lat, 4) - (-0.5315)) < TOL, \
        "ERROR: 2nd geometric_heliocentric_position() test doesn't match"

    assert abs(round(r, 5) - 19.86964) < TOL, \
        "ERROR: 3rd geometric_heliocentric_position() test doesn't match"


def test_uranus_orbital_elements_mean_equinox():
    """Tests the orbital_elements_mean_equinox() method of Uranus class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Uranus.orbital_elements_mean_equinox(epoch)

    assert abs(round(l, 6) - 235.517526) < TOL, \
        "ERROR: 1st orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(a, 8) - 19.21844604) < TOL, \
        "ERROR: 2nd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(e, 7) - 0.0463634) < TOL, \
        "ERROR: 3rd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(i, 6) - 0.77372) < TOL, \
        "ERROR: 4th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(ome, 5) - 74.34776) < TOL, \
        "ERROR: 5th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(arg, 6) - 99.630865) < TOL, \
        "ERROR: 6th orbital_elements_mean_equinox() test doesn't match"


def test_uranus_orbital_elements_j2000():
    """Tests the orbital_elements_j2000() method of Uranus class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Uranus.orbital_elements_j2000(epoch)

    assert abs(round(l, 6) - 234.602641) < TOL, \
        "ERROR: 1st orbital_elements_j2000() test doesn't match"

    assert abs(round(a, 8) - 19.21844604) < TOL, \
        "ERROR: 2nd orbital_elements_j2000() test doesn't match"

    assert abs(round(e, 7) - 0.0463634) < TOL, \
        "ERROR: 3rd orbital_elements_j2000() test doesn't match"

    assert abs(round(i, 6) - 0.772094) < TOL, \
        "ERROR: 4th orbital_elements_j2000() test doesn't match"

    assert abs(round(ome, 5) - 74.05468) < TOL, \
        "ERROR: 5th orbital_elements_j2000() test doesn't match"

    assert abs(round(arg, 6) - 99.009058) < TOL, \
        "ERROR: 6th orbital_elements_j2000() test doesn't match"
