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
from pymeeus.Mars import Mars
from pymeeus.Epoch import Epoch


# Mars class

def test_mars_geometric_heliocentric_position():
    """Tests the geometric_heliocentric_position() method of Mars class"""

    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Mars.geometric_heliocentric_position(epoch)

    assert abs(round(lon.to_positive(), 4) - 2.0015) < TOL, \
        "ERROR: 1st geometric_heliocentric_position() test doesn't match"

    assert abs(round(lat, 4) - (-1.3683)) < TOL, \
        "ERROR: 2nd geometric_heliocentric_position() test doesn't match"

    assert abs(round(r, 5) - 1.39306) < TOL, \
        "ERROR: 3rd geometric_heliocentric_position() test doesn't match"
