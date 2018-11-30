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


def test_mars_orbital_elements_mean_equinox():
    """Tests the orbital_elements_mean_equinox() method of Mars class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mars.orbital_elements_mean_equinox(epoch)

    assert abs(round(l, 6) - 288.855211) < TOL, \
        "ERROR: 1st orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(a, 8) - 1.52367934) < TOL, \
        "ERROR: 2nd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(e, 7) - 0.0934599) < TOL, \
        "ERROR: 3rd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(i, 6) - 1.849338) < TOL, \
        "ERROR: 4th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(ome, 5) - 50.06365) < TOL, \
        "ERROR: 5th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(arg, 6) - 287.202108) < TOL, \
        "ERROR: 6th orbital_elements_mean_equinox() test doesn't match"


def test_mars_orbital_elements_j2000():
    """Tests the orbital_elements_j2000() method of Mars class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mars.orbital_elements_j2000(epoch)

    assert abs(round(l, 6) - 287.94027) < TOL, \
        "ERROR: 1st orbital_elements_j2000() test doesn't match"

    assert abs(round(a, 8) - 1.52367934) < TOL, \
        "ERROR: 2nd orbital_elements_j2000() test doesn't match"

    assert abs(round(e, 7) - 0.0934599) < TOL, \
        "ERROR: 3rd orbital_elements_j2000() test doesn't match"

    assert abs(round(i, 6) - 1.844381) < TOL, \
        "ERROR: 4th orbital_elements_j2000() test doesn't match"

    assert abs(round(ome, 5) - 49.36464) < TOL, \
        "ERROR: 5th orbital_elements_j2000() test doesn't match"

    assert abs(round(arg, 6) - 286.98617) < TOL, \
        "ERROR: 6th orbital_elements_j2000() test doesn't match"


def test_mars_conjunction():
    """Tests the conjunction() method of Mars class"""

    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mars.conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st conjunction() test doesn't match"

    assert abs(round(m, 0) - 12) < TOL, \
        "ERROR: 2nd conjunction() test doesn't match"

    assert abs(round(d, 4) - 27.0898) < TOL, \
        "ERROR: 3rd conjunction() test doesn't match"


def test_mars_superior_conjunction():
    """Tests the opposition() method of Mars class"""

    epoch = Epoch(2729, 10, 1.0)
    oppo = Mars.opposition(epoch)
    y, m, d = oppo.get_date()

    assert abs(round(y, 0) - 2729) < TOL, \
        "ERROR: 1st opposition() test doesn't match"

    assert abs(round(m, 0) - 9) < TOL, \
        "ERROR: 2nd opposition() test doesn't match"

    assert abs(round(d, 4) - 9.1412) < TOL, \
        "ERROR: 3rd opposition() test doesn't match"


def test_mars_station_longitude_1():
    """Tests the station_longitude_1() method of Mars class"""

    epoch = Epoch(1997, 3, 1.0)
    sta1 = Mars.station_longitude_1(epoch)
    y, m, d = sta1.get_date()

    assert abs(round(y, 0) - 1997) < TOL, \
        "ERROR: 1st station_longitude_1() test doesn't match"

    assert abs(round(m, 0) - 2) < TOL, \
        "ERROR: 2nd station_longitude_1() test doesn't match"

    assert abs(round(d, 4) - 6.033) < TOL, \
        "ERROR: 3rd station_longitude_1() test doesn't match"


def test_mars_station_longitude_2():
    """Tests the station_longitude_2() method of Mars class"""

    epoch = Epoch(1997, 3, 1.0)
    sta2 = Mars.station_longitude_2(epoch)
    y, m, d = sta2.get_date()

    assert abs(round(y, 0) - 1997) < TOL, \
        "ERROR: 1st station_longitude_2() test doesn't match"

    assert abs(round(m, 0) - 4) < TOL, \
        "ERROR: 2nd station_longitude_2() test doesn't match"

    assert abs(round(d, 4) - 27.7553) < TOL, \
        "ERROR: 3rd station_longitude_2() test doesn't match"
