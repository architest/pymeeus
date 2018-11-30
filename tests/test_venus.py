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
from pymeeus.Venus import Venus
from pymeeus.Epoch import Epoch


# Venus class

def test_venus_geometric_heliocentric_position():
    """Tests the geometric_heliocentric_position() method of Venus class"""

    epoch = Epoch(1992, 12, 20.0)
    lon, lat, r = Venus.geometric_heliocentric_position(epoch, tofk5=False)

    assert abs(round(lon.to_positive(), 5) - 26.11412) < TOL, \
        "ERROR: 1st geometric_heliocentric_position() test doesn't match"

    assert abs(round(lat, 4) - (-2.6206)) < TOL, \
        "ERROR: 2nd geometric_heliocentric_position() test doesn't match"

    assert abs(round(r, 6) - 0.724602) < TOL, \
        "ERROR: 3rd geometric_heliocentric_position() test doesn't match"


def test_venus_orbital_elements_mean_equinox():
    """Tests the orbital_elements_mean_equinox() method of Venus class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Venus.orbital_elements_mean_equinox(epoch)

    assert abs(round(l, 6) - 338.646306) < TOL, \
        "ERROR: 1st orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(a, 8) - 0.72332982) < TOL, \
        "ERROR: 2nd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(e, 7) - 0.0067407) < TOL, \
        "ERROR: 3rd orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(i, 6) - 3.395319) < TOL, \
        "ERROR: 4th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(ome, 5) - 77.27012) < TOL, \
        "ERROR: 5th orbital_elements_mean_equinox() test doesn't match"

    assert abs(round(arg, 6) - 55.211257) < TOL, \
        "ERROR: 6th orbital_elements_mean_equinox() test doesn't match"


def test_venus_orbital_elements_j2000():
    """Tests the orbital_elements_j2000() method of Venus class"""

    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Venus.orbital_elements_j2000(epoch)

    assert abs(round(l, 6) - 337.731227) < TOL, \
        "ERROR: 1st orbital_elements_j2000() test doesn't match"

    assert abs(round(a, 8) - 0.72332982) < TOL, \
        "ERROR: 2nd orbital_elements_j2000() test doesn't match"

    assert abs(round(e, 7) - 0.0067407) < TOL, \
        "ERROR: 3rd orbital_elements_j2000() test doesn't match"

    assert abs(round(i, 6) - 3.394087) < TOL, \
        "ERROR: 4th orbital_elements_j2000() test doesn't match"

    assert abs(round(ome, 5) - 76.49782) < TOL, \
        "ERROR: 5th orbital_elements_j2000() test doesn't match"

    assert abs(round(arg, 6) - 55.068476) < TOL, \
        "ERROR: 6th orbital_elements_j2000() test doesn't match"


def test_venus_inferior_conjunction():
    """Tests the inferior_conjunction() method of Venus class"""

    epoch = Epoch(1882, 12, 1.0)
    conjunction = Venus.inferior_conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1882) < TOL, \
        "ERROR: 1st inferior_conjunction() test doesn't match"

    assert abs(round(m, 0) - 12) < TOL, \
        "ERROR: 2nd inferior_conjunction() test doesn't match"

    assert abs(round(d, 1) - 6.7) < TOL, \
        "ERROR: 3rd inferior_conjunction() test doesn't match"


def test_venus_superior_conjunction():
    """Tests the superior_conjunction() method of Venus class"""

    epoch = Epoch(1993, 10, 1.0)
    conjunction = Venus.superior_conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1994) < TOL, \
        "ERROR: 1st superior_conjunction() test doesn't match"

    assert abs(round(m, 0) - 1) < TOL, \
        "ERROR: 2nd superior_conjunction() test doesn't match"

    assert abs(round(d, 2) - 17.05) < TOL, \
        "ERROR: 3rd superior_conjunction() test doesn't match"


def test_venus_western_elongation():
    """Tests the western_elongation() method of Venus class"""

    epoch = Epoch(2019, 1, 1.0)
    time, elongation = Venus.western_elongation(epoch)
    y, m, d = time.get_date()

    assert abs(round(y, 0) - 2019) < TOL, \
        "ERROR: 1st western_elongation() test doesn't match"

    assert abs(round(m, 0) - 1) < TOL, \
        "ERROR: 2nd western_elongation() test doesn't match"

    assert abs(round(d, 4) - 6.1895) < TOL, \
        "ERROR: 3rd western_elongation() test doesn't match"

    assert abs(round(elongation, 4) - 46.9571) < TOL, \
        "ERROR: 4th western_elongation() test doesn't match"


def test_venus_eastern_elongation():
    """Tests the eastern_elongation() method of Venus class"""

    epoch = Epoch(2019, 10, 1.0)
    time, elongation = Venus.eastern_elongation(epoch)
    y, m, d = time.get_date()

    assert abs(round(y, 0) - 2020) < TOL, \
        "ERROR: 1st eastern_elongation() test doesn't match"

    assert abs(round(m, 0) - 3) < TOL, \
        "ERROR: 2nd eastern_elongation() test doesn't match"

    assert abs(round(d, 4) - 24.9179) < TOL, \
        "ERROR: 3rd eastern_elongation() test doesn't match"

    assert abs(round(elongation, 4) - 46.078) < TOL, \
        "ERROR: 3rd eastern_elongation() test doesn't match"


def test_venus_station_longitude_1():
    """Tests the station_longitude_1() method of Venus class"""

    epoch = Epoch(2018, 12, 1.0)
    sta1 = Venus.station_longitude_1(epoch)
    y, m, d = sta1.get_date()

    assert abs(round(y, 0) - 2018) < TOL, \
        "ERROR: 1st station_longitude_1() test doesn't match"

    assert abs(round(m, 0) - 10) < TOL, \
        "ERROR: 2nd station_longitude_1() test doesn't match"

    assert abs(round(d, 4) - 5.7908) < TOL, \
        "ERROR: 3rd station_longitude_1() test doesn't match"


def test_venus_station_longitude_2():
    """Tests the station_longitude_2() method of Venus class"""

    epoch = Epoch(2018, 12, 1.0)
    sta2 = Venus.station_longitude_2(epoch)
    y, m, d = sta2.get_date()

    assert abs(round(y, 0) - 2018) < TOL, \
        "ERROR: 1st station_longitude_2() test doesn't match"

    assert abs(round(m, 0) - 11) < TOL, \
        "ERROR: 2nd station_longitude_2() test doesn't match"

    assert abs(round(d, 4) - 16.439) < TOL, \
        "ERROR: 3rd station_longitude_2() test doesn't match"
