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


def test_mercury_inferior_conjunction():
    """Tests the inferior_conjunction() method of Mercury class"""

    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury.inferior_conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st inferior_conjunction() test doesn't match"

    assert abs(round(m, 0) - 11) < TOL, \
        "ERROR: 2nd inferior_conjunction() test doesn't match"

    assert abs(round(d, 4) - 6.1449) < TOL, \
        "ERROR: 3rd inferior_conjunction() test doesn't match"

    epoch = Epoch(1631, 10, 1.0)
    conjunction = Mercury.inferior_conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1631) < TOL, \
        "ERROR: 4th inferior_conjunction() test doesn't match"

    assert abs(round(m, 0) - 11) < TOL, \
        "ERROR: 5th inferior_conjunction() test doesn't match"

    assert abs(round(d, 3) - 7.306) < TOL, \
        "ERROR: 6th inferior_conjunction() test doesn't match"


def test_mercury_superior_conjunction():
    """Tests the superior_conjunction() method of Mercury class"""

    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury.superior_conjunction(epoch)
    y, m, d = conjunction.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st superior_conjunction() test doesn't match"

    assert abs(round(m, 0) - 8) < TOL, \
        "ERROR: 2nd superior_conjunction() test doesn't match"

    assert abs(round(d, 4) - 29.3301) < TOL, \
        "ERROR: 3rd superior_conjunction() test doesn't match"


def test_mercury_western_elongation():
    """Tests the western_elongation() method of Mercury class"""

    epoch = Epoch(1993, 11, 1.0)
    time, elongation = Mercury.western_elongation(epoch)
    y, m, d = time.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st western_elongation() test doesn't match"

    assert abs(round(m, 0) - 11) < TOL, \
        "ERROR: 2nd western_elongation() test doesn't match"

    assert abs(round(d, 4) - 22.6386) < TOL, \
        "ERROR: 3rd western_elongation() test doesn't match"

    assert abs(round(elongation, 4) - 19.7506) < TOL, \
        "ERROR: 4th western_elongation() test doesn't match"


def test_mercury_eastern_elongation():
    """Tests the eastern_elongation() method of Mercury class"""

    epoch = Epoch(1990, 8, 1.0)
    time, elongation = Mercury.eastern_elongation(epoch)
    y, m, d = time.get_date()

    assert abs(round(y, 0) - 1990) < TOL, \
        "ERROR: 1st eastern_elongation() test doesn't match"

    assert abs(round(m, 0) - 8) < TOL, \
        "ERROR: 2nd eastern_elongation() test doesn't match"

    assert abs(round(d, 4) - 11.8514) < TOL, \
        "ERROR: 3rd eastern_elongation() test doesn't match"

    assert abs(round(elongation, 4) - 27.4201) < TOL, \
        "ERROR: 4th eastern_elongation() test doesn't match"


def test_mercury_station_longitude_1():
    """Tests the station_longitude_1() method of Mercury class"""

    epoch = Epoch(1993, 10, 1.0)
    sta1 = Mercury.station_longitude_1(epoch)
    y, m, d = sta1.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st station_longitude_1() test doesn't match"

    assert abs(round(m, 0) - 10) < TOL, \
        "ERROR: 2nd station_longitude_1() test doesn't match"

    assert abs(round(d, 4) - 25.9358) < TOL, \
        "ERROR: 3rd station_longitude_1() test doesn't match"


def test_mercury_station_longitude_2():
    """Tests the station_longitude_2() method of Mercury class"""

    epoch = Epoch(1993, 10, 1.0)
    sta2 = Mercury.station_longitude_2(epoch)
    y, m, d = sta2.get_date()

    assert abs(round(y, 0) - 1993) < TOL, \
        "ERROR: 1st station_longitude_2() test doesn't match"

    assert abs(round(m, 0) - 11) < TOL, \
        "ERROR: 2nd station_longitude_2() test doesn't match"

    assert abs(round(d, 4) - 15.0724) < TOL, \
        "ERROR: 3rd station_longitude_2() test doesn't match"
