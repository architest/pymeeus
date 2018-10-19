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
from pymeeus.Sun import Sun
from pymeeus.Epoch import Epoch


# Sun module

def test_sun_true_longitude_coarse():
    """Tests the true_longitude_coarse() method of Sun class"""

    epoch = Epoch(1992, 10, 13)
    true_lon, r = Sun.true_longitude_coarse(epoch)

    assert true_lon.dms_str(n_dec=0) == "199d 54' 36.0''", \
        "ERROR: 1st true_longitude_coarse() test, 'true_lon' doesn't match"

    assert abs(round(r, 5) - 0.99766) < TOL, \
        "ERROR: 2nd true_longitude_coarse() test, 'r' value doesn't match"


def test_sun_apparent_longitude_coarse():
    """Tests apparent_longitude_coarse() method of Sun class"""

    epoch = Epoch(1992, 10, 13)
    alon, r = Sun.apparent_longitude_coarse(epoch)

    assert alon.dms_str(n_dec=0) == "199d 54' 32.0''", \
        "ERROR: 1st apparent_longitude_coarse() test, 'alon' doesn't match"


def test_sun_apparent_rightascension_declination_coarse():
    """Tests apparent_rightascension_declination_coarse() method of Sun
    class"""

    epoch = Epoch(1992, 10, 13)
    ra, delta, r = Sun.apparent_rightascension_declination_coarse(epoch)

    assert ra.ra_str(n_dec=1) == "13h 13' 31.4''", \
        "ERROR: 1st rightascension_declination_coarse() test doesn't match"

    assert delta.dms_str(n_dec=0) == "-7d 47' 6.0''", \
        "ERROR: 2nd rightascension_declination_coarse() test doesn't match"


def test_sun_geometric_geocentric_position():
    """Tests the geometric_geocentric_position() method of Earth class"""

    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = Sun.geometric_geocentric_position(epoch, toFK5=False)

    assert abs(round(lon.to_positive(), 6) - 199.907297) < TOL, \
        "ERROR: 1st geometric_geocentric_position() test, 'lon' doesn't match"

    assert lat.dms_str(n_dec=3) == "0.744''", \
        "ERROR: 2nd geometric_geocentric_position() test, 'lat' doesn't match"

    assert abs(round(r, 8) - 0.99760852) < TOL, \
        "ERROR: 3rd geometric_geocentric_position() test, 'r' doesn't match"


def test_sun_apparent_geocentric_position():
    """Tests the apparent_geocentric_position() method of Earth class"""

    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = Sun.apparent_geocentric_position(epoch)

    assert lon.to_positive().dms_str(n_dec=3) == "199d 54' 21.548''", \
        "ERROR: 1st apparent_geocentric_position() test, 'lon' doesn't match"

    assert lat.dms_str(n_dec=3) == "0.721''", \
        "ERROR: 2nd apparent_geocentric_position() test, 'lat' doesn't match"

    assert abs(round(r, 8) - 0.99760852) < TOL, \
        "ERROR: 3rd apparent_geocentric_position() test, 'r' doesn't match"
