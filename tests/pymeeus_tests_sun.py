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
from pymeeus.Sun import sun_apparent_rightascension_declination_coarse, \
        sun_true_longitude_coarse, sun_apparent_longitude_coarse
from pymeeus.Epoch import Epoch


# Coordinates module

def test_sun_true_longitude_coarse():
    """Tests the sun_true_longitude_coarse() method of Coordinates module"""

    epoch = Epoch(1992, 10, 13)
    true_lon, r = sun_true_longitude_coarse(epoch)

    assert true_lon.dms_str(n_dec=0) == "199d 54' 36.0''", \
        "ERROR: 1st sun_true_longitude_coarse() test, 'true_lon' doesn't match"

    assert abs(round(r, 5) - 0.99766) < TOL, \
        "ERROR: 2nd sun_true_longitude_coarse() test, 'r' value doesn't match"


def test_sun_apparent_longitude_coarse():
    """Tests sun_apparent_longitude_coarse() method of Coordinates module"""

    epoch = Epoch(1992, 10, 13)
    alon, r = sun_apparent_longitude_coarse(epoch)

    assert alon.dms_str(n_dec=0) == "199d 54' 32.0''", \
        "ERROR: 1st sun_apparent_longitude_coarse() test, 'alon' doesn't match"


def test_sun_apparent_rightascension_declination_coarse():
    """Tests sun_apparent_rightascension_declination_coarse() method of
    Coordinates module"""

    epoch = Epoch(1992, 10, 13)
    ra, delta, r = sun_apparent_rightascension_declination_coarse(epoch)

    assert ra.ra_str(n_dec=1) == "13h 13' 31.4''", \
        "ERROR: 1st sun_rightascension_declination_coarse() test doesn't match"

    assert delta.dms_str(n_dec=0) == "-7d 47' 6.0''", \
        "ERROR: 2nd sun_rightascension_declination_coarse() test doesn't match"
