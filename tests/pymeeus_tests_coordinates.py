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
from pymeeus.Coordinates import mean_obliquity, true_obliquity, \
        nutation_longitude, nutation_obliquity, precession_equatorial, \
        precession_ecliptical, motion_in_space, equatorial2ecliptical, \
        ecliptical2equatorial, equatorial2horizontal, horizontal2equatorial, \
        equatorial2galactic, galactic2equatorial, ecliptic_horizon, \
        diurnal_path_horizon
from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch, JDE2000
from math import cos


# Coordinates module

def test_coordinates_mean_obliquity():
    """Tests the mean_obliquity() method of Coordinates module"""

    e0 = mean_obliquity(1987, 4, 10)
    a = e0.dms_tuple()
    assert abs(a[0] - 23.0) < TOL, \
        "ERROR: 1st mean_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 26.0) < TOL, \
        "ERROR: 2nd mean_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 27.407) < TOL, \
        "ERROR: 3rd mean_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th mean_obliquity() test, 'sign' value doesn't match"


def test_coordinates_true_obliquity():
    """Tests the true_obliquity() method of Coordinates module"""

    epsilon = true_obliquity(1987, 4, 10)
    a = epsilon.dms_tuple()
    assert abs(a[0] - 23.0) < TOL, \
        "ERROR: 1st true_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 26.0) < TOL, \
        "ERROR: 2nd true_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 36.849) < TOL, \
        "ERROR: 3rd true_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th true_obliquity() test, 'sign' value doesn't match"


def test_coordinates_nutation_longitude():
    """Tests the nutation_longitude() method of Coordinates module"""

    dpsi = nutation_longitude(1987, 4, 10)
    a = dpsi.dms_tuple()
    assert abs(a[0] - 0.0) < TOL, \
        "ERROR: 1st nutation_longitude() test, 'degrees' value doesn't match"

    assert abs(a[1] - 0.0) < TOL, \
        "ERROR: 2nd nutation_longitude() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 3.788) < TOL, \
        "ERROR: 3rd nutation_longitude() test, 'seconds value doesn't match"

    assert abs(a[3] - (-1.0)) < TOL, \
        "ERROR: 4th nutation_longitude() test, 'sign' value doesn't match"


def test_coordinates_nutation_obliquity():
    """Tests the nutation_obliquity() method of Coordinates module"""

    depsilon = nutation_obliquity(1987, 4, 10)
    a = depsilon.dms_tuple()
    assert abs(a[0] - 0.0) < TOL, \
        "ERROR: 1st nutation_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 0.0) < TOL, \
        "ERROR: 2nd nutation_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 9.443) < TOL, \
        "ERROR: 3rd nutation_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th nutation_obliquity() test, 'sign' value doesn't match"


def test_coordinates_precession_equatorial():
    """Tests the precession_equatorial() method of Coordinates module"""

    start_epoch = JDE2000
    final_epoch = Epoch(2028, 11, 13.19)
    alpha0 = Angle(2, 44, 11.986, ra=True)
    delta0 = Angle(49, 13, 42.48)
    pm_ra = Angle(0, 0, 0.03425, ra=True)
    pm_dec = Angle(0, 0, -0.0895)

    alpha, delta = precession_equatorial(start_epoch, final_epoch, alpha0,
                                         delta0, pm_ra, pm_dec)

    assert alpha.ra_str(False, 3) == "2:46:11.331", \
        "ERROR: 1st precession_equatorial test, right ascension doesn't match"

    assert delta.dms_str(False, 2) == "49:20:54.54", \
        "ERROR: 2nd precession_equatorial() test, 'declination' doesn't match"


def test_coordinates_precession_ecliptical():
    """Tests the precession_ecliptical() method of Coordinates module"""

    start_epoch = JDE2000
    final_epoch = Epoch(-214, 6, 30.0)
    lon0 = Angle(149.48194)
    lat0 = Angle(1.76549)

    lon, lat = precession_ecliptical(start_epoch, final_epoch, lon0, lat0)

    assert abs(round(lon(), 3) - 118.704) < TOL, \
        "ERROR: 1st precession_ecliptical() test, 'longitude' doesn't match"

    assert abs(round(lat(), 3) - 1.615) < TOL, \
        "ERROR: 2nd precession_ecliptical() test, 'latitude' doesn't match"


def test_coordinates_motion_in_space():
    """Tests the motion_in_space() method of Coordinates module"""

    ra = Angle(6, 45, 8.871, ra=True)
    dec = Angle(-16.716108)
    pm_ra = Angle(0, 0, -0.03847, ra=True)
    pm_dec = Angle(0, 0, -1.2053)
    dist = 2.64
    vel = -7.6

    alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -2000.0)

    assert alpha.ra_str(False, 2) == "6:46:25.09", \
        "ERROR: 1st motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-16:3:0.8", \
        "ERROR: 2nd motion_in_space() test, 'declination' doesn't match"

    alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -3000.0)

    assert alpha.ra_str(False, 2) == "6:47:2.67", \
        "ERROR: 3rd motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-15:43:12.3", \
        "ERROR: 4th motion_in_space() test, 'declination' doesn't match"

    alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -12000.0)

    assert alpha.ra_str(False, 2) == "6:52:25.72", \
        "ERROR: 5th motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-12:50:6.7", \
        "ERROR: 6th motion_in_space() test, 'declination' doesn't match"


def test_coordinates_equatorial2ecliptical():
    """Tests the equatorial2ecliptical() method of Coordinates module"""

    ra = Angle(7, 45, 18.946, ra=True)
    dec = Angle(28, 1, 34.26)
    epsilon = Angle(23.4392911)
    lon, lat = equatorial2ecliptical(ra, dec, epsilon)

    assert abs(round(lon(), 5) - 113.21563) < TOL, \
        "ERROR: 1st equatorial2ecliptical() test, 'longitude' doesn't match"

    assert abs(round(lat(), 5) - 6.68417) < TOL, \
        "ERROR: 2nd equatorial2ecliptical() test, 'latitude' doesn't match"


def test_coordinates_ecliptical2equatorial():
    """Tests the ecliptical2equatorial() method of Coordinates module"""

    lon = Angle(113.21563)
    lat = Angle(6.68417)
    epsilon = Angle(23.4392911)
    ra, dec = ecliptical2equatorial(lon, lat, epsilon)

    assert ra.ra_str(n_dec=3) == "7h 45' 18.946''", \
        "ERROR: 1st ecliptical2equatorial() test, 'ra' doesn't match"

    assert dec.dms_str(n_dec=2) == "28d 1' 34.26''", \
        "ERROR: 2nd ecliptical2equatorial() test, 'declination' doesn't match"


def test_coordinates_equatorial2horizontal():
    """Tests the equatorial2horizontal() method of Coordinates module"""

    lon = Angle(77, 3, 56)
    lat = Angle(38, 55, 17)
    ra = Angle(23, 9, 16.641, ra=True)
    dec = Angle(-6, 43, 11.61)
    theta0 = Angle(8, 34, 57.0896, ra=True)
    eps = Angle(23, 26, 36.87)
    delta = Angle(0, 0, ((-3.868*cos(eps.rad()))/15.0), ra=True)
    theta0 += delta
    h = theta0 - lon - ra
    azi, ele = equatorial2horizontal(h, dec, lat)

    assert abs(round(azi, 3) - 68.034) < TOL, \
        "ERROR: 1st equatorial2horizontal() test, 'azimuth' doesn't match"

    assert abs(round(ele, 3) - 15.125) < TOL, \
        "ERROR: 2nd equatorial2horizontal() test, 'elevation' doesn't match"


def test_coordinates_horizontal2equatorial():
    """Tests the horizontal2equatorial() method of Coordinates module"""

    azi = Angle(68.0337)
    ele = Angle(15.1249)
    lat = Angle(38, 55, 17)
    h, dec = horizontal2equatorial(azi, ele, lat)

    assert abs(round(h, 4) - 64.3521) < TOL, \
        "ERROR: 1st horizontal2equatorial() test, 'hour angle' doesn't match"

    assert dec.dms_str(n_dec=0) == "-6d 43' 12.0''", \
        "ERROR: 2nd horizontal2equatorial() test, 'declination' match"


def test_coordinates_equatorial2galactic():
    """Tests the equatorial2galactic() method of Coordinates module"""

    ra = Angle(17, 48, 59.74, ra=True)
    dec = Angle(-14, 43, 8.2)
    lon, lat = equatorial2galactic(ra, dec)

    assert abs(round(lon, 4) - 12.9593) < TOL, \
        "ERROR: 1st equatorial2galactic() test, 'longitude' doesn't match"

    assert abs(round(lat, 4) - 6.0463) < TOL, \
        "ERROR: 2nd equatorial2galactic() test, 'latitude' doesn't match"


def test_coordinates_galactic2equatorial():
    """Tests the galactic2equatorial() method of Coordinates module"""

    lon = Angle(12.9593)
    lat = Angle(6.0463)
    ra, dec = galactic2equatorial(lon, lat)

    assert ra.ra_str(n_dec=1) == "17h 48' 59.7''", \
        "ERROR: 1st galactic2equatorial() test, 'ra' doesn't match"

    assert dec.dms_str(n_dec=0) == "-14d 43' 8.0''", \
        "ERROR: 2nd galactic2equatorial() test, 'declination' doesn't match"


def test_coordinates_ecliptic_horizon():
    """Tests the ecliptic_horizon() method of Coordinates module"""

    sidereal_time = Angle(5.0, ra=True)
    lat = Angle(51.0)
    epsilon = Angle(23.44)
    lon1, lon2, i = ecliptic_horizon(sidereal_time, lat, epsilon)

    assert lon1.dms_str(n_dec=1) == "169d 21' 29.9''", \
        "ERROR: 1st ecliptic_horizon() test, 'lon1' doesn't match"

    assert lon2.dms_str(n_dec=1) == "349d 21' 29.9''", \
        "ERROR: 2nd ecliptic_horizon() test, 'lon2' doesn't match"

    assert abs(round(i, 0) - 62.0) < TOL, \
        "ERROR: 3rd ecliptic_horizon() test, 'i' angle doesn't match"


def test_coordinates_diurnal_path_horizon():
    """Tests the diurnal_path_horizon() method of Coordinates module"""

    dec = Angle(23.44)
    lat = Angle(40.0)
    j = diurnal_path_horizon(dec, lat)

    assert j.dms_str(n_dec=1) == "45d 31' 28.4''", \
        "ERROR: 1st diurnal_path_horizon() test, 'j' angle doesn't match"
