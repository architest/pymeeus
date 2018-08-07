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
from pymeeus.Earth import Earth, IAU76
from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch, JDE2000


# Earth class

def test_earth_constructor():
    """Tests the constructor of Earth class"""

    a = Earth()
    assert abs(a._ellip._a - 6378137.0) < TOL, \
        "ERROR: 1st constructor test, 'a' value doesn't match"

    assert abs(a._ellip._f - (1.0/298.257223563)) < TOL, \
        "ERROR: 2nd constructor test, 'f' value doesn't match"

    assert abs(a._ellip._omega - 7292115e-11) < TOL, \
        "ERROR: 3rd constructor test, 'omega' value doesn't match"

    a = Earth(ellipsoid=IAU76)
    assert abs(a._ellip._a - 6378140.0) < TOL, \
        "ERROR: 4th constructor test, 'a' value doesn't match"

    assert abs(a._ellip._f - (1.0/298.257)) < TOL, \
        "ERROR: 5th constructor test, 'f' value doesn't match"

    assert abs(a._ellip._omega - 7.292114992e-5) < TOL, \
        "ERROR: 6th constructor test, 'omega' value doesn't match"


def test_earth_rho():
    """Tests the rho() method of Earth class"""

    a = Earth()
    assert abs(a.rho(0.0) - 1.0) < TOL, \
        "ERROR: 1st rho() test, output doesn't match"


def test_earth_rho_sinphi():
    """Tests the rho_sinphi() method of Earth class"""

    lat = Angle(33, 21, 22.0)
    e = Earth(ellipsoid=IAU76)
    assert abs(round(e.rho_sinphi(lat, 1706), 6) - 0.546861) < TOL, \
        "ERROR: 1st rho_sinphi() test, output doesn't match"


def test_earth_rho_cosphi():
    """Tests the rho_cosphi() method of Earth class"""

    lat = Angle(33, 21, 22.0)
    e = Earth(ellipsoid=IAU76)
    assert abs(round(e.rho_cosphi(lat, 1706), 6) - 0.836339) < TOL, \
        "ERROR: 1st rho_cosphi() test, output doesn't match"


def test_earth_rp():
    """Tests the rp() method of Earth class"""

    e = Earth(ellipsoid=IAU76)
    assert abs(round(e.rp(42.0), 1) - 4747001.2) < TOL, \
        "ERROR: 1st rp() test, output doesn't match"


def test_earth_rm():
    """Tests the rm() method of Earth class"""

    e = Earth(ellipsoid=IAU76)
    assert abs(round(e.rm(42.0), 1) - 6364033.3) < TOL, \
        "ERROR: 1st rm() test, output doesn't match"


def test_earth_linear_velocity():
    """Tests the linear_velocity() method of Earth class"""

    e = Earth(ellipsoid=IAU76)
    assert abs(round(e.linear_velocity(42.0), 2) - 346.16) < TOL, \
        "ERROR: 1st linear_velocity() test, output doesn't match"


def test_earth_distance():
    """Tests the distance() method of Earth class"""

    e = Earth(ellipsoid=IAU76)
    lon1 = Angle(-2, 20, 14.0)
    lat1 = Angle(48, 50, 11.0)
    lon2 = Angle(77, 3, 56.0)
    lat2 = Angle(38, 55, 17.0)
    dist, error = e.distance(lon1, lat1, lon2, lat2)

    assert abs(round(dist, 0) - 6181628.0) < TOL, \
        "ERROR: 1st distance() test, output doesn't match"

    assert abs(round(error, 0) - 69.0) < TOL, \
        "ERROR: 2nd distance() test, output doesn't match"

    lon1 = Angle(-2.09)
    lat1 = Angle(41.3)
    lon2 = Angle(73.99)
    lat2 = Angle(40.75)
    dist, error = e.distance(lon1, lat1, lon2, lat2)

    assert abs(round(dist, 0) - 6176760.0) < TOL, \
        "ERROR: 3rd distance() test, output doesn't match"

    assert abs(round(error, 0) - 69.0) < TOL, \
        "ERROR: 4th distance() test, output doesn't match"


def test_earth_mean_obliquity():
    """Tests the mean_obliquity() method of Earth class"""

    e0 = Earth.mean_obliquity(1987, 4, 10)
    a = e0.dms_tuple()
    assert abs(a[0] - 23.0) < TOL, \
        "ERROR: 1st mean_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 26.0) < TOL, \
        "ERROR: 2nd mean_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 27.407) < TOL, \
        "ERROR: 3rd mean_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th mean_obliquity() test, 'sign' value doesn't match"


def test_earth_true_obliquity():
    """Tests the true_obliquity() method of Earth class"""

    epsilon = Earth.true_obliquity(1987, 4, 10)
    a = epsilon.dms_tuple()
    assert abs(a[0] - 23.0) < TOL, \
        "ERROR: 1st true_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 26.0) < TOL, \
        "ERROR: 2nd true_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 36.849) < TOL, \
        "ERROR: 3rd true_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th true_obliquity() test, 'sign' value doesn't match"


def test_earth_nutation_longitude():
    """Tests the nutation_longitude() method of Earth class"""

    dpsi = Earth.nutation_longitude(1987, 4, 10)
    a = dpsi.dms_tuple()
    assert abs(a[0] - 0.0) < TOL, \
        "ERROR: 1st nutation_longitude() test, 'degrees' value doesn't match"

    assert abs(a[1] - 0.0) < TOL, \
        "ERROR: 2nd nutation_longitude() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 3.788) < TOL, \
        "ERROR: 3rd nutation_longitude() test, 'seconds value doesn't match"

    assert abs(a[3] - (-1.0)) < TOL, \
        "ERROR: 4th nutation_longitude() test, 'sign' value doesn't match"


def test_earth_nutation_obliquity():
    """Tests the nutation_obliquity() method of Earth class"""

    depsilon = Earth.nutation_obliquity(1987, 4, 10)
    a = depsilon.dms_tuple()
    assert abs(a[0] - 0.0) < TOL, \
        "ERROR: 1st nutation_obliquity() test, 'degrees' value doesn't match"

    assert abs(a[1] - 0.0) < TOL, \
        "ERROR: 2nd nutation_obliquity() test, 'minutes' value doesn't match"

    assert abs(round(a[2], 3) - 9.443) < TOL, \
        "ERROR: 3rd nutation_obliquity() test, 'seconds value doesn't match"

    assert abs(a[3] - 1.0) < TOL, \
        "ERROR: 4th nutation_obliquity() test, 'sign' value doesn't match"


def test_earth_precession_equatorial():
    """Tests the precession_equatorial() method of Earth class"""

    start_epoch = JDE2000
    final_epoch = Epoch(2028, 11, 13.19)
    alpha0 = Angle(2, 44, 11.986, ra=True)
    delta0 = Angle(49, 13, 42.48)
    pm_ra = Angle(0, 0, 0.03425, ra=True)
    pm_dec = Angle(0, 0, -0.0895)

    alpha, delta = Earth.precession_equatorial(start_epoch, final_epoch,
                                               alpha0, delta0, pm_ra, pm_dec)

    assert alpha.ra_str(False, 3) == "2:46:11.331", \
        "ERROR: 1st precession_equatorial test, right ascension doesn't match"

    assert delta.dms_str(False, 2) == "49:20:54.54", \
        "ERROR: 2nd precession_equatorial() test, 'declination' doesn't match"


def test_earth_precession_ecliptical():
    """Tests the precession_ecliptical() method of Earth class"""

    start_epoch = JDE2000
    final_epoch = Epoch(-214, 6, 30.0)
    lon0 = Angle(149.48194)
    lat0 = Angle(1.76549)

    lon, lat = Earth.precession_ecliptical(start_epoch, final_epoch,
                                           lon0, lat0)

    assert abs(round(lon(), 3) - 118.704) < TOL, \
        "ERROR: 1st precession_ecliptical() test, 'longitude' doesn't match"

    assert abs(round(lat(), 3) - 1.615) < TOL, \
        "ERROR: 2nd precession_ecliptical() test, 'latitude' doesn't match"


def test_earth_motion_in_space():
    """Tests the motion_in_space() method of Earth class"""

    ra = Angle(6, 45, 8.871, ra=True)
    dec = Angle(-16.716108)
    pm_ra = Angle(0, 0, -0.03847, ra=True)
    pm_dec = Angle(0, 0, -1.2053)
    dist = 2.64
    vel = -7.6

    alpha, delta = Earth.motion_in_space(ra, dec, dist, vel,
                                         pm_ra, pm_dec, -2000.0)

    assert alpha.ra_str(False, 2) == "6:46:25.09", \
        "ERROR: 1st motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-16:3:0.8", \
        "ERROR: 2nd motion_in_space() test, 'declination' doesn't match"

    alpha, delta = Earth.motion_in_space(ra, dec, dist, vel,
                                         pm_ra, pm_dec, -3000.0)

    assert alpha.ra_str(False, 2) == "6:47:2.67", \
        "ERROR: 3rd motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-15:43:12.3", \
        "ERROR: 4th motion_in_space() test, 'declination' doesn't match"

    alpha, delta = Earth.motion_in_space(ra, dec, dist, vel,
                                         pm_ra, pm_dec, -12000.0)

    assert alpha.ra_str(False, 2) == "6:52:25.72", \
        "ERROR: 5th motion_in_space() test, 'right ascension' doesn't match"

    assert delta.dms_str(False, 1) == "-12:50:6.7", \
        "ERROR: 6th motion_in_space() test, 'declination' doesn't match"
