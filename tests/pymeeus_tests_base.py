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


from math import sqrt, pi, degrees, radians, sin
import pymeeus.base


TOLERANCE = 1E-10


# Declare some objects to be used later
i_ra = pymeeus.base.Interpolation()
i_angles1 = pymeeus.base.Interpolation()
i_angles2 = pymeeus.base.Interpolation()
i_sine = pymeeus.base.Interpolation()
cf1 = pymeeus.base.CurveFitting()
cf2 = pymeeus.base.CurveFitting()
cf3 = pymeeus.base.CurveFitting()
cf4 = pymeeus.base.CurveFitting()


def setup():
    """This function is used to set up the environment for the tests"""
    # Set up a interpolation object which uses Right Ascension
    y0 = pymeeus.base.Angle(10, 18, 48.732, ra=True)
    y1 = pymeeus.base.Angle(10, 23, 22.835, ra=True)
    y2 = pymeeus.base.Angle(10, 27, 57.247, ra=True)
    y3 = pymeeus.base.Angle(10, 32, 31.983, ra=True)

    i_ra.set([8.0, 10.0, 12.0, 14.0], [y0, y1, y2, y3])

    # Set up a couple interpolation objects with Angles
    y0 = pymeeus.base.Angle(0, -28, 13.4)
    y1 = pymeeus.base.Angle(0, 6, 46.3)
    y2 = pymeeus.base.Angle(0, 38, 23.2)

    i_angles1.set([26.0, 27.0, 28.0], [y0, y1, y2])

    y0 = pymeeus.base.Angle(-1, 11, 21.23)
    y1 = pymeeus.base.Angle(0, -28, 12.31)
    y2 = pymeeus.base.Angle(0, 16, 7.02)
    y3 = pymeeus.base.Angle(1, 1, 0.13)
    y4 = pymeeus.base.Angle(1, 45, 46.33)

    i_angles2.set([25.0, 26.0, 27.0, 28.0, 29.0], [y0, y1, y2, y3, y4])

    # Set up an interpolation object with 6 interpolation table entries, based
    # on sine function
    i_sine.set([29.43, 30.97, 27.69, 28.11, 31.58, 33.05],
               [0.4913598528, 0.5145891926, 0.4646875083, 0.4711658342,
                0.5236885653, 0.5453707057])

    # Set up a few CurveFitting objects
    cf1.set([73.0, 38.0, 35.0, 42.0, 78.0, 68.0, 74.0, 42.0, 52.0, 54.0, 39.0,
             61.0, 42.0, 49.0, 50.0, 62.0, 44.0, 39.0, 43.0, 54.0, 44.0, 37.0],
            [90.4, 125.3, 161.8, 143.4, 52.5, 50.8, 71.5, 152.8, 131.3, 98.5,
             144.8, 78.1, 89.5, 63.9, 112.1, 82.0, 119.8, 161.2, 208.4, 111.6,
             167.1, 162.1])

    cf2.set([0.2982, 0.2969, 0.2918, 0.2905, 0.2707, 0.2574, 0.2485, 0.2287,
             0.2238, 0.2156, 0.1992, 0.1948, 0.1931, 0.1889, 0.1781, 0.1772,
             0.1770, 0.1755, 0.1746],
            [10.92, 11.01, 10.99, 10.78, 10.87, 10.80, 10.75, 10.14, 10.21,
             9.97, 9.69, 9.57, 9.66, 9.63, 9.65, 9.44, 9.44, 9.32, 9.20])

    cf3.set([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            [-9.372, -3.821, 0.291, 3.730, 5.822, 8.324, 9.083, 6.957, 7.006,
             0.365, -1.722])

    cf4.set([3, 20, 34, 50, 75, 88, 111, 129, 143, 160, 183, 200, 218, 230,
             248, 269, 290, 303, 320, 344],
            [0.0433, 0.2532, 0.3386, 0.3560, 0.4983, 0.7577, 1.4585, 1.8628,
             1.8264, 1.2431, -0.2043, -1.2431, -1.8422, -1.8726, -1.4889,
             -0.8372, -0.4377, -0.3640, -0.3508, -0.2126])


def teardown():
    pass


# Angle class

def test_angle_constructor():
    """Tests the constructor of Angle class"""

    a = pymeeus.base.Angle(23.44694444)
    assert abs(a._deg - 23.44694444) < TOLERANCE, \
        "ERROR: 1st constructor test, degrees value doesn't match"

    b = pymeeus.base.Angle(-23.44694444)
    assert abs(b._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 2nd constructor test, degrees value doesn't match"

    c = pymeeus.base.Angle(383.44694444)
    assert abs(c._deg - 23.44694444) < TOLERANCE, \
        "ERROR: 3rd constructor test, degrees value doesn't match"

    d = pymeeus.base.Angle(-383.44694444)
    assert abs(d._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 4th constructor test, degrees value doesn't match"

    e = pymeeus.base.Angle(-23.0, 26.0, 48.999983999)
    assert abs(e._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 5th constructor test, degrees value doesn't match"

    f = pymeeus.base.Angle(23.0, -30.0)
    assert abs(f._deg - (-23.5)) < TOLERANCE, \
        "ERROR: 6th constructor test, degrees value doesn't match"

    g = pymeeus.base.Angle((-23.0, 26.0, 48.999983999))
    assert abs(g._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 7th constructor test, degrees value doesn't match"

    h = pymeeus.base.Angle([-23.0, 26.0, 48.999983999])
    assert abs(h._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 8th constructor test, degrees value doesn't match"

    i = pymeeus.base.Angle(1.0, radians=True)
    assert abs(i._deg - 57.29577951308232) < TOLERANCE, \
        "ERROR: 9th constructor test, degrees value doesn't match"

    j = pymeeus.base.Angle((23.0, 26.0, 48.999983999, -1.0))
    assert abs(j._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 10th constructor test, degrees value doesn't match"

    k = pymeeus.base.Angle(23.0, 26.0, 48.999983999, -7.4)
    assert abs(k._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 11th constructor test, degrees value doesn't match"

    m = pymeeus.base.Angle([23.0, -26.0, 48.999983999, -4.5])
    assert abs(m._deg - (-23.44694444)) < TOLERANCE, \
        "ERROR: 12th constructor test, degrees value doesn't match"


def test_angle_set_radians():
    """Tests the set_radians() method of Angle class"""

    a = pymeeus.base.Angle()

    a.set_radians(pi)                               # Input is in radians
    assert abs(a() - 180.0) < TOLERANCE, \
        "ERROR: 1st set_radians() test, degrees value doesn't match"


def test_angle_set_ra():
    """Tests the set_ra() method of Angle class"""

    a = pymeeus.base.Angle()

    a.set_ra(9.248833333333)                        # Input is in RA
    assert abs(a._deg - 138.7325) < TOLERANCE, \
        "ERROR: 1st set_ra() test, degrees value doesn't match"

    a.set_ra(-9.248833333333)
    assert abs(a._deg - (-138.7325)) < TOLERANCE, \
        "ERROR: 2nd set_ra() test, degrees value doesn't match"

    a.set_ra(9, 14, 55.8)
    assert abs(a._deg - 138.7325) < TOLERANCE, \
        "ERROR: 3rd set_ra() test, degrees value doesn't match"

    a.set_ra((9, 14, 55.8, -1.0))
    assert abs(a._deg - (-138.7325)) < TOLERANCE, \
        "ERROR: 4th set_ra() test, degrees value doesn't match"


def test_angle_deg2dms():
    """Tests deg2dms() static method of Angle class"""

    (d, m, s, sign) = pymeeus.base.Angle.deg2dms(23.44694444)
    assert abs(d - 23.0) < TOLERANCE, \
        "ERROR: In 1st deg2dms() test, degrees value doesn't match"
    assert abs(m - 26.0) < TOLERANCE, \
        "ERROR: In 1st deg2dms() test, minutes value doesn't match"
    assert abs(s - 48.9999839999999) < TOLERANCE, \
        "ERROR: In 1st deg2dms() test, seconds value doesn't match"
    assert abs(sign - 1.0) < TOLERANCE, \
        "ERROR: In 1st deg2dms() test, sign value doesn't match"

    (d, m, s, sign) = pymeeus.base.Angle.deg2dms(-23.44694444)
    assert abs(d - 23.0) < TOLERANCE, \
        "ERROR: In 2nd deg2dms() test, degrees value doesn't match"
    assert abs(m - 26.0) < TOLERANCE, \
        "ERROR: In 2nd deg2dms() test, minutes value doesn't match"
    assert abs(s - 48.9999839999999) < TOLERANCE, \
        "ERROR: In 2nd deg2dms() test, seconds value doesn't match"
    assert abs(sign - (-1.0)) < TOLERANCE, \
        "ERROR: In 2nd deg2dms() test, sign value doesn't match"


def test_angle_dms2deg():
    """Tests dms2deg() static method of Angle class"""

    d = pymeeus.base.Angle.dms2deg(23.0, 26.0, 48.999984)
    assert abs(d - 23.44694444) < TOLERANCE, \
        "ERROR: In 1st dms2deg() test, degrees value doesn't match"

    d = pymeeus.base.Angle.dms2deg(-23.0, 26.0, 48.999984)
    assert abs(d - (-23.44694444)) < TOLERANCE, \
        "ERROR: In 2nd dms2deg() test, degrees value doesn't match"

    d = pymeeus.base.Angle.dms2deg(0.0, -26.0, 48.999984)
    assert abs(d - (-0.44694444)) < TOLERANCE, \
        "ERROR: In 3rd dms2deg() test, degrees value doesn't match"


def test_angle_reduce_deg():
    """Tests reduce_deg() static method of Angle class"""

    d = pymeeus.base.Angle.reduce_deg(745.67)
    assert abs(d - 25.67) < TOLERANCE, \
        "ERROR: In 1st reduce_deg() test, degrees value doesn't match"

    d = pymeeus.base.Angle.reduce_deg(-360.86)
    assert abs(d - (-0.86)) < TOLERANCE, \
        "ERROR: In 2nd reduce_deg() test, degrees value doesn't match"


def test_angle_reduce_dms():
    """Tests reduce_dms() static method of Angle class"""

    (d, m, s, sign) = pymeeus.base.Angle.reduce_dms(383.0, 26.0, 48.999984)
    assert abs(d - 23.0) < TOLERANCE, \
        "ERROR: In 1st reduce_dms() test, degrees value doesn't match"

    assert abs(m - 26.0) < TOLERANCE, \
        "ERROR: In 1st reduce_dms() test, minutes value doesn't match"

    assert abs(s - 48.999984) < TOLERANCE, \
        "ERROR: In 1st reduce_dms() test, seconds value doesn't match"

    assert abs(sign - 1.0) < TOLERANCE, \
        "ERROR: In 1st reduce_dms() test, sign value doesn't match"

    (d, m, s, sign) = pymeeus.base.Angle.reduce_dms(-1103.6, 86.5, 168.84)
    assert abs(d - 25.0) < TOLERANCE, \
        "ERROR: In 2nd reduce_dms() test, degrees value doesn't match"

    assert abs(m - 5.0) < TOLERANCE, \
        "ERROR: In 2nd reduce_dms() test, minutes value doesn't match"

    assert abs(s - 18.8399999997) < TOLERANCE, \
        "ERROR: In 2nd reduce_dms() test, seconds value doesn't match"

    assert abs(sign - (-1.0)) < TOLERANCE, \
        "ERROR: In 2nd reduce_dms() test, sign value doesn't match"

    (d, m, s, sign) = pymeeus.base.Angle.reduce_dms(0.0, -206.71)
    assert abs(d - 3.0) < TOLERANCE, \
        "ERROR: In 3rd reduce_dms() test, degrees value doesn't match"

    assert abs(m - 26.0) < TOLERANCE, \
        "ERROR: In 3rd reduce_dms() test, minutes value doesn't match"

    assert abs(s - 42.6) < TOLERANCE, \
        "ERROR: In 3rd reduce_dms() test, seconds value doesn't match"

    assert abs(sign - (-1.0)) < TOLERANCE, \
        "ERROR: In 3rd reduce_dms() test, sign value doesn't match"


def test_angle_dms_str():
    """Tests dms_str() method of Angle class"""
    a = pymeeus.base.Angle(0, -46.25, 0.0)

    result = a.dms_str()
    assert result == "-46' 15.0''", \
        "ERROR: In 1st dms_str() test, the output value doesn't match"

    result = a.dms_str(False)
    assert result == "0:-46:15.0", \
        "ERROR: In 2nd dms_str() test, the output value doesn't match"


def test_angle_ra_str():
    """Tests ra_str() method of Angle class"""
    a = pymeeus.base.Angle(138.75)

    result = a.ra_str()
    assert result == "9h 15' 0.0''", \
        "ERROR: In 1st ra_str() test, the output value doesn't match"

    result = a.ra_str(False)
    assert result == "9:15:0.0", \
        "ERROR: In 2nd ra_str() test, the output value doesn't match"


def test_angle_get_ra():
    """Tests get_ra() method of Angle class"""
    a = pymeeus.base.Angle(138.75)

    assert abs(a.get_ra() - 9.25) < TOLERANCE, \
        "ERROR: In 1st get_ra() test, the output value doesn't match"


def test_angle_call():
    """Tests the __call__() method of Angle class"""

    type_ok = False
    a = pymeeus.base.Angle(40, -46.25, 0.0)
    if isinstance(a(), (int, float)):       # Test the returned type
        type_ok = True

    assert type_ok, "ERROR: In 1st __call__() test, type doesn't match"

    assert abs(a() - (-40.770833333333)) < TOLERANCE, \
        "ERROR: In 2nd __call__() test, degrees value doesn't match"


def test_angle_str():
    """Tests the __str__() method of Angle class"""

    type_ok = False
    a = pymeeus.base.Angle(40, -46.25, 0.0)
    if isinstance(a.__str__(), str):                  # Test the returned type
        type_ok = True

    assert type_ok, "ERROR: In 1st __str__() test, type doesn't match"

    assert a.__str__() == "-40.7708333333", \
        "ERROR: In 2nd __str__() test, degrees value doesn't match"


def test_angle_rad():
    """Tests the rad() method of Angle class"""
    a = pymeeus.base.Angle(180.0)

    assert abs(a.rad() - pi) < TOLERANCE, \
        "ERROR: In 1st rad() test, radians value doesn't match"


def test_angle_to_positive():
    """Tests the to_positive() method"""
    a = pymeeus.base.Angle(-87.32)
    b = pymeeus.base.Angle(87.32)

    assert abs(a.to_positive()() - 272.68) < TOLERANCE, \
        "ERROR: In 1st to_positive() test, value doesn't match"

    assert abs(b.to_positive()() - 87.32) < TOLERANCE, \
        "ERROR: In 2nd to_positive() test, value doesn't match"


def test_angle_ne():
    """Tests the 'is not equal' operator of Angles"""
    # NOTE: Test 'is not equal' also tests 'is equal' operator
    # Default tolerance for Angles is 1E-10
    a = pymeeus.base.Angle(152.7)
    b = pymeeus.base.Angle(152.7000000001)

    assert (a != b), \
        "ERROR: In 1st __ne__() test, Angles are different but taken as equal"

    a = pymeeus.base.Angle(-13, 30)
    b = pymeeus.base.Angle(-13.50000000001)

    assert not (a != b), \
        "ERROR: In 2nd __ne__() test, Angles are equal but taken as different"


def test_angle_ge():
    """Tests the 'is greater or equal' operator of Angles"""
    # NOTE: Test of 'is greater or equal' also test 'is less than' operator
    a = pymeeus.base.Angle(152.7)
    b = pymeeus.base.Angle(152.70000001)

    assert not (a >= b), \
        "ERROR: In 1st __ge__() test, Angles values don't match operator"

    a = pymeeus.base.Angle(-13, 30)
    b = pymeeus.base.Angle(-13.5)

    assert (a >= b), \
        "ERROR: In 2nd __ne__() test, Angles values don't match operator"


def test_angle_le():
    """Tests the 'is less or equal' operator of Angles"""
    # NOTE: Test of 'is less or equal' also test 'is greater than' operator
    a = pymeeus.base.Angle(152.7)
    b = pymeeus.base.Angle(152.70000001)

    assert (a <= b), \
        "ERROR: In 1st __le__() test, Angles values don't match operator"

    a = pymeeus.base.Angle(-13, 30)
    b = pymeeus.base.Angle(-13.5)

    assert (a <= b), \
        "ERROR: In 2nd __le__() test, Angles values don't match operator"


def test_angle_neg():
    """Tests the negation of Angles"""
    a = pymeeus.base.Angle(152.7)

    b = -a

    assert abs(b() - (-152.7)) < TOLERANCE, \
        "ERROR: In 1st __neg__() test, degrees value doesn't match"

    a = pymeeus.base.Angle(-13, 30)
    b = -a

    assert abs(b() - 13.5) < TOLERANCE, \
        "ERROR: In 2nd __neg__() test, degrees value doesn't match"


def test_angle_abs():
    """Tests the absolute value of Angles"""
    a = pymeeus.base.Angle(152.7)

    b = abs(a)

    assert abs(b() - 152.7) < TOLERANCE, \
        "ERROR: In 1st __abs__() test, degrees value doesn't match"

    a = pymeeus.base.Angle(-13, 30)
    b = abs(a)

    assert abs(b() - 13.5) < TOLERANCE, \
        "ERROR: In 2nd __abs__() test, degrees value doesn't match"


def test_angle_mod():
    """Tests the module of Angles"""
    a = pymeeus.base.Angle(152.7)

    b = a % 50

    assert abs(b() - 2.7) < TOLERANCE, \
        "ERROR: In 1st __mod__() test, degrees value doesn't match"

    a = pymeeus.base.Angle(-13, 30)
    b = a % 7

    assert abs(b() - (-6.5)) < TOLERANCE, \
        "ERROR: In 2nd __mod__() test, degrees value doesn't match"


def test_angle_add():
    """Tests the addition between Angles"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(13, 30)
    c = a + b

    assert abs(c() - 193.5) < TOLERANCE, \
        "ERROR: In 1st __add__() test, degrees value doesn't match"

    b.set(-13, 30)
    c = a + b

    assert abs(c() - 166.5) < TOLERANCE, \
        "ERROR: In 2nd __add__() test, degrees value doesn't match"

    c = a + 11.5

    assert abs(c() - 191.5) < TOLERANCE, \
        "ERROR: In 3rd __add__() test, degrees value doesn't match"


def test_angle_sub():
    """Tests the subtraction between Angles"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(13, 30)
    c = a - b

    assert abs(c() - 166.5) < TOLERANCE, \
        "ERROR: In 1st __sub__() test, degrees value doesn't match"

    b.set(-13, 30)
    c = a - b

    assert abs(c() - 193.5) < TOLERANCE, \
        "ERROR: In 2nd __sub__() test, degrees value doesn't match"

    c = a - 11.5

    assert abs(c() - 168.5) < TOLERANCE, \
        "ERROR: In 3rd __sub__() test, degrees value doesn't match"


def test_angle_mul():
    """Tests the multiplication between Angles"""
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(5.0)
    c = a * b

    assert abs(c() - 30.0) < TOLERANCE, \
        "ERROR: In 1st __mul__() test, degrees value doesn't match"

    b.set(-5.0)
    c = a * b

    assert abs(c() - (-30.0)) < TOLERANCE, \
        "ERROR: In 2nd __mul__() test, degrees value doesn't match"

    c = a * 2.5

    assert abs(c() - 15.0) < TOLERANCE, \
        "ERROR: In 3rd __mul__() test, degrees value doesn't match"


def test_angle_div():
    """Tests the division between Angles"""
    # NOTE: This also tests method self.__truediv__()
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(6.0)
    c = a / b

    assert abs(c() - 25.0) < TOLERANCE, \
        "ERROR: In 1st __div__() test, degrees value doesn't match"

    b.set(-6.0)
    c = a / b

    assert abs(c() - (-25.0)) < TOLERANCE, \
        "ERROR: In 2nd __div__() test, degrees value doesn't match"

    c = a / 1.5

    assert abs(c() - 100.0) < TOLERANCE, \
        "ERROR: In 3rd __div__() test, degrees value doesn't match"


def test_angle_pow():
    """Tests the power operation in Angles"""
    a = pymeeus.base.Angle(13, 30)
    b = a ** 3

    assert abs(b() - 300.375) < TOLERANCE, \
        "ERROR: In 1st __pow__() test, degrees value doesn't match"

    b.set(-2.0)
    c = a ** b

    assert abs(c() - 0.005486968449931413) < TOLERANCE, \
        "ERROR: In 2nd __pow__() test, degrees value doesn't match"


def test_angle_imod():
    """Tests the accumulative module between Angles"""
    a = pymeeus.base.Angle(152.7)

    a %= 50

    assert abs(a() - 2.7) < TOLERANCE, \
        "ERROR: In 1st __imod__() test, degrees value doesn't match"

    a = pymeeus.base.Angle(-13, 30)
    a %= 7

    assert abs(a() - (-6.5)) < TOLERANCE, \
        "ERROR: In 2nd __imod__() test, degrees value doesn't match"


def test_angle_iadd():
    """Tests the accumulative addition between Angles"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(13, 30)
    a += b

    assert abs(a() - 193.5) < TOLERANCE, \
        "ERROR: In 1st __iadd__() test, degrees value doesn't match"

    b.set(-10, 30)
    a += b

    assert abs(a() - 183.0) < TOLERANCE, \
        "ERROR: In 2nd __iadd__() test, degrees value doesn't match"

    a += 37.5

    assert abs(a() - 220.5) < TOLERANCE, \
        "ERROR: In 3rd __iadd__() test, degrees value doesn't match"


def test_angle_isub():
    """Tests the accumulative subtraction between Angles"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(13, 30)
    a -= b

    assert abs(a() - 166.5) < TOLERANCE, \
        "ERROR: In 1st __isub__() test, degrees value doesn't match"

    b.set(-10, 30)
    a -= b

    assert abs(a() - 177.0) < TOLERANCE, \
        "ERROR: In 2nd __isub__() test, degrees value doesn't match"

    a -= 37.5

    assert abs(a() - 139.5) < TOLERANCE, \
        "ERROR: In 3rd __isub__() test, degrees value doesn't match"


def test_angle_imul():
    """Tests the accumulative multiplication between Angles"""
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(5.0)
    a *= b

    assert abs(a() - 30.0) < TOLERANCE, \
        "ERROR: In 1st __imul__() test, degrees value doesn't match"

    b.set(-5.0)
    a *= b

    assert abs(a() - (-150.0)) < TOLERANCE, \
        "ERROR: In 2nd __imul__() test, degrees value doesn't match"

    a *= 2.5

    assert abs(a() - (-15.0)) < TOLERANCE, \
        "ERROR: In 3rd __imul__() test, degrees value doesn't match"


def test_angle_idiv():
    """Tests the accumulative division between Angles"""
    # NOTE: This also tests method self.__itruediv__()
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(6.0)
    a /= b

    assert abs(a() - 25.0) < TOLERANCE, \
        "ERROR: In 1st __idiv__() test, degrees value doesn't match"

    b.set(-20.0)
    a /= b

    assert abs(a() - (-1.25)) < TOLERANCE, \
        "ERROR: In 2nd __idiv__() test, degrees value doesn't match"

    a /= 1.5

    assert abs(a() - (-0.833333333333333)) < TOLERANCE, \
        "ERROR: In 3rd __idiv__() test, degrees value doesn't match"


def test_angle_ipow():
    """Tests the accumulative power operation in Angles"""
    a = pymeeus.base.Angle(13, 30)
    a **= 3

    assert abs(a() - 300.375) < TOLERANCE, \
        "ERROR: In 1st __ipow__() test, degrees value doesn't match"

    b = pymeeus.base.Angle(-2.0)
    a **= b

    assert abs(a() - 1.108338532999e-05) < TOLERANCE, \
        "ERROR: In 2nd __ipow__() test, degrees value doesn't match"


def test_angle_rmod():
    """Tests the module operation between Angles by the right"""
    a = pymeeus.base.Angle(25.0)
    b = pymeeus.base.Angle(163.0)
    c = a.__rmod__(b)

    assert abs(c() - 13.0) < TOLERANCE, \
        "ERROR: In 1st __rmod__() test, degrees value doesn't match"

    b.set(-78.0)
    c = a.__rmod__(b)

    assert abs(c() - (-3.0)) < TOLERANCE, \
        "ERROR: In 2nd __rmod__() test, degrees value doesn't match"

    c = 31.5 % a

    assert abs(c() - 6.5) < TOLERANCE, \
        "ERROR: In 3rd __rmod__() test, degrees value doesn't match"


def test_angle_radd():
    """Tests the addition between Angles by the right"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(13, 30)
    c = a.__radd__(b)

    assert abs(c() - 193.5) < TOLERANCE, \
        "ERROR: In 1st __radd__() test, degrees value doesn't match"

    b.set(-13, 30)
    c = a.__radd__(b)

    assert abs(c() - 166.5) < TOLERANCE, \
        "ERROR: In 2nd __radd__() test, degrees value doesn't match"

    c = 11.5 + a

    assert abs(c() - 191.5) < TOLERANCE, \
        "ERROR: In 3rd __radd__() test, degrees value doesn't match"


def test_angle_rsub():
    """Tests the subtraction between Angles by the right"""
    a = pymeeus.base.Angle(180.0)
    b = pymeeus.base.Angle(-13, 30)
    c = a.__rsub__(b)

    assert abs(c() - (-193.5)) < TOLERANCE, \
        "ERROR: In 1st __rsub__() test, degrees value doesn't match"

    b.set(13, 30)
    c = a.__rsub__(b)

    assert abs(c() - (-166.5)) < TOLERANCE, \
        "ERROR: In 2nd __rsub__() test, degrees value doesn't match"

    c = 11.5 - a

    assert abs(c() - (-168.5)) < TOLERANCE, \
        "ERROR: In 3rd __rsub__() test, degrees value doesn't match"


def test_angle_rmul():
    """Tests the multiplication between Angles by the right"""
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(5.0)
    c = a.__rmul__(b)

    assert abs(c() - 30.0) < TOLERANCE, \
        "ERROR: In 1st __rmul__() test, degrees value doesn't match"

    b.set(-5.0)
    c = a.__rmul__(b)

    assert abs(c() - (-30.0)) < TOLERANCE, \
        "ERROR: In 2nd __rmul__() test, degrees value doesn't match"

    c = 2.5 * a

    assert abs(c() - 15.0) < TOLERANCE, \
        "ERROR: In 3rd __rmul__() test, degrees value doesn't match"


def test_angle_rdiv():
    """Tests the division between Angles by the right"""
    a = pymeeus.base.Angle(150.0)
    b = pymeeus.base.Angle(5.0)
    c = b.__rdiv__(a)

    assert abs(c() - 30.0) < TOLERANCE, \
        "ERROR: In 1st __rdiv__() test, degrees value doesn't match"

    b.set(-5.0)
    a = b.__rdiv__(c)

    assert abs(a() - (-6.0)) < TOLERANCE, \
        "ERROR: In 2nd __rdiv__() test, degrees value doesn't match"

    c = -24.0 / a

    assert abs(c() - 4.0) < TOLERANCE, \
        "ERROR: In 3rd __rdiv__() test, degrees value doesn't match"


def test_angle_rpow():
    """Tests the power operation between Angles by the right"""
    a = pymeeus.base.Angle(15.0)
    b = pymeeus.base.Angle(3.0)
    c = b.__rpow__(a)

    assert abs(c() - 135.0) < TOLERANCE, \
        "ERROR: In 1st __rpow__() test, degrees value doesn't match"

    b.set(-2.0)
    a = b.__rpow__(c)

    assert abs(a() - 5.48697e-05) < TOLERANCE, \
        "ERROR: In 2nd __rpow__() test, degrees value doesn't match"

    c = -10.0 ** b

    assert abs(c() - (-0.01)) < TOLERANCE, \
        "ERROR: In 3rd __rpow__() test, degrees value doesn't match"


def test_angle_float():
    """Tests the 'float()' operation on Angles"""
    a = pymeeus.base.Angle(15, 30)

    assert abs(float(a) - 15.5) < TOLERANCE, \
        "ERROR: In 1st __float__() test, degrees value doesn't match"


def test_angle_int():
    """Tests the 'int()' operation on Angles"""
    a = pymeeus.base.Angle(15, 30)

    assert abs(int(a) - 15) < TOLERANCE, \
        "ERROR: In 1st __int__() test, degrees value doesn't match"


def test_angle_round():
    """Tests the 'round()' operation on Angles"""
    a = pymeeus.base.Angle(1.0, radians=True)

    # NOTE: The 'float(round(x))' hack makes this test work in Python 2 and 3
    assert abs(float(round(a)) - 57.0) < TOLERANCE, \
        "ERROR: In 1st __round__() test, degrees value doesn't match"

    assert abs(float(round(a, 3)) - 57.296) < TOLERANCE, \
        "ERROR: In 2nd __round__() test, degrees value doesn't match"

    assert abs(float(round(a, 7)) - 57.2957795) < TOLERANCE, \
        "ERROR: In 3rd __round__() test, degrees value doesn't match"


# Interpolation class

def test_interpolation_constructor():
    """Tests the constructor of Interpolation class"""

    i = pymeeus.base.Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
    assert i._x == [1, 2, 3, 4, 5, 6], \
        "ERROR: 1st constructor test, 'x' values don't match"

    assert i._y == [2, 4, 6, 8, 10, 12], \
        "ERROR: 2nd constructor test, 'y' values don't match"

    j = pymeeus.base.Interpolation([3, -8, 1, 12, 2, 5, 8])
    assert j._x == [0, 1, 2, 3, 4, 5, 6], \
        "ERROR: 3rd constructor test, 'x' values don't match"

    assert j._y == [3, -8, 1, 12, 2, 5, 8], \
        "ERROR: 4th constructor test, 'y' values don't match"

    k = pymeeus.base.Interpolation(3, -8, 1, 12, 2, 5, 8)
    assert k._x == [1, 2, 3], \
        "ERROR: 5th constructor test, 'x' values don't match"

    assert k._y == [12, 5, -8], \
        "ERROR: 6th constructor test, 'y' values don't match"


def test_interpolation_call():
    """Tests the __call__() method of Interpolation class"""

    m = pymeeus.base.Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])

    assert abs(m(-0.8) - (-0.52)) < TOLERANCE, \
        "ERROR: In 1st __call__() test, output value doesn't match"

    assert abs(m(0.7) - 2.93) < TOLERANCE, \
        "ERROR: In 2nd __call__() test, output value doesn't match"

    assert abs(m(-1.0) - (-2.0)) < TOLERANCE, \
        "ERROR: In 3rd __call__() test, output value doesn't match"

    m = pymeeus.base.Interpolation([-3.0, 0.0, 2.5], [12.0, -3.0, -1.75])

    assert abs(m(-2.0) - 5.0) < TOLERANCE, \
        "ERROR: In 4th __call__() test, output value doesn't match"

    assert abs(m(2.5) - (-1.75)) < TOLERANCE, \
        "ERROR: In 5th __call__() test, output value doesn't match"

    # This interpolation test uses Right Ascension
    a = pymeeus.base.Angle(i_ra(11.0))
    h, m, s, sign = a.ra_tuple()
    assert abs(h - 10.0) < TOLERANCE and \
        abs(m - 25.0) < TOLERANCE and \
        abs(s - 40.0014375) < TOLERANCE and \
        abs(sign == 1.0), \
        "ERROR: In 6th __call__() test, output value doesn't match"

    # Test with 6 interpolation table entries, based on sine function
    assert abs(i_sine(30.0) - 0.5) < TOLERANCE, \
        "ERROR: In 7th __call__() test, output value doesn't match"


def test_interpolation_derivative():
    """Tests the derivative() method of Interpolation class"""

    m = pymeeus.base.Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])

    assert abs(m.derivative(-1.0) - 8.0) < TOLERANCE, \
        "ERROR: In 1st derivative() test, output value doesn't match"

    assert abs(m.derivative(0.0) - 2.0) < TOLERANCE, \
        "ERROR: In 2nd derivative() test, output value doesn't match"

    assert abs(m.derivative(0.5) - (-1.0)) < TOLERANCE, \
        "ERROR: In 3rd derivative() test, output value doesn't match"

    m = pymeeus.base.Interpolation([-3.0, 0.0, 2.5], [12.0, -3.0, -1.75])

    assert abs(m.derivative(-3.0) - (-8.0)) < TOLERANCE, \
        "ERROR: In 4th derivative() test, output value doesn't match"

    assert abs(m.derivative(0.0) - (-2.0)) < TOLERANCE, \
        "ERROR: In 5th derivative() test, output value doesn't match"

    assert abs(m.derivative(2.5) - 3.0) < TOLERANCE, \
        "ERROR: In 6th derivative() test, output value doesn't match"

    # Do test with an interpolation object with 6 table entries, based on sine
    # We need to adjust the result because degrees were used instead of radians
    res = degrees(i_sine.derivative(30.0))
    assert abs(res - sqrt(3.0)/2.0) < TOLERANCE, \
        "ERROR: In 7th derivative() test, output value doesn't match"


def test_interpolation_root():
    """Tests the root() method of Interpolation class"""

    m = pymeeus.base.Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])

    assert abs(m.root() - (-0.7207592200561265)) < TOLERANCE, \
        "ERROR: In 1st root() test, output value doesn't match"

    m = pymeeus.base.Interpolation([-3.0, 0.0, 2.5], [12.0, -3.0, -1.75])

    assert abs(m.root(-2.0, 0.0) - (-1.0)) < TOLERANCE, \
        "ERROR: In 2nd root() test, output value doesn't match"

    assert abs(m.root() - (-1.0)) < TOLERANCE, \
        "ERROR: In 3rd root() test, output value doesn't match"

    m = pymeeus.base.Interpolation([-3.0, 0.0, 2.5, 3.5],
                                   [12.0, -3.0, -1.75, 2.25])

    assert abs(m.root(0.0, 3.15) - 3.0) < TOLERANCE, \
        "ERROR: In 4th root() test, output value doesn't match"

    # Let's do some tests with Angles
    assert abs(i_angles1.root() - 26.798732705) < TOLERANCE, \
        "ERROR: In 5th root() test, output value doesn't match"

    assert abs(i_angles2.root() - 26.6385869469) < TOLERANCE, \
        "ERROR: In 6th root() test, output value doesn't match"


def test_interpolation_minmax():
    """Tests the minmax() method of Interpolation class"""

    m = pymeeus.base.Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])

    assert abs(m.minmax() - 0.3333333333) < TOLERANCE, \
        "ERROR: In 1st minmax() test, output value doesn't match"

    m = pymeeus.base.Interpolation([-3.0, 0.0, 2.5], [12.0, -3.0, -1.75])

    assert abs(m.minmax() - 1.0) < TOLERANCE, \
        "ERROR: In 2nd minmax() test, output value doesn't match"

    m = pymeeus.base.Interpolation([12.0, 16.0, 20.0],
                                   [1.3814294, 1.3812213, 1.3812453])

    assert abs(m.minmax() - 17.5863851788) < TOLERANCE, \
        "ERROR: In 3rd minmax() test, output value doesn't match"

    assert abs(m(m.minmax()) - 1.38120304666) < TOLERANCE, \
        "ERROR: In 4th minmax() test, output value doesn't match"


# CurveFitting class

def test_curvefitting_constructor():
    """Tests the constructor of CurveFitting class"""

    i = pymeeus.base.CurveFitting([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
    assert i._x == [5, 3, 6, 1, 2, 4], \
        "ERROR: 1st constructor test, 'x' values don't match"

    assert i._y == [10, 6, 12, 2, 4, 8], \
        "ERROR: 2nd constructor test, 'y' values don't match"

    j = pymeeus.base.CurveFitting([3, -8, 1, 12, 2, 5, 8])
    assert j._x == [0, 1, 2, 3, 4, 5, 6], \
        "ERROR: 3rd constructor test, 'x' values don't match"

    assert j._y == [3, -8, 1, 12, 2, 5, 8], \
        "ERROR: 4th constructor test, 'y' values don't match"

    k = pymeeus.base.CurveFitting(3, -8, 1, 12, 2, 5, 8)
    assert k._x == [3, 1, 2], \
        "ERROR: 5th constructor test, 'x' values don't match"

    assert k._y == [-8, 12, 5], \
        "ERROR: 6th constructor test, 'y' values don't match"


def test_curvefitting_correlation_coeff():
    """Tests the correlation_coeff() method of CurveFitting class"""

    r = cf1.correlation_coeff()
    assert abs(round(r, 3) - (-0.767)) < TOLERANCE, \
        "ERROR: 1st correlation_coeff() test, 'r' value doesn't match"


def test_curvefitting_linear_fitting():
    """Tests the linear_fitting() method of CurveFitting class"""

    a, b = cf1.linear_fitting()
    assert abs(round(a, 2) - (-2.49)) < TOLERANCE, \
        "ERROR: In 1st linear_fitting() test, 'a' value doesn't match"

    assert abs(round(b, 2) - 244.18) < TOLERANCE, \
        "ERROR: In 2nd linear_fitting() test, 'b' value doesn't match"

    a, b = cf2.linear_fitting()
    assert abs(round(a, 2) - 13.67) < TOLERANCE, \
        "ERROR: In 3rd linear_fitting() test, 'a' value doesn't match"

    assert abs(round(b, 2) - 7.03) < TOLERANCE, \
        "ERROR: In 4th linear_fitting() test, 'b' value doesn't match"


def test_curvefitting_quadratic_fitting():
    """Tests the quadratic_fitting() method of CurveFitting class"""

    a, b, c = cf3.quadratic_fitting()
    assert abs(round(a, 2) - (-2.22)) < TOLERANCE, \
        "ERROR: In 1st quadratic_fitting() test, 'a' value doesn't match"

    assert abs(round(b, 2) - 3.76) < TOLERANCE, \
        "ERROR: In 2nd quadratic_fitting() test, 'b' value doesn't match"

    assert abs(round(c, 2) - 6.64) < TOLERANCE, \
        "ERROR: In 3rd quadratic_fitting() test, 'c' value doesn't match"


def test_curvefitting_general_fitting():
    """Tests the general_fitting() method of CurveFitting class"""

    # Let's define the three functions to be used for fitting
    def sin1(x): return sin(radians(x))

    def sin2(x): return sin(radians(2.0*x))

    def sin3(x): return sin(radians(3.0*x))

    a, b, c = cf4.general_fitting(sin1, sin2, sin3)
    assert abs(round(a, 2) - 1.2) < TOLERANCE, \
        "ERROR: In 1st general_fitting() test, 'a' value doesn't match"

    assert abs(round(b, 2) - (-0.77)) < TOLERANCE, \
        "ERROR: In 2nd general_fitting() test, 'b' value doesn't match"

    assert abs(round(c, 2) - 0.39) < TOLERANCE, \
        "ERROR: In 3rd general_fitting() test, 'c' value doesn't match"

    cf5 = pymeeus.base.CurveFitting([0, 1.2, 1.4, 1.7, 2.1, 2.2])

    a, b, c = cf5.general_fitting(sqrt)
    assert abs(round(a, 3) - 1.016) < TOLERANCE, \
        "ERROR: In 4th general_fitting() test, 'a' value doesn't match"

    assert abs(round(b, 3) - 0.0) < TOLERANCE, \
        "ERROR: In 5th general_fitting() test, 'b' value doesn't match"

    assert abs(round(c, 3) - 0.0) < TOLERANCE, \
        "ERROR: In 6th general_fitting() test, 'c' value doesn't match"
