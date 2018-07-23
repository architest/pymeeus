Angle examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Create an Angle object, providing degrees, minutes and seconds::

    a = Angle(-23.0, 26.0, 48.999983999)

First we print using the ``__call__()`` method (note the extra parentheses)::

    print_me("The angle 'a()' is", a())

    # The angle 'a()' is: -23.44694444

Second we print using the ``__str__()`` method (no extra parentheses needed)::

    print_me("The angle 'a' is", a)

    # The angle 'a' is: -23.44694444

Use the copy constructor::

    b = Angle(a)

    print_me("Angle 'b', which is a copy of 'a', is", b)

    # Angle 'b', which is a copy of 'a', is: -23.44694444

Use the static ``deg2dms()`` method to carry out conversions::

    d, m, s, sign = Angle.deg2dms(23.44694444)

    val = "{}d {}' {}''".format(int(sign*d), m, s)

    print_me("{Deg}d {Min}' {Sec}''", val)

    # {Deg}d {Min}' {Sec}'': 23d 26' 48.999984''

We can print Angle ``a`` directly in sexagesimal format.

- In *fancy* format::

    print_me("{Deg}d {Min}' {Sec}''", a.dms_str())

    # {Deg}d {Min}' {Sec}'': -23d 26' 48.999983999''

- In plain format::

    print_me("{Deg}:{Min}:{Sec}", a.dms_str(False))

    # {Deg}:{Min}:{Sec}: -23:26:48.999983999

Print directly as a tuple::

    a = Angle(23.44694444)

    print_me("a.dms_tuple()", a.dms_tuple())

    # a.dms_tuple(): (23, 26, 48.999983999997596, 1.0)

    print_me("a.ra_tuple()", a.ra_tuple())

    # a.ra_tuple(): (1, 33, 47.26666559999941, 1.0)

Redefine Angle ``a`` several times::

    a.set(-0.44694444)

    print_me("   a.dms_str()", a.dms_str())

    #    a.dms_str(): -26' 48.999984''

    a.set(0, 0, -46.31)

    print_me("   a.dms_str(False)", a.dms_str(False))

    #    a.dms_str(False): 0:0:-46.31

We can use decimals in degrees/minutes. They are converted automatically::

    a.set(0, -46.25, 0.0)

    print_me("   a.dms_str()", a.dms_str())

    #    a.dms_str(): -46' 15.0''

    a.set(0, 0, 0.0)

    print_me("   a.dms_str()", a.dms_str())

    #    a.dms_str(): 0d 0' 0.0''

We can define the angle as in radians. It will be converted to degrees::

    b = Angle(pi, radians=True)

    print_me("b = Angle(pi, radians=True); print(b)", b)

    # b = Angle(pi, radians=True); print(b): 180.0

And we can easily carry out the *degrees to radians* conversion::

    print_me("print(b.rad())", b.rad())

    # print(b.rad()): 3.14159265359

We can also specify the angle as a Right Ascension. Angle can be given as a Right Ascension: Hours, Minutes, Seconds::

    a.set_ra(9, 14, 55.8)

    print_me("   print(a)", a)

    #    print(a): 138.7325

    b = Angle(9, 14, 55.8, ra=True)

    print_me("   print(b)", b)

    #    print(b): 138.7325


We can print the Angle as Right Ascension, as a float and as string::

    a = Angle(138.75)

    print_me("   print(a.get_ra())", a.get_ra())

    #    print(a.get_ra()): 9.25

    print_me("   print(a.ra_str())", a.ra_str())

    #    print(a.ra_str()): 9h 15' 0.0''

    print_me("   print(a.ra_str(False))", a.ra_str(False))

    #    print(a.ra_str(False)): 9:15:0.0


Use the ``to_positive()`` method to get the positive version of an angle::

    a = Angle(-87.32)

    print_me("   print(a.to_positive())", a.to_positive())

    #    print(a.to_positive()): 272.68


Call the ``__repr__()`` method to get a string defining the current object. This string can then be fed to the ``eval()`` function to generate the object::

    print_me("print(b.__repr__())", b.__repr__())

    # print(b.__repr__()): Angle(138.7325)

    c = eval(repr(b))

    print_me("c = eval(repr(b)); print(c)", c)

    # c = eval(repr(b)); print(c): 138.7325

Let's now work with some useful operators and functions::

    print_me("c", c)

    # c: 138.7325

- Negate an angle::

    d = Angle(13, 30)

    print_me("d", d)

    # d: 13.5

    e = -d

    print_me("   e = -d", e)

    #   e = -d: -13.5

- Get the absolute value of an angle::

    e = abs(e)

    print_me("   e = abs(e)", e)

    #    e = abs(e): 13.5

- Module operation on an angle::

    d = Angle(17.0)

    print_me("d", d)

    # d: 17.0

    e = c % d

    print_me("   e = c % d", e)

    #    e = c % d: 2.7325


- Convert the angle to an integer::

    d = Angle(13.95)

    print_me("d", d)

    # d: 13.95

    print_me("   int(d)", int(d))

    #    int(d): 13

    d = Angle(-4.95)

    print_me("d", d)

    # d: -4.95

    print_me("   int(d)", int(d))

    #   int(d): -4

- Convert the angle to a float::

    print_me("   float(d)", float(d))

    #    float(d): -4.95

- Round the angle to a float::

    e = Angle(-4.951648)

    print_me("e", e)

    # e: -4.951648

    print_me("   round(e)", round(e))

    #    round(e): -5.0

    print_me("   round(e, 2)", round(e, 2))

    #    round(e, 2): -4.95

    print_me("   round(e, 3)", round(e, 3))

    #    round(e, 3): -4.952

    print_me("   round(e, 4)", round(e, 4))

    #    round(e, 4): -4.9516

- Comparison operators::

    print_me("   d == e", d == e)

    #    d == e: False

    print_me("   d != e", d != e)

    #    d != e: True

    print_me("   d > e", d > e)

    #    d > e: True

    print_me("   c >= 180.0", c >= 180.0)

    #    c >= 180.0: False

    print_me("   c < 180.0", c < 180.0)

    #    c < 180.0: True

    print_me("   c <= 180.0", c <= 180.0)

    #    c <= 180.0: True

- It is very easy to add Angles to obtain a new Angle::

    e = c + d

    print_me("   c + d", e)

    #    c + d: 133.7825

  We can also directly add a decimal angle::

    e = c + 11.5

    print_me("   c + 11.5", e)

    #    c + 11.5: 150.2325

  Types allowed are int, float and Angle::

    print('e = c + "32.5"')

    # e = c + "32.5"

    try:

        e = c + "32.5"

    except TypeError:

        print("TypeError!: Valid types are int, float, and Angle, not string!")

    # TypeError!: Valid types are int, float, and Angle, not string!


- Subtraction::

    e = c - d

    print_me("   c - d", e)

    #    c - d: 143.6825

- Multiplication::

    c.set(150.0)

    d.set(5.0)

    print_me("c", c)

    # c: 150.0

    print_me("d", d)

    # d: 5.0

    e = c * d

    print_me("   c * d", e)

    #    c * d: 30.0

- Division::

    c.set(150.0)

    d.set(6.0)

    print_me("d", d)

    # d: 6.0

    e = c / d

    print_me("   c / d", e)

    #    c / d: 25.0


  Division by zero is not allowed::

    d.set(0.0)

    print_me("d", d)

    # d: 0.0

    print('e = c / d')

    # e = c / d

    try:

        e = c / d

    except ZeroDivisionError:

        print("ZeroDivisionError!: Division by zero is not allowed!")

    # ZeroDivisionError!: Division by zero is not allowed!

- Power::

    d.set(2.2)

    print_me("d", d)

    # d: 2.2

    e = c ** d

    print_me("   c ** d", e)

    #    c ** d: 91.5733670999


- Accumulative module operation::

    d.set(17.0)

    print_me("d", d)

    # d: 17.0

    e %= d

    print_me("   e %= d", e)

    #    e %= d: 6.57336709993

- Accumulative addition::

    c += d

    print_me("   c += d", c)

    #    c += d: 167.0

- Accumulative subtraction::

    print_me("b", b)

    # b: 138.7325

    c -= b

    print_me("   c -= b", c)

    #    c -= b: 28.2675

- Accumulative multiplication::

    print_me("b", b)

    # b: 138.7325

    c *= b

    print_me("   c *= b", c)

    #    c *= b: 321.62094375

- Accumulative division::

    print_me("b", b)

    # 138.7325

    d.set(6.0)

    print_me("d", d)

    # d: 6.0

    b /= d

    print_me("   b /= d", b)

    #    b /= d: 23.1220833333

- Accumulative power::

    d.set(2.2)

    print_me("d", d)

    # d: 2.2

    c = abs(c)

    print_me("   c = abs(c)", c)

    #    c = abs(c): 321.62094375

    c **= d

    print_me("   c **= d", c)

    #    c **= d: 254.307104203


The same operations, but by the right side::

    e = 3.5 + b

    print_me("   e = 3.5 + b", e)

    #    e = 3.5 + b: 26.6220833333

    e = 3.5 - b

    print_me("   e = 3.5 - b", e)

    #    e = 3.5 - b: -19.6220833333

    e = 3.5 * b

    print_me("   e = 3.5 * b", e)

    #    e = 3.5 * b: 80.9272916667

    e = 3.5 / b

    print_me("   e = 3.5 / b", e)

    #    e = 3.5 / b: 0.151370443119

    e = 3.5 ** b

    print_me("   e = 3.5 ** b", e)

    #    e = 3.5 ** b: 260.783691406

