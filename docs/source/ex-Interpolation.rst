Interpolation examples
**********************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Declare an Interpolation object::

    i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])

    print(i)

    # X: [1, 2, 3, 4, 5, 6]

    # Y: [2, 4, 6, 8, 10, 12]

.. note::

   a. They are ordered in 'x'
   b. The extra value in 'x' was dropped

Use the copy constructor. We can easily make a copy of an Interpolation object::

    j = Interpolation(i)

    print(j)

    # X: [1, 2, 3, 4, 5, 6]

    # Y: [2, 4, 6, 8, 10, 12]

    j = Interpolation([0.0, 1.0, 3.0], [-1.0, -2.0, 2.0])

    print(j)

    # X: [0.0, 1.0, 3.0]

    # Y: [-1.0, -2.0, 2.0]

    print_me("j(2)", j(2))

    # j(2): -1.0

    print_me("j(0.5)", j(0.5))

    # j(0.5): -1.75

- Test with a value already in the data table::

    print_me("j(1)", j(1))

    # j(1): -2.0

Get the number of interpolation points internally stored::

    print_me("Number or interpolation points in 'j'", len(j))

    # Number or interpolation points in 'j': 3

We can interpolate Angles too::

    k = Interpolation([27.0, 27.5, 28.0, 28.5, 29.0],

                      [Angle(0, 54, 36.125), Angle(0, 54, 24.606),

                       Angle(0, 54, 15.486), Angle(0, 54, 8.694),

                       Angle(0, 54, 4.133)])

    print_me("k(28.278)", Angle(k(28.278)).dms_str())

    # k(28.278): 54' 11.4279073579''

Let's work with a new Interpolation object::

    m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])

    print(m)

    # X: [-1.0, 0.0, 1.0]

    # Y: [-2.0, 3.0, 2.0]

- Get some interpolated values::

    print_me("m(-0.5)", m(-0.5))

    # m(-0.5): 1.25

    print_me("m(0.5)", m(0.5))

    # m(0.5): 3.25

- Get derivatives::

    print_me("m'(-1.0)", m.derivative(-1.0))

    # m'(-1.0): 8.0

    print_me("m'(-0.5)", m.derivative(-0.5))

    # m'(-0.5): 5.0

    print_me("m'(0.0)", m.derivative(0.0))

    # m'(0.0): 2.0

    print_me("m'(0.5)", m.derivative(0.5))

    # m'(0.5): -1.0

    print_me("m'(1.0)", m.derivative(1.0))

    # m'(1.0): -4.0

- Get the root within the interval::

    print_me("m.root()", m.root())

    # m.root(): -0.720759220056

- Get the extremum within the interval::

    print_me("m.minmax()", m.minmax())

    # m.minmax(): 0.333333333333

Let's work now with the interpolation of sine function::

    m = Interpolation([29.43, 30.97, 27.69, 28.11, 31.58, 33.05],

                      [0.4913598528, 0.5145891926, 0.4646875083,

                       0.4711658342, 0.5236885653, 0.5453707057])

    print_me("sin(29.5)\t", m(29.5))

    # sin(29.5)	: 0.492423560118

    print_me("sin(30.0)\t", m(30.0))

    # sin(30.0)	: 0.500000000018

    print_me("sin(30.5)\t", m(30.5))

    # sin(30.5)	: 0.507538362978

Derivatives must be adjusted because degrees were used instead of radians::

    print_me("sin'(29.5)\t", degrees(m.derivative(29.5)))

    # sin'(29.5)	: 0.870355696916

    print_me("sin'(30.0)\t", degrees(m.derivative(30.0)))

    # sin'(30.0)	: 0.866025403791

    print_me("sqrt(3.0)/2.0\t", sqrt(3.0)/2.0)

    # sqrt(3.0)/2.0	: 0.866025403784

    print_me("sin'(30.5)\t", degrees(m.derivative(30.5)))

    # sin'(30.5)	: 0.861629160353
