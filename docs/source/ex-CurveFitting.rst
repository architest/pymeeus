CurveFitting examples
*********************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Now let's work with the CurveFitting class. First, create a CurveFitting object::

    cf1 = CurveFitting([73.0, 38.0, 35.0, 42.0, 78.0, 68.0, 74.0, 42.0, 52.0,

                        54.0, 39.0, 61.0, 42.0, 49.0, 50.0, 62.0, 44.0, 39.0,

                        43.0, 54.0, 44.0, 37.0],

                       [90.4, 125.3, 161.8, 143.4, 52.5, 50.8, 71.5, 152.8,

                        131.3, 98.5, 144.8, 78.1, 89.5, 63.9, 112.1, 82.0,

                        119.8, 161.2, 208.4, 111.6, 167.1, 162.1])


Let's use ``linear_fitting()``::

    a, b = cf1.linear_fitting()

    print("   a = {}\tb = {}".format(round(a, 2), round(b, 2)))

    #    a = -2.49	b = 244.18

Use the copy constructor::

    cf2 = CurveFitting(cf1)

    a, b = cf2.linear_fitting()

    print("   a = {}\tb = {}".format(round(a, 2), round(b, 2)))

    #    a = -2.49	b = 244.18

Get the number of value pairs internally stored::

    print_me("Number of value pairs inside 'cf2'", len(cf2))

    # Number of value pairs inside 'cf2': 22

Compute the correlation coefficient::

    r = cf1.correlation_coeff()

    print_me("   r", round(r, 3))

    #    r: -0.767

Define a new CurveFitting object::

    cf2 = CurveFitting([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],

                       [-9.372, -3.821, 0.291, 3.730, 5.822, 8.324, 9.083,

                        6.957, 7.006, 0.365, -1.722])

Now use ``quadratic_fitting()``::

    a, b, c = cf2.quadratic_fitting()

    # Original curve: y = -2.0*x*x + 3.5*x + 7.0 + noise

    print("   a = {}\tb = {}\tc = {}".format(round(a, 2), round(b, 2), round(c, 2)))

    #    a = -2.22	b = 3.76	c = 6.64

Yet another CurveFitting object::

    cf4 = CurveFitting([3, 20, 34, 50, 75, 88, 111, 129, 143, 160, 183, 200,

                        218, 230, 248, 269, 290, 303, 320, 344],

                       [0.0433, 0.2532, 0.3386, 0.3560, 0.4983, 0.7577, 1.4585,

                        1.8628, 1.8264, 1.2431, -0.2043, -1.2431, -1.8422,

                        -1.8726, -1.4889, -0.8372, -0.4377, -0.3640, -0.3508,

                        -0.2126])

Let's define the three functions to be used for fitting::

    def sin1(x): return sin(radians(x))

    def sin2(x): return sin(radians(2.0*x))

    def sin3(x): return sin(radians(3.0*x))

Use ``general_fitting()`` here::

    a, b, c = cf4.general_fitting(sin1, sin2, sin3)

    # General fitting with f0 = sin(x), f1 = sin(2*x), f2 = sin(3*x)

    print("   a = {}\tb = {}\tc = {}".format(round(a, 2), round(b, 2), round(c, 2)))

    #    a = 1.2	b = -0.77	c = 0.39

A final example::

    cf5 = CurveFitting([0, 1.2, 1.4, 1.7, 2.1, 2.2])

    a, b, c = cf5.general_fitting(sqrt)

    # General fitting with f0 = sqrt(x), f1 = 0.0 and f2 = 0.0

    print("   a = {}\tb = {}\t\tc = {}".format(round(a, 3), round(b, 3), round(c, 3)))

    #   a = 1.016	b = 0.0		c = 0.0
