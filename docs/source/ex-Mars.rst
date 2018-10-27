Mars examples
*************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Mars.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 2.0015

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: -1.3683

    print_me("Radius vector", r)

    # Radius vector: 1.39306
