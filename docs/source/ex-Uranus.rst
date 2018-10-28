Uranus examples
***************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Uranus.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 30.5888

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: -0.5315

    print_me("Radius vector", r)

    # Radius vector: 19.86964
