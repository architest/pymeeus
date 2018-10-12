Venus examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(1992, 12, 20.0)

    lon, lat, r = Venus.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", round(lon.to_positive(), 5))

    # Geometric Heliocentric Longitude: 26.11428

    print_me("Geometric Heliocentric Latitude", round(lat, 4))

    # Geometric Heliocentric Latitude: -2.6207

    print_me("Radius vector", round(r, 6))

    # Radius vector: 0.724603
