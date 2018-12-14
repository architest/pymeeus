Pluto examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Pluto.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 232.740711423

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: 14.5878173017

    print_me("Radius vector", r)

    # Radius vector: 29.711110981

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Pluto.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 15h 31' 43.7''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: -4d 27' 28.8''
