Minor examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Let's compute the equatorial coordinates of comet Encke::

    a = 2.2091404

    e = 0.8502196

    i = Angle(11.94524)

    omega = Angle(334.75006)

    w = Angle(186.23352)

    t = Epoch(1990, 10, 28.54502)

    epoch = Epoch(1990, 10, 6.0)

    ra, dec = Minor.geocentric_position(a, e, i, omega, w, t, epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 10h 34' 13.7''

    print_me("Declination", dec.dms_str(n_dec=0))

    # Declination: 19d 9' 32.0''

Now compute the heliocentric ecliptical coordinates::

    a = 2.2091404

    e = 0.8502196

    i = Angle(11.94524)

    omega = Angle(334.75006)

    w = Angle(186.23352)

    t = Epoch(1990, 10, 28.54502)

    epoch = Epoch(1990, 10, 6.0)

    lon, lat = Minor.heliocentric_ecliptical_position(a, e, i, omega, w, t, epoch)

    print_me("Heliocentric ecliptical longitude", lon.dms_str(n_dec=1))

    # Heliocentric ecliptical longitude: 66d 51' 57.8''

    print_me("Heliocentric ecliptical latitude", lat.dms_str(n_dec=1))

    # Heliocentric ecliptical latitude: 11d 56' 14.4''
