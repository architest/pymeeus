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

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Uranus.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 19h 13' 48.7''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: -22d 46' 13.0''

    print_me("Elongation", elon.dms_str(n_dec=1))

    # Elongation: 18d 44' 18.7''

Print mean orbital elements for Uranus at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Uranus.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 235.517526

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 19.21844604

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0463634

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 0.77372

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 74.34776

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: 99.630865

Compute the time of the conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conj = Uranus.conjunction(epoch)

    y, m, d = conj.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Conjunction date", date)

    # Conjunction date: 1994/1/12.7365

Compute the time of the opposition close to 1780/12/1::

    epoch = Epoch(1780, 12, 1.0)

    oppo = Uranus.opposition(epoch)

    y, m, d = oppo.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Opposition date", date)

    # Opposition date: 1780/12/17.5998

Find the epoch of the Perihelion closer to 1780/1/1::

    epoch = Epoch(1780, 1, 1.0)

    e = Uranus.perihelion_aphelion(epoch)

    y, m, d = e.get_date()

    peri = str(y) + '/' + str(m) + '/' + str(int(d))

    print_me("The Perihelion closest to 1780/1/1 happened on", peri)

    # The Perihelion closest to 1780/1/1 happened on: 1798/2/26

Compute the time of passage through an ascending node::

    epoch = Epoch(2019, 1, 1)

    time, r = Uranus.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 2028/8/23.2

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 19.3201
