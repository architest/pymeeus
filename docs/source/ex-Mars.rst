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

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Mars.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 7h 48' 35.4''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: 24d 35' 33.9''

    print_me("Elongation", elon.dms_str(n_dec=1))

    # Elongation: 153d 35' 1.6''

Print mean orbital elements for Mars at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Mars.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 288.855211

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 1.52367934

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0934599

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 1.849338

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 50.06365

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: 287.202108

Compute the time of the conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conj = Mars.conjunction(epoch)

    y, m, d = conj.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Conjunction date", date)

    # Conjunction date: 1993/12/27.0898

Compute the time of the opposition close to 2729/10/1::

    epoch = Epoch(2729, 10, 1.0)

    oppo = Mars.opposition(epoch)

    y, m, d = oppo.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Opposition date", date)

    # Opposition date: 2729/9/9.1412

Compute the time of the station in longitude #1 close to 1997/3/1::

    epoch = Epoch(1997, 3, 1.0)

    sta1 = Mars.station_longitude_1(epoch)

    y, m, d = sta1.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #1", date)

    # Date of station in longitude #1: 1997/2/6.033

Compute the time of the station in longitude #2 close to 1997/3/1::

    epoch = Epoch(1997, 3, 1.0)

    sta2 = Mars.station_longitude_2(epoch)

    y, m, d = sta2.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #2", date)

    # Date of station in longitude #2: 1997/4/27.7553

Find the epoch of the Aphelion closer to 2032/1/1::

    epoch = Epoch(2032, 1, 1.0)

    e = Mars.perihelion_aphelion(epoch, perihelion=False)

    y, m, d, h, mi, s = e.get_full_date()

    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'

    print_me("The Aphelion closest to 2032/1/1 will happen on", peri)

    # The Aphelion closest to 2032/1/1 will happen on: 2032/10/24 at 22 hours

Compute the time of passage through an ascending node::

    epoch = Epoch(2019, 1, 1)

    time, r = Mars.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 2019/1/15.2

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 1.4709
