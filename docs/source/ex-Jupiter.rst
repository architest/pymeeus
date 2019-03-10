Jupiter examples
****************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Jupiter.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 241.5873

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: 0.8216

    print_me("Radius vector", r)

    # Radius vector: 5.36848

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Jupiter.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 12h 47' 9.6''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: -3d 41' 55.3''

    print_me("Elongation", elon.dms_str(n_dec=1))

    # Elongation: 76d 2' 26.0''

Print mean orbital elements for Jupiter at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Jupiter.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 222.433723

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 5.20260333

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0486046

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 1.29967

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 101.13309

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: -85.745532

Compute the time of the conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conj = Jupiter.conjunction(epoch)

    y, m, d = conj.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Conjunction date", date)

    # Conjunction date: 1993/10/18.3341

Compute the time of the opposition close to -6/9/1::

    epoch = Epoch(-6, 9, 1.0)

    oppo = Jupiter.opposition(epoch)

    y, m, d = oppo.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Opposition date", date)

    # Opposition date: -6/9/15.2865

Compute the time of the station in longitude #1 close to 2018/11/1::

    epoch = Epoch(2018, 11, 1.0)

    sta1 = Jupiter.station_longitude_1(epoch)

    y, m, d = sta1.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #1", date)

    # Date of station in longitude #1: 2018/3/9.1288

Compute the time of the station in longitude #2 close to 2018/11/1::

    epoch = Epoch(2018, 11, 1.0)

    sta2 = Jupiter.station_longitude_2(epoch)

    y, m, d = sta2.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #2", date)

    # Date of station in longitude #2: 2018/7/10.6679

Find the epoch of the Aphelion closer to 1981/6/1::

    epoch = Epoch(1981, 6, 1.0)

    e = Jupiter.perihelion_aphelion(epoch, perihelion=False)

    y, m, d, h, mi, s = e.get_full_date()

    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'

    print_me("The Aphelion closest to 1981/6/1 will happen on", peri)

    # The Aphelion closest to 1981/6/1 will happen on: 1981/7/28 at 6 hours

Compute the time of passage through an ascending node::

    epoch = Epoch(2019, 1, 1)

    time, r = Jupiter.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 2025/9/15.6

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 5.1729
