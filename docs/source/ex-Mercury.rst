Mercury examples
****************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Mercury.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 287.4887

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: -6.0086

    print_me("Radius vector", r)

    # Radius vector: 0.45113

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Mercury.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 16h 33' 59.3''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: -20d 53' 31.6''

    print_me("Elongation", elon.dms_str(n_dec=1))

    # Elongation: 18d 24' 29.8''

Print mean orbital elements for Mercury at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Mercury.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 203.494701

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 0.38709831

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.2056451

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 7.006171

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 49.10765

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: 29.367732

Compute the time of the inferior conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conjunction = Mercury.inferior_conjunction(epoch)

    y, m, d = conjunction.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Inferior conjunction date", date)

    # Inferior conjunction date: 1993/11/6.1449

Compute the time of the superior conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conjunction = Mercury.superior_conjunction(epoch)

    y, m, d = conjunction.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Superior conjunction date", date)

    # Superior conjunction date: 1993/8/29.3301

Compute the time and angle of the western elongation close to 1993/11/1::

    epoch = Epoch(1993, 11, 1.0)

    time, elongation = Mercury.western_elongation(epoch)

    y, m, d = time.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Western elongation date", date)

    # Western elongation date: 1993/11/22.6386

    elong = round(elongation, 4)

    print_me("Maximum western elongation angle", elong)

    # Maximum western elongation angle: 19.7506

Compute the time and angle of the eastern elongation close to 1990/8/1::

    epoch = Epoch(1990, 8, 1.0)

    time, elongation = Mercury.eastern_elongation(epoch)

    y, m, d = time.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Eastern elongation date", date)

    # Eastern elongation date: 1990/8/11.8514

    elong = round(elongation, 4)

    print_me("Maximum eastern elongation angle", elong)

    # Maximum eastern elongation angle: 27.4201

Compute the time of the station in longitude #1 close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    sta1 = Mercury.station_longitude_1(epoch)

    y, m, d = sta1.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #1", date)

    # Date of station in longitude #1: 1993/10/25.9358

Compute the time of the station in longitude #2 close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    sta2 = Mercury.station_longitude_2(epoch)

    y, m, d = sta2.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #2", date)

    # Date of station in longitude #2: 1993/11/15.0724

Find the epoch of the Perihelion closer to 2000/01/01::

    epoch = Epoch(2000, 1, 1.0)

    e = Mercury.perihelion_aphelion(epoch)

    y, m, d, h, mi, s = e.get_full_date()

    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'

    print_me("The Perihelion closest to 2000/1/1 happened on", peri)

    # The Perihelion closest to 2000/1/1 happened on: 2000/2/15 at 18 hours

Compute the time of passage through an ascending node::

    epoch = Epoch(2019, 1, 1)

    time, r = Mercury.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 2018/11/24.7

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 0.3143
