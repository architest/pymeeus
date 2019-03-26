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

Compute the geocentric position for 1992/12/20::

    epoch = Epoch(1992, 12, 20.0)

    ra, dec, elon = Venus.geocentric_position(epoch)

    print_me("Right ascension", ra.ra_str(n_dec=1))

    # Right ascension: 21h 4' 41.5''

    print_me("Declination", dec.dms_str(n_dec=1))

    # Declination: -18d 53' 16.8''

    print_me("Elongation", elon.dms_str(n_dec=1))

    # Elongation: 44d 46' 8.9''

Print mean orbital elements for Venus at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Venus.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 338.646306

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 0.72332982

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0067407

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 3.395319

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 77.27012

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: 55.211257

Compute the time of the inferior conjunction close to 1882/12/1.0::

    epoch = Epoch(1882, 12, 1.0)

    conjunction = Venus.inferior_conjunction(epoch)

    y, m, d = conjunction.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Inferior conjunction date", date)

    # Inferior conjunction date: 1882/12/6.6912

Compute the time of the superior conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conjunction = Venus.superior_conjunction(epoch)

    y, m, d = conjunction.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Superior conjunction date", date)

    # Superior conjunction date: 1994/1/17.0465

Compute the time and angle of the western elongation close to 2019/1/1::

    epoch = Epoch(2019, 1, 1.0)

    time, elongation = Venus.western_elongation(epoch)

    y, m, d = time.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Western elongation date", date)

    # Western elongation date: 2019/1/6.1895

    elong = round(elongation, 4)

    print_me("Maximum western elongation angle", elong)

    # Maximum western elongation angle: 46.9571

Compute the time and angle of the eastern elongation close to 2019/10/1::

    epoch = Epoch(2019, 10, 1.0)

    time, elongation = Venus.eastern_elongation(epoch)

    y, m, d = time.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Eastern elongation date", date)

    # Eastern elongation date: 2020/3/24.9179

    elong = round(elongation, 4)

    print_me("Maximum eastern elongation angle", elong)

    # Maximum eastern elongation angle: 46.078

Compute the time of the station in longitude #1 close to 2018/12/1::

    epoch = Epoch(2018, 12, 1.0)

    sta1 = Venus.station_longitude_1(epoch)

    y, m, d = sta1.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #1", date)

    # Date of station in longitude #1: 2018/10/5.7908

Compute the time of the station in longitude #2 close to 2018/12/1::

    epoch = Epoch(2018, 12, 1.0)

    sta2 = Venus.station_longitude_2(epoch)

    y, m, d = sta2.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Date of station in longitude #2", date)

    # Date of station in longitude #2: 2018/11/16.439

Find the epoch of the Perihelion closer to 1978/10/15::

    epoch = Epoch(1978, 10, 15.0)

    e = Venus.perihelion_aphelion(epoch)

    y, m, d, h, mi, s = e.get_full_date()

    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'

    print_me("The Perihelion closest to 1978/10/15 happened on", peri)

    # The Perihelion closest to 1978/10/15 happened on: 1978/12/31 at 4 hours

Compute the time of passage through an ascending node::

    epoch = Epoch(1979, 1, 1)

    time, r = Venus.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 1978/11/27.4

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 0.7205

Compute the (approximate) illuminated fraction of Venus disk for an Epoch::

    epoch = Epoch(1992, 12, 20)

    k = Venus.illuminated_fraction(epoch)

    print_me("Approximate illuminated fraction of Venus", round(k, 2))

    # Approximate illuminated fraction of Venus: 0.64

Compute the magnitude of Venus::

    sun_dist = 0.724604

    earth_dist = 0.910947

    phase_angle = Angle(72.96)

    m = Venus.magnitude(sun_dist, earth_dist, phase_angle)

    print_me("Venus' magnitude", round(m, 1))

    # Venus' magnitude: -3.8
