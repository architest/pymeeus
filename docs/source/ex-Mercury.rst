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
