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
