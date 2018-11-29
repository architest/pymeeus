Saturn examples
***************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

We can compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(2018, 10, 27.0)

    lon, lat, r = Saturn.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 279.5108

    print_me("Geometric Heliocentric Latitude", lat)

    # Geometric Heliocentric Latitude: 0.6141

    print_me("Radius vector", r)

    # Radius vector: 10.06266

Print mean orbital elements for Saturn at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Saturn.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 131.196871

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 9.55490779

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0553209

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 2.486426

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 114.23974

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: -19.896331

Compute the time of the conjunction close to 1993/10/1::

    epoch = Epoch(1993, 10, 1.0)

    conj = Saturn.conjunction(epoch)

    y, m, d = conj.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Conjunction date", date)

    # Conjunction date: 1994/2/21.7347

Compute the time of the opposition close to -6/9/1::

    epoch = Epoch(-6, 9, 1.0)

    oppo = Saturn.opposition(epoch)

    y, m, d = oppo.get_date()

    d = round(d, 4)

    date = "{}/{}/{}".format(y, m, d)

    print_me("Opposition date", date)

    # Opposition date: -6/9/14.3709
