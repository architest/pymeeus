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
