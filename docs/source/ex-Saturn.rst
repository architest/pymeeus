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
