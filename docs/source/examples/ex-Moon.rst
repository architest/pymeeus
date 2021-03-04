Moon examples
*************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Let's now compute the Moon geocentric ecliptical position for a given epoch::

    epoch = Epoch(1992, 4, 12.0)

    Lambda, Beta, Delta, ppi = Moon.geocentric_ecliptical_pos(epoch)

    print_me("Longitude (Lambda)", round(Lambda, 6))

    # Longitude (Lambda): 133.162655

    print_me("Latitude (Beta)", round(Beta, 6))

    # Latitude (Beta): -3.229126

    print_me("Distance (Delta)", round(Delta, 1))

    # Distance (Delta): 368409.7

    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))

    # Equatorial horizontal parallax (Pi): 0.99199

Now let's compute the apparent ecliptical position::

    epoch = Epoch(1992, 4, 12.0)

    Lambda, Beta, Delta, ppi = Moon.apparent_ecliptical_pos(epoch)

    print_me("Longitude (Lambda)", round(Lambda, 6))

    # Longitude (Lambda): 133.167264

    print_me("Latitude (Beta)", round(Beta, 6))

    # Latitude (Beta): -3.229126

    print_me("Distance (Delta)", round(Delta, 1))

    # Distance (Delta): 368409.7

    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))

    # Equatorial horizontal parallax (Pi): 0.99199

Get the apparent equatorial position::

    epoch = Epoch(1992, 4, 12.0)

    ra, dec, Delta, ppi = Moon.apparent_equatorial_pos(epoch)

    print_me("Right Ascension (ra)", round(ra, 6))

    # Right Ascension (ra): 134.688469

    print_me("Declination (dec)", round(dec, 6))

    # Declination (dec): 13.768367

    print_me("Distance (Delta)", round(Delta, 1))

    # Distance (Delta): 368409.7

    print_me("Equatorial horizontal parallax (Pi)", round(ppi, 6))

    # Equatorial horizontal parallax (Pi): 0.99199

Compute the longitude of the Moon's mean ascending node::

    epoch = Epoch(1913, 5, 27.0)

    Omega = Moon.longitude_mean_ascending_node(epoch)

    print_me("Longitude of the mean ascending node", round(Omega, 1))

    # Longitude of the mean ascending node: 0.0

    epoch = Epoch(1959, 12, 7.0)

    Omega = Moon.longitude_mean_ascending_node(epoch)

    print_me("Longitude of the mean ascending node", round(Omega, 1))

    # Longitude of the mean ascending node: 180.0

Get the longitude of the Moonś true ascending node::

    epoch = Epoch(1913, 5, 27.0)

    Omega = Moon.longitude_true_ascending_node(epoch)

    print_me("Longitude of the true ascending node", round(Omega, 4))

    # Longitude of the true ascending node: 0.8763

Compute the longitude of the Moon's mean perigee::

    epoch = Epoch(2021, 3, 5.0)

    Pi = Moon.longitude_mean_perigee(epoch)

    print_me("Longitude of the mean perigee", round(Pi, 5))

    # Longitude of the mean perigee: 224.89194

