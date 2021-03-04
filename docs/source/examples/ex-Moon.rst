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

