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

Compute the approximate illuminated fraction of the Moon's disk::

    epoch = Epoch(1992, 4, 12.0)

    k = Moon.illuminated_fraction_disk(epoch)

    print_me("Approximate illuminated fraction of Moon's disk", round(k, 2))

    # Approximate illuminated fraction of Moon's disk: 0.68

Compute the position angle of the bright limb of the Moon::

    epoch = Epoch(1992, 4, 12.0)

    xi = Moon.position_bright_limb(epoch)

    print_me("Position angle of the bright limb of the Moon", round(xi, 1))

    # Position angle of the bright limb of the Moon: 285.0

Calculate the instant of a New Moon::

    epoch = Epoch(1977, 2, 15.0)

    new_moon = Moon.moon_phase(epoch, target="new")

    y, m, d, h, mi, s = new_moon.get_full_date()

    print("New Moon: {}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s)))

    # New Moon: 1977/2/18 3:37:42

Calculate the time of a Last Quarter::

    epoch = Epoch(2044, 1, 15.0)

    new_moon = Moon.moon_phase(epoch, target="last")

    y, m, d, h, mi, s = new_moon.get_full_date()

    print("Last Quarter: {}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s)))

    # Last Quarter: 2044/1/21 23:48:17

Compute the time and parallax of apogee::

    epoch = Epoch(1988, 10, 1.0)

    apogee, parallax = Moon.moon_perigee_apogee(epoch, target="apogee")

    y, m, d, h, mi, s = apogee.get_full_date()

    print("Apogee epoch: {}/{}/{} {}:{}".format(y, m, d, h, mi))

    # Apogee epoch: 1988/10/7 20:30

    print_me("Equatorial horizontal parallax", parallax.dms_str(n_dec=3))

    # Equatorial horizontal parallax: 54' 0.679''

Compute the time of passage by the ascending node::

    epoch = Epoch(1987, 5, 15.0)

    passage = Moon.moon_passage_nodes(epoch, target="ascending")

    y, m, d, h, mi, s = passage.get_full_date()

    mi += s/60.0

    print("Passage by the ascending node: {}/{}/{} {}:{}".format(y, m, d, h, round(mi)))

    # Passage by the ascending node: 1987/5/23 6:26

Compute the epoch and amplitude of maximum southern declination::

    epoch = Epoch(2049, 4, 15.0)

    epo, dec = Moon.moon_maximum_declination(epoch, target='southern')

    y, m, d, h, mi, s = epo.get_full_date()

    print("Epoch of maximum declination: {}/{}/{} {}:{}".format(y, m, d, h, mi))

    # Epoch of maximum declination: 2049/4/21 14:0

    print_me("Amplitude of maximum declination", dec.dms_str(n_dec=0))

    # Amplitude of maximum declination: -22d 8' 18.0''

Compute the librations of the Moon::

    epoch = Epoch(1992, 4, 12.0)

    lopt, bopt, lphys, bphys, ltot, btot = Moon.moon_librations(epoch)

    print_me("Optical libration in longitude", round(lopt, 3))

    # Optical libration in longitude: -1.206

    print_me("Optical libration in latitude", round(bopt, 3))

    # Optical libration in latitude: 4.194

    print_me("Physical libration in longitude", round(lphys, 3))

    # Physical libration in longitude: -0.025

    print_me("Physical libration in latitude", round(bphys, 3))

    # Physical libration in latitude: 0.006

    print_me("Total libration in longitude", round(lphys, 2))

    # Total libration in longitude: -1.23

    print_me("Total libration in latitude", round(bphys, 3))

    # Total libration in latitude: 4.2

Let's calculate the position angle of the Moon's axis of rotation::

    epoch = Epoch(1992, 4, 12.0)

    p = Moon.moon_position_angle_axis(epoch)

    print_me("Position angle of Moon's axis of rotation", round(p, 2))

    # Position angle of Moon's axis of rotation: 15.08
