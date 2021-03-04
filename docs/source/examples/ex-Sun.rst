Sun examples
************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

It is possible to compute an approximation of the Sun's **true** ecliptical longitude::

    epoch = Epoch(1992, 10, 13)

    true_lon, r = Sun.true_longitude_coarse(epoch)

    print_me("Sun's approximate true longitude", true_lon.dms_str(n_dec=0))

    # Sun's approximate true longitude: 199d 54' 36.0''

    print_me("Sun's radius vector", round(r, 5))

    # Sun's radius vector: 0.99766

Now let's compute the Sun's approximate **apparent** ecliptical longitude::

    epoch = Epoch(1992, 10, 13)

    app_lon, r = Sun.apparent_longitude_coarse(epoch)

    print_me("Sun's approximate apparent longitude", app_lon.dms_str(n_dec=0))

    # Sun's approximate apparent longitude: 199d 54' 32.0''

And now is the turn for the **apparent** right ascension and declination::

    epoch = Epoch(1992, 10, 13)

    ra, delta, r = Sun.apparent_rightascension_declination_coarse(epoch)

    print_me("Sun's apparent right ascension", ra.ra_str(n_dec=1))

    # Sun's apparent right ascension: 13h 13' 31.4''

    print_me("Sun's apparent declination", delta.dms_str(n_dec=0))

    # Sun's apparent declination: -7d 47' 6.0''

Now, let's compute Sun's true (**geometric**) position again, but more accurately::

    epoch = Epoch(1992, 10, 13.0)

    l, b, r = Sun.geometric_geocentric_position(epoch, toFK5=False)

    print_me("Geometric Geocentric Longitude", round(l.to_positive(), 6))

    # Geometric Geocentric Longitude: 199.907297

    print_me("Geometric Geocentric Latitude", b.dms_str(n_dec=3))

    # Geometric Geocentric Latitude: 0.744''

    print_me("Radius vector", round(r, 8))

    # Radius vector: 0.99760852

Compute Sun's **apparent** postion accurately::

    epoch = Epoch(1992, 10, 13.0)

    l, b, r = Sun.apparent_geocentric_position(epoch)

    print_me("Apparent Geocentric Longitude", l.to_positive().dms_str(n_dec=3))

    # Apparent Geocentric Longitude: 199d 54' 21.548''

    print_me("Apparent Geocentric Latitude", b.dms_str(n_dec=3))

    # Apparent Geocentric Latitude; 0.721''

    print_me("Radius vector", round(r, 8))

    # Radius vector: 0.99760852

We can compute rectangular coordinates referred to mean equinox of date::

    epoch = Epoch(1992, 10, 13.0)

    x, y, z = Sun.rectangular_coordinates_mean_equinox(epoch)

    print_me("X", round(x, 7))

    # X: -0.9379963

    print_me("Y", round(y, 6))

    # Y: -0.311654

    print_me("Z", round(z, 7))

    # Z: -0.1351207

Now, compute rectangular coordinates w.r.t. standard equinox J2000.0::

    epoch = Epoch(1992, 10, 13.0)

    x, y, z = Sun.rectangular_coordinates_j2000(epoch)

    print_me("X", round(x, 8))

    # X: -0.93740485

    print_me("Y", round(y, 8))

    # Y: -0.3131474

    print_me("Z", round(z, 8))

    # Z: -0.12456646

Compute rectangular coordinates w.r.t. mean equinox of B1950.0::

    epoch = Epoch(1992, 10, 13.0)

    x, y, z = Sun.rectangular_coordinates_b1950(epoch)

    print_me("X", round(x, 8))

    # X: -0.94149557

    print_me("Y", round(y, 8))

    # Y: -0.30259922

    print_me("Z", round(z, 8))

    # Z: -0.11578695

And compute rectangular coordinates w.r.t. an arbitrary mean equinox::

    epoch = Epoch(1992, 10, 13.0)

    e_equinox = Epoch(2467616.0)

    x, y, z = Sun.rectangular_coordinates_equinox(epoch, e_equinox)

    print_me("X", round(x, 8))

    # X: -0.93373777

    print_me("Y", round(y, 8))

    # Y: -0.32235109

    print_me("Z", round(z, 8))

    # Z: -0.12856709

It is possible to compute the date of equinoxes and solstices::

    epoch = Sun.get_equinox_solstice(1962, target="summer")

    y, m, d, h, mi, s = epoch.get_full_date()

    print("The summer solstice of 1962:")

    print("{}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))

    # 1962/6/21 21:24:42.0

The equation of time, i.e., the difference between apparent and mean time, can be easily computed::

    epoch = Epoch(1992, 10, 13.0)

    m, s = Sun.equation_of_time(epoch)

    print("Equation of time difference: {} min {} secs".format(m, round(s, 1)))

    # Equation of time difference: 13 min 42.6 secs

Compute the ephemeris of physical observations of the Sun using Carrington's formulas::

    epoch = Epoch(1992, 10, 13)

    p, b0, l0 = Sun.ephemeris_physical_observations(epoch)

    print("Ephemeris of physical observations of the Sun:")

    print_me("P ", round(p, 2))

    # P : 26.27

    print_me("B0", round(b0, 2))

    # B0: 5.99

    print_me("L0", round(l0, 2))

    # L0: 238.63

Get the epoch when the Carrington's synodic rotation No. 'number' starts::

    epoch = Sun.beginning_synodic_rotation(1699)

    print_me("Epoch for Carrington's synodic rotation No. 1699", round(epoch(), 3))

    # Epoch for Carrington's synodic rotation No. 1699: 2444480.723
