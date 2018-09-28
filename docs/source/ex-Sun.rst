Sun examples
************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

It is possible to compute an approximation of the Sun's **true** ecliptical longitude::

    epoch = Epoch(1992, 10, 13)

    true_lon, r = sun_true_longitude_coarse(epoch)

    print_me("Sun's approximate true longitude", true_lon.dms_str(n_dec=0))

    # Sun's approximate true longitude: 199d 54' 36.0''

    print_me("Sun's radius vector", round(r, 5))

    # Sun's radius vector: 0.99766

Now let's compute the Sun's approximate **apparent** ecliptical longitude::

    epoch = Epoch(1992, 10, 13)

    app_lon, r = sun_apparent_longitude_coarse(epoch)

    print_me("Sun's approximate apparent longitude", app_lon.dms_str(n_dec=0))

    # Sun's approximate apparent longitude: 199d 54' 32.0''

And now is the turn for the **apparent** right ascension and declination::

    epoch = Epoch(1992, 10, 13)

    ra, delta, r = sun_apparent_rightascension_declination_coarse(epoch)

    print_me("Sun's apparent right ascension", ra.ra_str(n_dec=1))

    # Sun's apparent right ascension: 13h 13' 31.4''

    print_me("Sun's apparent declination", delta.dms_str(n_dec=0))

    # Sun's apparent declination: -7d 47' 6.0''
