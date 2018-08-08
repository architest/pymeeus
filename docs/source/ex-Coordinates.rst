Coordinates examples
********************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

It follows a series of important parameters related to the angle between Earth's rotation axis and the ecliptic.

- The mean angle between Earth rotation axis and ecliptic axis is the **mean obliquity**::

    e0 = Earth.mean_obliquity(1987, 4, 10)

    print_me("Mean obliquity for 1987/4/10", e0.dms_str())

    # Mean obliquity for 1987/4/10: 23d 26' 27.4066466278''

- If we take into account the nutation effect on the obliquity, we get the **true obliquity**::

    epsilon = Earth.true_obliquity(1987, 4, 10)

    print_me("True obliquity for 1987/4/10", epsilon.dms_str())

    # True obliquity for 1987/4/10: 23d 26' 36.8491882378''

    epsilon = Earth.true_obliquity(2018, 7, 29)

    print_me("True obliquity for 2018/7/29", epsilon.dms_str())

    # True obliquity for 2018/7/29: 23d 26' 7.21570241139''

- The nutation effect is separated in two components: One parallel to the ecliptic (nutation in longitude) and other perpendicular to the ecliptic (nutation in obliquity)::

    dpsi = Earth.nutation_longitude(1987, 4, 10)

    print_me("Nutation in longitude for 1987/4/10", dpsi.dms_str(n_dec=3))

    # Nutation in longitude for 1987/4/10: -3.788''

    depsilon = Earth.nutation_obliquity(1987, 4, 10)

    print_me("Nutation in obliquity for 1987/4/10", depsilon.dms_str(n_dec=3))

    # Nutation in obliquity for 1987/4/10: 9.443''

- We can compute the effects of precession on the equatorial coordinates of a given star, taking also into account its proper motion::

    start_epoch = JDE2000

    final_epoch = Epoch(2028, 11, 13.19)

    alpha0 = Angle(2, 44, 11.986, ra=True)

    delta0 = Angle(49, 13, 42.48)

    print_me("Initial right ascension", alpha0.ra_str(n_dec=3))

    # Initial right ascension: 2h 44' 11.986''

    print_me("Initial declination", delta0.dms_str(n_dec=2))

    # Initial declination: 49d 13' 42.48''

    pm_ra = Angle(0, 0, 0.03425, ra=True)

    pm_dec = Angle(0, 0, -0.0895)

    alpha, delta = Earth.precession_equatorial(start_epoch, final_epoch,

                                               alpha0, delta0, pm_ra, pm_dec)

    print_me("Final right ascension", alpha.ra_str(n_dec=3))

    # Final right ascension: 2h 46' 11.331''

    print_me("Final declination", delta.dms_str(n_dec=2))

    # Final declination: 49d 20' 54.54''

Something similar can also be done with the ecliptical coordinates::

    start_epoch = JDE2000

    final_epoch = Epoch(-214, 6, 30.0)

    lon0 = Angle(149.48194)

    lat0 = Angle(1.76549)

    print_me("Initial ecliptical longitude", round(lon0(), 5))

    # Initial ecliptical longitude: 149.48194

    print_me("Initial ecliptical latitude", round(lat0(), 5))

    # Initial ecliptical latitude: 1.76549

    lon, lat = Earth.precession_ecliptical(start_epoch, final_epoch, lon0, lat0)

    print_me("Final ecliptical longitude", round(lon(), 3))

    # Final ecliptical longitude: 118.704

    print_me("Final ecliptical latitude", round(lat(), 3))

    # Final ecliptical latitude: 1.615

Additionally, module ``Coordinates`` provides a function to compute the true movement of a star through the sky relative to the Sun::

    ra = Angle(6, 45, 8.871, ra=True)

    dec = Angle(-16.716108)

    pm_ra = Angle(0, 0, -0.03847, ra=True)

    pm_dec = Angle(0, 0, -1.2053)

    dist = 2.64

    vel = -7.6

    alpha, delta = motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -4000.0)

    print(alpha.ra_str(False, 2))

    # 6:47:39.91

    print(delta.dms_str(False, 1))

    # -15:23:30.6
