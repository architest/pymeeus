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

This module ``Coordinates`` also provides a series of functions to convert between equatorial, ecliptical, horizontal and galactic coordinates.

- Equatorial to ecliptical coordinates::

    ra = Angle(7, 45, 18.946, ra=True)

    dec = Angle(28, 1, 34.26)

    epsilon = Angle(23.4392911)

    lon, lat = equatorial2ecliptical(ra, dec, epsilon)

    print_me("Equatorial to ecliptical. Longitude", round(lon(), 5))

    # Equatorial to ecliptical. Longitude: 113.21563

    print_me("Equatorial to ecliptical. Latitude", round(lat(), 5))

    # Equatorial to ecliptical. Latitude: 6.68417

- Ecliptical to equatorial coordinates::

    lon = Angle(113.21563)

    lat = Angle(6.68417)

    epsilon = Angle(23.4392911)

    ra, dec = ecliptical2equatorial(lon, lat, epsilon)

    print_me("Ecliptical to equatorial. Right ascension", ra.ra_str(n_dec=3))

    # Ecliptical to equatorial. Right ascension: 7h 45' 18.946''

    print_me("Ecliptical to equatorial. Declination", dec.dms_str(n_dec=2))

    # Ecliptical to equatorial. Declination: 28d 1' 34.26''

- Equatorial to horizontal coordinates::

    lon = Angle(77, 3, 56)

    lat = Angle(38, 55, 17)

    ra = Angle(23, 9, 16.641, ra=True)

    dec = Angle(-6, 43, 11.61)

    theta0 = Angle(8, 34, 57.0896, ra=True)

    eps = Angle(23, 26, 36.87)

    # Compute correction to convert from mean to apparent sidereal time

    delta = Angle(0, 0, ((-3.868*cos(eps.rad()))/15.0), ra=True)

    theta0 += delta

    h = theta0 - lon - ra

    azi, ele = equatorial2horizontal(h, dec, lat)

    print_me("Equatorial to horizontal: Azimuth", round(azi, 3))

    # Equatorial to horizontal: Azimuth: 68.034

    print_me("Equatorial to horizontal: Elevation", round(ele, 3))

    # Equatorial to horizontal: Elevation: 15.125

- Horizontal to equatorial coordinates::

    azi = Angle(68.0337)

    ele = Angle(15.1249)

    lat = Angle(38, 55, 17)

    h, dec = horizontal2equatorial(azi, ele, lat)

    print_me("Horizontal to equatorial. Hour angle", round(h, 4))

    # Horizontal to equatorial. Hour angle: 64.3521

    print_me("Horizontal to equatorial. Declination", dec.dms_str(n_dec=0))

    # Horizontal to equatorial. Declination: -6d 43' 12.0''

- Equatorial to galactic coordinates::

    ra = Angle(17, 48, 59.74, ra=True)

    dec = Angle(-14, 43, 8.2)

    lon, lat = equatorial2galactic(ra, dec)

    print_me("Equatorial to galactic. Longitude", round(lon, 4))

    # Equatorial to galactic. Longitude: 12.9593

    print_me("Equatorial to galactic. Latitude", round(lat, 4))

    # Equatorial to galactic. Latitude: 6.0463

- Galactic to equatorial coordinates::

    lon = Angle(12.9593)

    lat = Angle(6.0463)

    ra, dec = galactic2equatorial(lon, lat)

    print_me("Galactic to equatorial. Right ascension", ra.ra_str(n_dec=1))

    # Galactic to equatorial. Right ascension: 17h 48' 59.7''

    print_me("Galactic to equatorial. Declination", dec.dms_str(n_dec=0))

    # Galactic to equatorial. Declination: -14d 43' 8.0''

In addition, there is a function to compute the ecliptic longitudes of the two points of the ecliptic which are on the horizon, as well as the angle between the ecliptic and the horizon::

    sidereal_time = Angle(5.0, ra=True)

    lat = Angle(51.0)

    epsilon = Angle(23.44)

    lon1, lon2, i = ecliptic_horizon(sidereal_time, lat, epsilon)

    print_me("Longitude of ecliptic point #1 on the horizon", lon1.dms_str(n_dec=1))

    # Longitude of ecliptic point #1 on the horizon: 169d 21' 29.9''

    print_me("Longitude of ecliptic point #2 on the horizon", lon2.dms_str(n_dec=1))

    # Longitude of ecliptic point #2 on the horizon: 349d 21' 29.9''

    print_me("Angle between the ecliptic and the horizon", round(i, 0))

    # Angle between the ecliptic and the horizon: 62.0

Also, it is possible to compute the angle of the diurnal path of a celestial body relative to the horizon at the time of rising and setting::

    dec = Angle(23.44)

    lat = Angle(40.0)

    j = diurnal_path_horizon(dec, lat)

    print_me("Diurnal path vs. horizon angle at time of rising and setting", j.dms_str(n_dec=1))

    # Diurnal path vs. horizon angle at time of rising and setting: 45d 31' 28.4''
