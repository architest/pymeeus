Earth examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

An important concept are the reference ellipsoids, comprising information about
the Earth global model we are going to use.

A very important reference ellipsoid is **WGS84**, predefined here::

    print_me("WGS84", WGS84)

    # WGS84: 6378137.0:0.00335281066475:7.292115e-05

    # First field is equatorial radius, second field is the flattening, and the

    # third field is the angular rotation velocity, in radians per second

Let's print the semi-minor axis (polar radius)::

    print_me("Polar radius, b", WGS84.b())

    # Polar radius, b: 6356752.31425

And now, let's print the eccentricity of Earth's meridian::

    print_me("Eccentricity, e", WGS84.e())

    # Eccentricity, e: 0.0818191908426

We create an Earth object with a given reference ellipsoid. By default, it is
**WGS84**, but we can use another::

    e = Earth(IAU76)

Print the parameters of reference ellipsoid being used::

    print_me("'e' Earth object parameters", e)

    # 'e' Earth object parameters: 6378140.0:0.0033528131779:7.292114992e-05

Compute the distance to the center of the Earth from a given point at sea
level, and at a certain latitude. It is given as a fraction of equatorial
radius::

    lat = Angle(65, 45, 30.0)               # We can use an Angle for this

    print_me("Distance to Earth's center, from latitude 65d 45' 30''", e.rho(lat))

    # Distance to Earth's center, from latitude 65d 45' 30'': 0.997216343095

Parameters *rho\*sin(lat)* and *rho\*cos(lat)* are useful for different
astronomical applications::

    height = 650.0

    print_me("rho*sin(lat)", e.rho_sinphi(lat, height))

    # rho*sin(lat): 0.908341718779

    print_me("rho*cos(lat)", e.rho_cosphi(lat, height))

    # rho*cos(lat): 0.411775501279

Compute the radius of the parallel circle at a given latitude::

    print_me("Radius of parallel circle at latitude 65d 45' 30'' (meters)", e.rp(lat))

    # Radius of parallel circle at latitude 65d 45' 30'' (meters): 2626094.91467

Compute the radius of curvature of the Earth's meridian at given latitude::

    print_me("Radius of Earth's meridian at latitude 65d 45' 30'' (meters)", e.rm(lat))

    # Radius of Earth's meridian at latitude 65d 45' 30'' (meters): 6388705.74543

It is easy to compute the linear velocity at different latitudes::

    print_me("Linear velocity at the Equator (meters/second)", e.linear_velocity(0.0))

    # Linear velocity at the Equator (meters/second): 465.101303151

    print_me("Linear velocity at latitude 65d 45' 30'' (meters/second)", e.linear_velocity(lat))

    # Linear velocity at latitude 65d 45' 30'' (meters/second): 191.497860977

Finally, let's compute the distance between two points on the Earth:

- Bangkok: 13d 14' 09'' North, 100d 29' 39'' East
- Buenos Aires: 34d 36' 12'' South,  58d 22' 54'' West

.. note:: We will consider that positions 'East' and 'South' are negative

Here we will take advantage of facilities provided by ``Angle`` class::

    lon_ban = Angle(-100, 29, 39.0)

    lat_ban = Angle(13, 14, 9.0)

    lon_bai = Angle(58, 22, 54.0)

    lat_bai = Angle(-34, 36, 12.0)

    dist, error = e.distance(lon_ban, lat_ban, lon_bai, lat_bai)

    print_me("The distance between Bangkok and Buenos Aires is (km)", round(dist/1000.0, 2))

    # The distance between Bangkok and Buenos Aires is (km): 16832.89

    print_me("The approximate error of the estimation is (meters)", round(error, 0))

    # The approximate error of the estimation is (meters): 189.0

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
