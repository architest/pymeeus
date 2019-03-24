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

And now, let's compute the distance between two points on the Earth:

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

Let's now compute the geometric heliocentric position for a given epoch::

    epoch = Epoch(1992, 10, 13.0)

    lon, lat, r = Earth.geometric_heliocentric_position(epoch)

    print_me("Geometric Heliocentric Longitude", lon.to_positive())

    # Geometric Heliocentric Longitude: 19.9072721503

    print_me("Geometric Heliocentric Latitude", lat.dms_str(n_dec=3))

    # Geometric Heliocentric Latitude: -0.721''

    print_me("Radius vector", r)

    # Radius vector: 0.997608520236

And now, compute the apparent heliocentric position for the same epoch::

    epoch = Epoch(1992, 10, 13.0)

    lon, lat, r = Earth.apparent_heliocentric_position(epoch)

    print_me("Apparent Heliocentric Longitude", lon.to_positive())

    # Apparent Heliocentric Longitude: 19.9059856939

    print_me("Apparent Heliocentric Latitude", lat.dms_str(n_dec=3))

    # Apparent Heliocentric Latitude: -0.721''

    print_me("Radius vector", r)

    # Radius vector: 0.997608520236

Print mean orbital elements for Earth at 2065.6.24::

    epoch = Epoch(2065, 6, 24.0)

    l, a, e, i, ome, arg = Earth.orbital_elements_mean_equinox(epoch)

    print_me("Mean longitude of the planet", round(l, 6))

    # Mean longitude of the planet: 272.716028

    print_me("Semimajor axis of the orbit (UA)", round(a, 8))

    # Semimajor axis of the orbit (UA): 1.00000102

    print_me("Eccentricity of the orbit", round(e, 7))

    # Eccentricity of the orbit: 0.0166811

    print_me("Inclination on plane of the ecliptic", round(i, 6))

    # Inclination on plane of the ecliptic: 0.0

    print_me("Longitude of the ascending node", round(ome, 5))

    # Longitude of the ascending node: 174.71534

    print_me("Argument of the perihelion", round(arg, 6))

    # Argument of the perihelion: -70.651889

Find the epoch of the Perihelion closer to 2008/02/01::

    epoch = Epoch(2008, 2, 1.0)

    e = Earth.perihelion_aphelion(epoch)

    y, m, d, h, mi, s = e.get_full_date()

    peri = str(y) + '/' + str(m) + '/' + str(d) + ' ' + str(h) + ':' + str(mi)

    print_me("The Perihelion closest to 2008/2/1 happened on", peri)

    # The Perihelion closest to 2008/2/1 happened on: 2008/1/2 23:53

Compute the time of passage through an ascending node::

    epoch = Epoch(2019, 1, 1)

    time, r = Earth.passage_nodes(epoch)

    y, m, d = time.get_date()

    d = round(d, 1)

    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))

    # Time of passage through ascending node: 2019/3/15.0

    print("Radius vector at ascending node: {}".format(round(r, 4)))

    # Radius vector at ascending node: 0.9945

Compute the parallax correction::

    right_ascension = Angle(22, 38, 7.25, ra=True)

    declination = Angle(-15, 46, 15.9)

    latitude = Angle(33, 21, 22)

    distance = 0.37276

    hour_angle = Angle(288.7958)

    top_ra, top_dec = Earth.parallax_correction(right_ascension, declination,

                                                latitude, distance, hour_angle)

    print_me("Corrected topocentric right ascension: ", top_ra.ra_str(n_dec=2))

    # Corrected topocentric right ascension: : 22h 38' 8.54''

    print_me("Corrected topocentric declination", top_dec.dms_str(n_dec=1))

    # Corrected topocentric declination: -15d 46' 30.0''

Compute the parallax correction in ecliptical coordinates::

    longitude = Angle(181, 46, 22.5)

    latitude = Angle(2, 17, 26.2)

    semidiameter = Angle(0, 16, 15.5)

    obs_lat = Angle(50, 5, 7.8)

    obliquity = Angle(23, 28, 0.8)

    sidereal_time = Angle(209, 46, 7.9)

    distance = 0.0024650163

    topo_lon, topo_lat, topo_diam = \

        Earth.parallax_ecliptical(longitude, latitude, semidiameter, obs_lat,

                                  obliquity, sidereal_time, distance)

    print_me("Corrected topocentric longitude", topo_lon.dms_str(n_dec=1))

    # Corrected topocentric longitude: 181d 48' 5.0''

    print_me("Corrected topocentric latitude", topo_lat.dms_str(n_dec=1))

    # Corrected topocentric latitude: 1d 29' 7.1''

    print_me("Corrected topocentric semidiameter", topo_diam.dms_str(n_dec=1))

    # Corrected topocentric semidiameter: 16' 25.5''
