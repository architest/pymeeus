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
