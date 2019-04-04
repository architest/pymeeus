PyMeeus
=======

    **Library of astronomical algorithms in Python**.

PyMeeus is a Python implementation of the astronomical algorithms
described in the classical book "Astronomical Algorithms, 2nd Edition,
Willmann-Bell Inc. (1998)" by Jean Meeus.

There are great astronomical libraries out there. For instance, if
you're looking for high precision and speed you should take a look at
`libnova <http://libnova.sourceforge.net/>`__. For a set of python
modules aimed at professional astronomers, you should look at
`Astropy <http://www.astropy.org/>`__. On the other hand, the advantages
of PyMeeus are its simplicity, ease of use, ease of reading, ease of
installation (it has the minimum amount of dependencies) and abundant
documentation.

Installation
------------

The easiest way of installing PyMeeus is using pip:

.. code:: sh

    pip install pymeeus

Or, for a per-user installation:

.. code:: sh

    pip install --user pymeeus

If you prefer Python3, you can use:

.. code:: sh

    pip3 install --user pymeeus

Meta
----

Author: Dagoberto Salazar

Distributed under the GNU Lesser General Public License v3 (LGPLv3). See
``LICENSE.txt`` and ``COPYING.LESSER`` for more information.

Documentation: https://pymeeus.readthedocs.io/en/latest/

GitHub: https://github.com/architest/pymeeus

If you have Sphinx installed, you can generate your own, latest
documentation going to directory 'docs' and issuing:

.. code:: sh

    make html

Then the HTML documentation pages can be found in 'build/html'.

Contributing
------------

The preferred method to contribute is through forking and pull requests:

1. Fork it (https://github.com/architest/pymeeus/fork)
2. Create your feature branch (``git checkout -b feature/fooBar``)
3. Commit your changes (``git commit -am 'Add some fooBar'``)
4. Push to the branch (``git push origin feature/fooBar``)
5. Create a new Pull Request

Please bear in mind that PyMeeus follows the PEP8 style guide for Python
code `(PEP8) <https://www.python.org/dev/peps/pep-0008/?>`__. We suggest
you install and use a linter like
`Flake8 <http://flake8.pycqa.org/en/latest/>`__ before contributing.

Additionally, PyMeeus makes heavy use of automatic tests. As a general
rule, every function or method added must have a corresponding test in
the proper place in ``tests`` directory.

Finally, documentation is also a big thing here. Add proper and abundant
documentation to your new code. This also includes in-line comments!!!.

What's new
----------

-  0.3.6

   -  Add method to compute rising and setting times of the Sun

-  0.3.5

   -  Add method magnitude() to planet classes

-  0.3.4

   -  Add method to compute the parallax correction to Earth class

-  0.3.3

   -  Add methods to compute the passage through the nodes

-  0.3.2

   -  Add methods to compute the perihelion and aphelion of all planets

-  0.3.1

   -  Fix errors in the elongation computation, add tests and examples
      of use of methods 'geocentric\_position()', and tests and examples
      for 'Pluto' class

-  0.3.0

   -  Added 'Pluto' class

-  0.2.11

   -  Added conjunction and opposition methods for Jupiter, Saturn,
      Uranus and Neptune

-  0.2.10

   -  Added 'geocentric\_position()' method to 'Minor' class, and added
      conjunction and opposition methods for Mercury, Venus and Mars.

-  0.2.9

   -  Added class 'Minor', as well as functions to compute velocity of
      an object and length of an orbit

-  0.2.8

   -  Added methods 'geocentric\_position()' to all the planets

-  0.2.7

   -  Added function 'kepler\_equation()' to 'Coordinates' module, and
      'orbital\_elements' methods to classes 'Mercury', 'Venus',
      'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus' and 'Neptune'

-  0.2.6

   -  Added classes 'Uranus' and 'Neptune', plus additional functions in
      Coordinates module'

-  0.2.5

   -  Added classes 'Jupiter' and 'Saturn'

-  0.2.4

   -  Minor bug fixing, added methods
      'ephemeris\_physical\_observations()' and
      'beginning\_synodic\_rotation()', and added classes 'Mercury' and
      'Mars'

-  0.2.3

   -  Added the complete list of VSOP87 parameters to Venus, method to
      compute dates of equinoxes and solstices, and the Equation of Time

-  0.2.2

   -  Added heliocentric position method for J2000.0 (Earth) and
      rectangular coordinates methods (Sun)

-  0.2.1

   -  Added Venus module and VSOP87-based positioning methods

-  0.2.0

   -  Added Sun module

-  0.1.10

   -  Added functions to compute if three objects are in a straight
      line, and the smallest diameter of a circle containing them.

-  0.1.9

   -  Added functions to compute relative position angles and
      conjunctions

-  0.1.8

   -  Added functions to compute angular separation to Coordinates
      module

-  0.1.7

   -  Added functions to compute atmospheric refraction to Coordinates
      module

-  0.1.6

   -  Added function 'times\_rise\_transit\_set()' to Coordinates module

-  0.1.5

   -  Added functions for parallactic angle, ecliptic points in the
      horizon, angle between north celestial pole and north pole of the
      ecliptic, and diurnal path vs. the horizon at time of rising or
      setting

-  0.1.4

   -  Added several conversion functions to Coordinates module

-  0.1.3

   -  Added Coordinates module

-  0.1.2

   -  Added precession and proper motion methods, and changed handling
      of Epoch class

-  0.1.1

   -  Added methods related to nutation corrections

-  0.1.0

   -  Earth class added

-  0.0.9

   -  Significant documentation improvements

-  0.0.8

   -  Epoch class finished

-  0.0.7

   -  Epoch class added

-  0.0.6

   -  CurveFitting class added

-  0.0.5

   -  Interpolation class added

-  0.0.4

   -  Angle class finished

-  0.0.3

   -  Removed unnecessary dependencies

-  0.0.2

   -  Documentation improvements

-  0.0.1

   -  Initial commit
