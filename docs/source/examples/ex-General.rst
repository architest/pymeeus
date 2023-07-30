Properly Using PyMeeus
**********************

It is very common to try to run PyMeeus like this:

.. code:: sh

   import pymeeus

   mydate = pymeeus.Epoch(1992, 10, 13.0)

But if you do that, you'll get an error like this:

.. code:: sh

   Traceback (most recent call last):
     File "/home/user/test/test.py", line 3, in <module>
       epoch = pymeeus.Epoch(1992, 10, 13.0)
   AttributeError: module 'pymeeus' has no attribute 'Epoch'

This issue points to a misunderstanding that is very common in the
Python world. The keyword ``import`` is used to import **MODULES**\ ...
but PyMeeus is **NOT** a module: It is a **LIBRARY** composed of
**MULTIPLE** modules (``Angle``, ``Epoch``, ``Coordinates``, etc). As of
today, the library Pymeeus has 19 different modules (if you look into
the directory where ``pip`` stores the library, you'll find one ".py"
file per module).

Therefore if you want to use, for example, the module ``Angle`` you
should use:

.. code:: sh

   import pymeeus.Angle

I.e., your *module* is ``pymeeus.Angle``, and not just ``Angle``.

But there is more! When you use ``import`` to fetch a module, you must
then use the *dot* notation to access the components of the module
(classes, functions, etc). For instance:

.. code:: sh

   import pymeeus.Angle

   i = pymeeus.Angle.Angle(11.94524)

In this case, you are telling the Python interpreter that you want to
use the class ``Angle`` (with parameter '11.94524') from the module
``Angle`` belonging to the library ``pymeeus``.

There is, however, a more practical (and common) way to handle modules
using the statement ``from <MODULE> import <COMPONENT>``. For instance:

.. code:: sh

   from pymeeus.Angle import Angle
   from pymeeus.Epoch import Epoch, JDE2000
   from math import sin, cos, tan, acos, atan2, sqrt, radians, log10

This way is preferred because, among other reasons, only the required
components are loaded into memory instead of the whole module. Also, now
the component is directly added to your execution environment, which
means that you no longer need to use the *dot* notation.

Therefore, the script at the beginning would become:

.. code:: sh

   from pymeeus.Epoch import Epoch

   mydate = Epoch(1992, 10, 13.0)

