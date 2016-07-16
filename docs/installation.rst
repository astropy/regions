.. _install:

************
Installation
************

* Python 2.7 and 3.4+ are supported.
* The only required dependency for ``regions`` is Astropy (version 1.2 or later).

The ``regions`` package works like most other Astropy affiliated packages.
Since it is planned to be merged into the Astropy core, we didn't put much
effort into writing up installation instructions for this separate package.

Stable version
==============

Install latest stable version from https://pypi.python.org/pypi/regions :

.. code-block:: bash

    pip install regions

(conda package coming soon, not available yet.)

To check if your install is OK, run the tests:

.. code-block:: bash

    python -c 'import regions; regions.test()'

Development version
===================

Install the latest development version:

.. code-block:: bash

    git clone https://github.com/astropy/regions
    cd regions
    python setup.py install
    python setup.py test
    python setup.py build_sphinx

Optional dependencies
=====================

The following packages are optional dependencies, install if needed:

* Shapely for advanced pixel region operations
* matplotlib and wcsaxes for plotting regions
* maybe https://github.com/spacetelescope/sphere
