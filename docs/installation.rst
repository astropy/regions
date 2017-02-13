.. include:: references.txt

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

Installing the latest stable version is possible either using pip or conda.

Using pip
---------

To install regions with `pip <http://www.pip-installer.org/en/latest/>`_
from `PyPI <https://pypi.python.org/pypi/regions>`_
simply run::

    pip install --no-deps regions

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

Using conda
-----------

To install regions with `Anaconda <https://www.continuum.io/downloads>`_
from the `astropy channel on anaconda.org <https://anaconda.org/astropy/regions>`__
simply run::

    conda install -c astropy regions


Testing installation
--------------------

To check if your install is OK, run the tests:

.. code-block:: bash

    python -c 'import regions; regions.test()'

Development version
===================

Install the latest development version from https://github.com/astropy/regions :

.. code-block:: bash

    git clone https://github.com/astropy/regions
    cd regions
    python setup.py install

To run the tests, you will need to make sure the `pytest-arraydiff
<https://pypi.python.org/pypi/pytest-arraydiff>`_ package is installed. Then,
run the tests with:

.. code-block:: bash

    python setup.py test

To build the documentation, do:

.. code-block:: bash

    python setup.py build_docs

Optional dependencies
=====================

The following packages are optional dependencies, install if needed:

* `shapely`_ for advanced pixel region operations
* `matplotlib`_ for plotting regions
* maybe `spherical_geometry`_ for polygons (not used yet)
