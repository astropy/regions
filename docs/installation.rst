.. include:: references.txt

.. _install:

************
Installation
************

The regions package requires the following packages:

* Python 2.7 or 3.4 and above
* `Numpy <http://www.numpy.org>`_ 1.9 or later
* `Astropy <http://www.astropy.org>`__ 1.2 or later

In addition, the following packages are needed for optional functionality:

* `Matplotlib <http://www.matplotlib.org>`__ 1.5 or later
* `Shapely <http://toblerity.org/shapely/manual.html>`__

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

Using pip
---------

To install regions with `pip <http://www.pip-installer.org/en/latest/>`_
from `PyPI <https://pypi.python.org/pypi/regions>`_, run::

    pip install regions --no-deps

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

To check if there are any issues with your installation, you can run the tests:

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
