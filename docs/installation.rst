.. include:: references.txt

.. _install:

************
Installation
************

The regions package requires the following packages:

* Python 3.6 or later
* `Numpy <http://www.numpy.org>`_ 1.16 or later
* `Astropy <http://www.astropy.org>`__ 2.0 or later

In addition, the following packages are needed for optional functionality:

* `Matplotlib <https://matplotlib.org>`__ 2.0 or later

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

Using pip
---------

To install regions with `pip <https://pip.pypa.io/en/latest/>`_
from `PyPI <https://pypi.python.org/pypi/regions>`_, run::

    pip install regions --no-deps

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

Using conda
-----------

To install regions with `Anaconda <https://www.anaconda.com/download/>`_
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
<https://pypi.python.org/pypi/pytest-arraydiff>`_ package is installed
(version v0.3 or newer). Then, run the tests with:

.. code-block:: bash

    python setup.py test

To build the documentation, do:

.. code-block:: bash

    python setup.py build_docs
