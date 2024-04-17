************
Installation
************

Requirements
============

Regions has the following strict requirements:

* `Python <https://www.python.org/>`_ 3.10 or later

* `NumPy <https://numpy.org/>`_ 1.23 or later

* `Astropy`_ 5.1 or later

Region also optionally depends on other packages for some features:

* `Matplotlib <https://matplotlib.org/>`_ 3.5 or later


Installing the latest released version
======================================

The latest released (stable) version of Regions can be installed either
with `pip`_ or `conda`_.

Using pip
---------

To install Regions with `pip`_, run::

    pip install regions

If you want to install Regions along with all of its optional
dependencies, you can instead do::

    pip install "regions[all]"

In most cases, this will install a pre-compiled version (called a
wheel) of Regions, but if you are using a very recent version of Python
or if you are installing Regions on a platform that is not common,
Regions will be installed from a source file. In this case you will
need a C compiler (e.g., ``gcc`` or ``clang``) to be installed for the
installation to succeed (see :ref:`building_source` prerequisites).

If you get a ``PermissionError``, this means that you do not have the
required administrative access to install new packages to your Python
installation.  In this case you may consider using the ``--user``
option to install the package into your home directory.  You can read
more about how to do this in the `pip documentation
<https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

Do **not** install Regions or other third-party packages using ``sudo``
unless you are fully aware of the risks.

Using conda
-----------

Regions can be installed with `conda`_ if you have installed
`Anaconda <https://www.anaconda.com/download>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.
To install Regions using the `conda-forge Anaconda channel
<https://anaconda.org/conda-forge/regions>`_, run::

    conda install -c conda-forge regions


.. _building_source:

Building from Source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python
and Numpy in order to build Regions from the source distribution. You
do not need to install any other specific build dependencies (such as
Cython) since these will be automatically installed into a temporary
build environment by `pip`_.

On Linux, using the package manager for your distribution will usually be
the easiest route.

On MacOS X you will need the `XCode`_ command-line tools, which can be
installed using::

    xcode-select --install

Follow the onscreen instructions to install the command-line tools
required.  Note that you do not need to install the full `XCode`_
distribution (assuming you are using MacOS X 10.9 or later).


Installing the development version
----------------------------------

Regions is being developed on `GitHub`_. The latest development version
of the Regions source code can be retrieved using git::

    git clone https://github.com/astropy/regions.git

Then to build and install Regions (with all of its optional
dependencies), run::

    cd regions
    pip install ".[all]"

If you wish to install the package in "editable" mode, instead include
the "-e" option::

    pip install -e ".[all]"

Alternatively, `pip`_ can be used to retrieve, build, and install the
latest development version from `GitHub`_::

    pip install "git+https://github.com/astropy/regions.git#egg=regions[all]"


Testing an installed Regions
============================

To test your installed version of Regions, you can run the test suite
using the `pytest`_ command. Running the test suite requires installing
the `pytest-astropy <https://github.com/astropy/pytest-astropy>`_ (0.11
or later) package.

To run the test suite, use the following command::

    pytest --pyargs photutils

Any test failures can be reported to the `Regions issue tracker
<https://github.com/astropy/regions/issues>`_.


.. _pip: https://pip.pypa.io/en/latest/
.. _conda: https://docs.conda.io/en/latest/
.. _GitHub: https://github.com/astropy/regions
.. _Xcode: https://developer.apple.com/xcode/
.. _pytest: https://docs.pytest.org/en/latest/
