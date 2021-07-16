.. include:: references.txt

###############
Astropy Regions
###############

**Regions** is an in-development `coordinated package`_ of `Astropy`_
for region handling.

To get an overview of available features, see :ref:`getting_started`.

The eventual goal is to merge this package into the Astropy core
package.

* Code : `GitHub repository`_
* Contributors : https://github.com/astropy/regions/graphs/contributors

.. warning::
    This ``regions`` package is in an early stage of development. It
    is neither feature complete nor API stable. That said, please
    have a look and try to use it for your applications. Feedback and
    contributions are welcome!


Getting Started
===============

.. toctree::
   :maxdepth: 1

   install
   getting_started
   contributing
   license
   changelog


User Documentation
==================

.. toctree::
   :maxdepth: 1

   shapes
   contains
   compound
   masks
   plotting
   region_io
   shapely
   api

+----------------------------------------+------------------------------------+
| ds9                                    | regions + matplotlib               |
+========================================+====================================+
| .. image::    _static/region_ds9.png   | .. plot:: plot_reg.py              |
|   :width: 300px                        |   :width: 300px                    |
|   :target:    _static/region_ds9.png   |                                    |
+----------------------------------------+------------------------------------+
