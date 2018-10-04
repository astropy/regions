.. include:: references.txt

.. warning::
    This ``regions`` package is in a very early stage of development.
    It is neither feature complete nor API stable!
    That said, please have a look and try to use it for your applications.
    Feedback and contributions are welcome!


#############################
Astropy Regions Documentation
#############################

This is an in-development `affiliated package`_ of `Astropy`_ for region handling.

To get an overview of available features, see :ref:`gs`.

The goal is to merge the functionality from `pyregion`_ and `photutils`_ apertures
and then after some time propose this package for inclusion in the Astropy core.

* Code : `Github repository`_
* Docs : `Region documentation`_
* Contributors : https://github.com/astropy/regions/graphs/contributors
* Releases: https://pypi.python.org/pypi/regions

User Documentation
==================

.. toctree::
   :maxdepth: 1

   installation
   getting_started
   shapes
   moc
   contains
   compound
   masks
   plotting
   ds9
   crtf
   fits_region
   shapely
   changelog

These plots are made using :download:`plot_image.reg`.

+----------------------------------------+----------------------------------------+
| ds9                                    | regions + matplotlib                   |
+========================================+========================================+
| .. image::    _static/region_ds9.png   | .. plot:: plot_reg.py                  |
|   :width: 300px                        |   :width: 300px                        |
|   :target:    _static/region_ds9.png   |                                        |
+----------------------------------------+----------------------------------------+

Advanced
========

.. toctree::
   :maxdepth: 1

   api
   development

Reporting Issues
================

If you have found a bug in Regions please report it by creating a
new issue on the `Regions GitHub issue tracker
<https://github.com/astropy/regions/issues>`_.

Please include an example that demonstrates the issue that will allow
the developers to reproduce and fix the problem. You may also be asked to
provide information about your operating system and a full Python
stack trace. The developers will walk you through obtaining a stack
trace if it is necessary.

Astropy Regions uses a package of utilities called `astropy-helpers
<https://github.com/astropy/astropy-helpers>`_ during building and
installation. If you have any build or installation issue mentioning
the ``astropy_helpers`` or ``ah_bootstrap`` modules please send a
report to the `astropy-helpers issue tracker
<https://github.com/astropy/astropy-helpers/issues>`_. If you are
unsure, then it's fine to report to the main Regions issue tracker.


Contributing
============

Like the `Astropy`_ project, Regions is made both by and for its
users. We accept contributions at all levels, spanning the gamut from
fixing a typo in the documentation to developing a major new feature.
We welcome contributors who will abide by the `Python Software
Foundation Code of Conduct
<https://www.python.org/psf/codeofconduct/>`_.
If you have a feature request or would like to contribute to ``regions``,
please go here: https://github.com/astropy/regions/

Regions follows the same workflow and coding guidelines as
`Astropy`_. The following pages will help you get started with
contributing fixes, code, or documentation (no git or GitHub
experience necessary):

* `How to make a code contribution <http://astropy.readthedocs.io/en/stable/development/workflow/development_workflow.html>`_

* `Coding Guidelines <http://docs.astropy.org/en/latest/development/codeguide.html>`_

* `Try the development version <http://astropy.readthedocs.io/en/stable/development/workflow/get_devel_version.html>`_

* `Developer Documentation <http://docs.astropy.org/en/latest/#developer-documentation>`_


Get Help
========

Besides github, you can `get help`_ from the community in a number of ways.
There is also a slack channel for regions hosted under the main astropy slack.