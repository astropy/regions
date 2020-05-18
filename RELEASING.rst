Release instructions
====================

Once the package is ready to release, use ``git tag`` to tag the
release::

    git tag -m <version> <version>

e.g::

    git tag -m v0.1 v0.1

You can also include the ``-s`` flag to sign the tag if you have
PGP keys set up. Then, push the tag to GitHub, e.g.::

    git push upstream v0.1

and the build should happen automatically on Azure pipelines. You can
follow the progress of the build here:

https://dev.azure.com/astropy-project/regions/_build

If there are any failures, you can always delete the tag, fix the
issues, tag the release again, and push the tag to GitHub.

See the `OpenAstronomy Azure Pipelines Templates Documentation
<https://openastronomy-azure-pipelines.readthedocs.io/en/latest/publish.html>`_
for more details about the Azure Pipelines set-up.
