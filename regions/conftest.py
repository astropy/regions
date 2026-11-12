# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Configuration file for the pytest test suite.
"""

import pytest
from astropy.wcs.wcsapi import BaseLowLevelWCS


def _shares_wcs(item):
    """
    Return whether a test uses a WCS object that is shared between the
    threads of pytest-run-parallel.

    A WCS is shared when it comes from a fixture, whose value all
    threads receive, or from a class attribute.
    """
    if any(name.endswith('wcs') for name in getattr(item, 'fixturenames', ())):
        return True
    cls = getattr(item, 'cls', None)
    return cls is not None and any(isinstance(value, BaseLowLevelWCS)
                                   for value in vars(cls).values())


def pytest_collection_modifyitems(items):
    """
    Run the tests that are not thread-safe in a single thread with
    pytest-run-parallel.

    Matplotlib is not thread-safe, and every test that uses it is
    skipped without it, so the skip reason identifies these tests.
    Every test that requests the ``request`` fixture uses it to call
    ``request.getfixturevalue``, which is not thread-safe. An astropy
    WCS object can return wrong results when several threads use it at
    once on a free-threaded build, so tests that share one between
    threads are also run in a single thread.
    """
    for item in items:
        reasons = [marker.kwargs.get('reason', '')
                   for marker in item.iter_markers('skipif')]
        if any('matplotlib' in reason for reason in reasons):
            reason = 'matplotlib is not thread-safe'
        elif 'request' in getattr(item, 'fixturenames', ()):
            reason = 'request.getfixturevalue is not thread-safe'
        elif _shares_wcs(item):
            reason = 'a shared astropy WCS object is not thread-safe'
        else:
            continue
        item.add_marker(pytest.mark.thread_unsafe(reason=reason))
