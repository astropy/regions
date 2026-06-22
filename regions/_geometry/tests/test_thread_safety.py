# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test that the Cython overlap functions return identical results when run
concurrently from multiple threads.

This catches data races and shared-state issues introduced by ``with
nogil:`` sections in the functions.
"""
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from regions._geometry import (circle_overlap_grid, ellipse_overlap_grid,
                               polygon_overlap_grid, rectangle_overlap_grid)
from regions._geometry.polygon_contains import (points_in_polygon,
                                                points_in_polygon_covers)

N_THREADS = 8
N_CALLS_PER_THREAD = 4


def _run_concurrent(fn, *, n_threads=N_THREADS,
                    n_calls=N_CALLS_PER_THREAD):
    expected = fn()
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futures = [ex.submit(fn)
                   for _ in range(n_threads * n_calls)]
        for fut in futures:
            np.testing.assert_array_equal(fut.result(), expected)


@pytest.mark.parametrize('use_exact', [1, 0])
def test_circle_overlap_grid_threadsafe(use_exact):
    def fn():
        return circle_overlap_grid(-15.0, 15.0, -15.0, 15.0, 200, 200, 8.0,
                                   use_exact, 5)
    _run_concurrent(fn)


@pytest.mark.parametrize('use_exact', [1, 0])
def test_ellipse_overlap_grid_threadsafe(use_exact):
    def fn():
        return ellipse_overlap_grid(-15.0, 15.0, -15.0, 15.0, 200, 200, 8.0,
                                    5.0, 0.7, use_exact, 5)
    _run_concurrent(fn)


@pytest.mark.parametrize('use_exact', [1, 0])
def test_rectangle_overlap_grid_threadsafe(use_exact):
    def fn():
        return rectangle_overlap_grid(-15.0, 15.0, -15.0, 15.0, 200, 200, 12.0,
                                      7.0, 0.5, use_exact, 5)
    _run_concurrent(fn)


@pytest.mark.parametrize('use_exact', [1, 0])
def test_polygon_overlap_grid_threadsafe(use_exact):
    n_vertices = 32
    angles = np.linspace(0.0, 2.0 * np.pi, n_vertices, endpoint=False)
    radii = 8.0 + 2.0 * np.sin(5.0 * angles)
    vx = np.ascontiguousarray(radii * np.cos(angles))
    vy = np.ascontiguousarray(radii * np.sin(angles))

    def fn():
        return polygon_overlap_grid(-15.0, 15.0, -15.0, 15.0, 200, 200,
                                    vx, vy, use_exact, 8)
    _run_concurrent(fn)


def test_mixed_functions_concurrent():
    """
    Run all geometry function types concurrently from many threads.

    This is a more aggressive smoke test that mixes different functions
    so that if any of them shares mutable state (e.g., a module-level
    static buffer), interference between calls would surface.
    """
    n_vertices = 16
    angles = np.linspace(0.0, 2.0 * np.pi, n_vertices, endpoint=False)
    vx = np.ascontiguousarray(7.0 * np.cos(angles))
    vy = np.ascontiguousarray(5.0 * np.sin(angles))

    expected_circ = circle_overlap_grid(-10.0, 10.0, -10.0, 10.0, 150, 150,
                                        6.0, 1, 5)
    expected_ell = ellipse_overlap_grid(-10.0, 10.0, -10.0, 10.0, 150, 150,
                                        7.0, 4.0, 0.3, 1, 5)
    expected_rect = rectangle_overlap_grid(-10.0, 10.0, -10.0, 10.0, 150, 150,
                                           9.0, 5.0, 0.4, 1, 5)
    expected_poly = polygon_overlap_grid(-10.0, 10.0, -10.0, 10.0, 150, 150,
                                         vx, vy, 1, 5)

    def task(which):
        if which == 0:
            return 'c', circle_overlap_grid(
                -10.0, 10.0, -10.0, 10.0, 150, 150, 6.0, 1, 5)
        if which == 1:
            return 'e', ellipse_overlap_grid(
                -10.0, 10.0, -10.0, 10.0, 150, 150, 7.0, 4.0, 0.3, 1, 5)
        if which == 2:
            return 'r', rectangle_overlap_grid(
                -10.0, 10.0, -10.0, 10.0, 150, 150, 9.0, 5.0, 0.4, 1, 5)
        return 'p', polygon_overlap_grid(
            -10.0, 10.0, -10.0, 10.0, 150, 150, vx, vy, 1, 5)

    expected_map = {'c': expected_circ, 'e': expected_ell,
                    'r': expected_rect, 'p': expected_poly}

    with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
        futures = [ex.submit(task, i % 4)
                   for i in range(N_THREADS * 8)]
        for future in futures:
            tag, result = future.result()
            assert_array_equal(result, expected_map[tag])


@pytest.mark.parametrize('func', [points_in_polygon, points_in_polygon_covers])
def test_points_in_polygon_threadsafe(func):
    n_vertices = 24
    angles = np.linspace(0.0, 2.0 * np.pi, n_vertices, endpoint=False)
    radii = 8.0 + 2.0 * np.sin(5.0 * angles)
    vx = np.ascontiguousarray(radii * np.cos(angles))
    vy = np.ascontiguousarray(radii * np.sin(angles))

    grid = np.linspace(-12.0, 12.0, 80)
    gx, gy = np.meshgrid(grid, grid)
    x = np.ascontiguousarray(gx.ravel())
    y = np.ascontiguousarray(gy.ravel())

    def fn():
        return func(x, y, vx, vy)
    _run_concurrent(fn)
