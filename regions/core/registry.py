# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides a RegionsRegistry class.
"""

__all__ = []


class RegionsRegistry:
    """
    Class to hold registry to read, write, parse, and serialize regions
    in various formats.
    """

    registry = {}

    @classmethod
    def register(cls, classname, methodname, filetype):
        def inner_wrapper(wrapped_func):
            key = (classname, methodname, filetype)
            if key in cls.registry:
                raise ValueError(f'{methodname} for {filetype} is already '
                                 f'registered for {classname}')
            cls.registry[key] = wrapped_func
            return wrapped_func
        return inner_wrapper
