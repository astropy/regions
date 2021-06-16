# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides a RegionsRegistry class.
"""

__all__ = []


class IORegistryError(Exception):
    """Custom error for registry errors."""
    pass


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

    @classmethod
    def get_identifiers(cls, classname):
        return [key for key in cls.registry.keys()
                if key[0] == classname and key[1] == 'identify']

    @classmethod
    def get_format(cls, filename, classname, methodname):
        format = None
        identifiers = cls.get_identifiers(classname)
        if identifiers:
            for identifier in identifiers:
                if cls.registry[identifier](methodname, filename):
                    format = identifier[2]
                    break  # finds the first valid filetype

        if format is None:
            raise IORegistryError('Format could not be identified based on '
                                  'the file name or contents, please provide '
                                  'a "format" argument.')
        return format

    @classmethod
    def read(cls, filename, classname, format=None, **kwargs):
        """
        Read in a regions file.
        """
        if format is None:
            format = cls.get_format(filename, classname, 'read')
        key = (classname, 'read', format)
        return cls.registry[key](filename, **kwargs)

    @classmethod
    def parse(cls, data, classname, format=None, **kwargs):
        """
        Parse a regions string or table.
        """
        key = (classname, 'parse', format)
        return cls.registry[key](data, **kwargs)

    @classmethod
    def write(cls, regions, filename, classname, format=None, **kwargs):
        """
        Write to a regions file.
        """
        if format is None:
            format = cls.get_format(filename, classname, 'write')
        key = (classname, 'write', format)
        return cls.registry[key](regions, filename, **kwargs)

    @classmethod
    def serialize(cls, regions, classname, format=None, **kwargs):
        """
        Serialize to a regions string or table.
        """
        key = (classname, 'serialize', format)
        return cls.registry[key](regions, **kwargs)
