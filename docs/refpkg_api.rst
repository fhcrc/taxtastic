Refpkg API
----------

Taxtastic provides a Python API for creating and manipulating refpkgs wihtout having to recourse to the command line.  It is entirely based around the ``Refpkg`` class in ``taxtastic.refpkg``.  Thus all scripts dealing with refpkgs will include code that looks something like::

    from taxtastic.refpkg import Refpkg

    r = Refpkg('/path/to/refpkg')

The rest of the API is implemented as methods of the ``Refpkg`` class.
