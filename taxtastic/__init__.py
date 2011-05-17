def _safeint(s):
    try:
        return int(s)
    except ValueError:
        return s

__version__ = "0.2.0"
__version_info__ = tuple([_safeint(num) for num in __version__.split('.')])


