from importlib.metadata import version, PackageNotFoundError

# https://setuptools-scm.readthedocs.io/en/latest/usage/#at-runtime
try:
    __version__ = version("taxtastic")
except PackageNotFoundError:
    # package is not installed
    pass
