import os

ver = os.path.join(os.path.dirname(__file__), 'data', 'ver')
if 'TAXTASTIC_VERSION' in os.environ:
    __version__ = os.environ['TAXTASTIC_VERSION']
elif os.path.isfile(ver):
    with open(ver) as f:
        __version__ = f.read().strip().replace('-', '+', 1).replace('-', '.')
        __version__ = __version__.lstrip('v')
else:
    __version__ = ''
