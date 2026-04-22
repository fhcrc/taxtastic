Installation
============

``taxtastic`` requires Python 3.8 or later. The simplest method of
installation is using `uv <https://docs.astral.sh/uv/>`_::

  uv tool install taxtastic

Or to install from the git repository::

  uv tool install git+https://github.com/fhcrc/taxtastic.git

For development we recommend installing into a virtual environment::

  python -m venv .env
  source .env/bin/activate
  pip install taxtastic

Or installing using the `uv tool install --editable` flag::

  uv tool install --editable taxtastic
