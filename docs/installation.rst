Installation
============

``taxtastic`` requires Python 2.7.  The simplest method of installing
is using `pip <http://pip-installer.org>`_::

    pip install taxtastic

If you don't have pip, download the latest version from
https://pypi.python.org/pypi/taxtastic, decompress the archive, enter
the directory containing the source code, and install with::

    python setup.py install

Note that either method above will attempt to download and install
several additional Python packages, which, in turn have their own
system requirements. For example, at least one of the dependencies
require that the Python libraries are available to compile C
extensions. On Ubuntu (or another system using apt), the Python
libraries can be installed with::

  sudo apt-get install python-dev

Although you can certainly install to the system, we recommend
installing ``taxtastic`` in a `virtualenv
<http://virtualenv.readthedocs.org/en/latest/virtualenv.html>`_ to
avoid conflicts with system packages that may already be installed (or
that may be installed later).

So to put this all together, the following will result in
``taxtastic`` and dependencies installed in a new virtualenv named
``taxit-env`` (tested on Ubuntu 12.04 and 14.04)::

  sudo apt-get update
  sudo apt-get install python-pip python-virtualenv python-dev
  virtualenv taxit-env
  source taxit-env/bin/activate
  pip install taxtastic
  taxit -h

