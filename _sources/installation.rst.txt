Installation
============

``taxtastic`` requires Python 2.7.  The simplest method of installing
is using `pip <http://pip-installer.org>`_::

  pip install taxtastic

We strongly recommend installation into a virtualenv. On a clean
Ubuntu 16.04 system, complete instructions for installing the
``taxtastic`` package and the ``taxit`` command line entry point in a
virtualenv are below. Note that python2.7 is not longer installed
by default in this OS::

  sudo apt-get update
  sudo apt-get install python2.7 python-virtualenv

Once python2 is installed, create a virtualenv and install ``taxtastic``::

  virtualenv taxtastic-env
  source taxtastic-env/bin/activate
  pip install -U pip
  pip install taxtastic

If you prefer to install from the git repository::

  git clone https://github.com/fhcrc/taxtastic.git
  cd taxtastic
  virtualenv taxtastic-env
  source taxtastic-env/bin/activate
  pip install .

If you want to live dangerously and install the package to the system
despite our pleas not to do so::

  sudo apt-get install python-pip
  sudo pip install taxtastic

If you are not familiar with python virtual environments, the
following post is helpful:
https://realpython.com/blog/python/python-virtual-environments-a-primer/
