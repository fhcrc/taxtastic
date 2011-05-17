==============
 Installation
==============

dependencies
------------

Your best bet for installing dependencies will probably be to follow
instructions provided by each project. Bare-bones installation
instructions for installation from source is provided below.

Python 2.7
~~~~~~~~~~

At the moment, only this version of Python is supported (we use the
``argparse`` module, which is new in 2.7).


SQLAlchemy
~~~~~~~~~~

SQLAlchemy provides the ORM database layer::

 wget "http://prdownloads.sourceforge.net/sqlalchemy/SQLAlchemy-0.6.4.tar.gz?download" && \
 tar -xf SQLAlchemy-0.6.4.tar.gz && \
 cd SQLAlchemy-0.6.4 && \
 python setup.py install

xlrd
~~~~

xlrd is required to read Excel files (this is not an absolute
dependency, but some optons in ``taxit`` depend on it)::

 wget http://pypi.python.org/packages/source/x/xlrd/xlrd-0.7.1.tar.gz && \
 tar -xf xlrd-0.7.1.tar.gz && \
 cd xlrd-0.7.1 && \
 python setup.py install

docutils
~~~~~~~~

docutils is required only to compile docs::

 wget "http://prdownloads.sourceforge.net/docutils/docutils-0.7.tar.gz?download" && \
 tar -xf docutils-0.7.tar.gz && \
 cd docutils-0.7 && \
 python setup.py install

taxtastic package
----------------

This project is hosted by github. Create
a local repository and install as follows::

 git clone git@github.com:fhcrc/taxtastic.git
 cd taxtastic
 python setup.py install

You can also download the sources without using git by visiting
https://github.com/fhcrc/taxtastic and clicking on the "Downloads"
button on the right.

