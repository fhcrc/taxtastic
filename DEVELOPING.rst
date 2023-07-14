======================
 Developing taxtastic
======================

Requirements
============

Note that building docs, publishing to pypi, etc require some
additional dependencies. In your active virtualenv::

  pip3 install -r requirements-dev.txt

Git workflow
============

We aspire to more or less use a `feature branch workflow
<https://www.atlassian.com/git/workflows#!workflow-feature-branch>`_
for development. Briefly (for those working on the main fork):

* Features or bugfixes should start life as a GitHub issue
* Work on the feature occurs in a "feature branch" named like
  '%i-brief-description' % issue_number
* When completed (with tests passing) the feature branch is merged
  into dev (a pull request at this point might be appropriate if you
  want to request a code review).
* When it's time for a release, dev is merged into master (as a
  result, the head of the `master` branch is always on a release
  version).

Unit tests
==========

Run all unit tests from an active virtualenv as follows::

  python -m unittest discover

Run a subset of tests by providing a pattern::

  % python -m unittest discover -k test_subcommands.TestAddNode.test_new_nodes03 -v
  test_new_nodes03 (tests.test_subcommands.TestAddNode) ... ok

  ----------------------------------------------------------------------
  Ran 1 test in 0.151s

  OK

Preparing a release
===================

1. check out dev and make sure it is up to date with GitHub
1. update CHANGES.rst (the topmost section header should indicate
   the new version number)
1. make a final commit to dev
1. `git push origin dev`
1. checkout master and merge dev
1. tag the commit to reflect the new version number:
   `git tag -a -m "some message" vX.Y.Z`
1. `git push origin master`
1. `git push --tags`
1. update PyPi (see below)

PyPi
====

If you have not done so create a ~/.pypirc file::

  python setup.py register

Proceed to build and upload (assuming python3 in an active virtualenv)::

  python setup.py clean
  pip install -e .
  rm -r dist
  python3 setup.py sdist bdist_wheel
  python3 -m twine upload dist/*

Building docs with Sphinx
=========================

It's best to create and activate a virtualenv first to provide all of
the requirements for building and publishing the documentation (see
above).

The Sphinx configuration uses the version defined in the package,
which in turn uses the git tag, so make sure that the tag is up to
date before building the docs::

  (cd docs && make html)

Note that the html directory needs to contain a file named `.nojekyll`
to prevent GitHub from ignoring pages with leading underscores (like
`_static`), so the Makefile adds one

Manual deployment should not be necessary since there is a GH action
to build and deploy. However, is manual deployment is necessary, use
`ghp-import`::

  ghp-import -p html
