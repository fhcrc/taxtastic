======================
 Developing taxtastic
======================

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

PyPI
====

Publishing to PyPI is handled automatically by the ``publish-pypi.yml``
GitHub Actions workflow when a release is created on GitHub.

Building docs with Sphinx
=========================

Install the package and Sphinx into an active virtualenv::

  pip install -e . sphinx

If subcommand help text has changed, regenerate it first::

  python dev/helptext.py -d docs/_helptext ./taxit.py

Then build the docs::

  cd docs && make html

Open ``html/index.html`` in a browser to review.

For iterative development, use ``sphinx-autobuild`` which watches for
changes and auto-reloads at ``http://127.0.0.1:8000``::

  pip install sphinx-autobuild
  sphinx-autobuild docs html

Deployment is handled automatically by the ``gh-pages.yml`` GitHub
Actions workflow on push to master.
