======================
 Developing taxtastic
======================

Requirements
============

Note that building docs, publishing to pypi, etc require some
additional dependencies. It's probably best to work in a virtualenv::

  virtualenv taxtastic-env
  pip install -r requirements.txt


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
1. update the docs (details below): `(cd docs && make html)`
1. publish the updated docs: `ghp-import -p html`
1. `git push origin master`
1. `git push --tags`
1. update PyPi (see below)

PyPi
====

If you have not done so create a ~/.pypirc file::

  python setup.py register

Proceed to build and upload::

  python setup.py clean
  rm -r dist
  python setup.py sdist bdist_wheel
  twine upload dist/*

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

Thankfully, publishing to the GitHub page for the bioscons repository
is easy using `ghp-import`::

  ghp-import -p html

