=============================
 unit tests for this package
=============================

Two shell scripts in this directory can be used to run some or all of
the unit tests defined in ./tests. The shell script `./test` can be
used to run all tests, with optional flags controlling verbosity::

  ./test -q # be as quiet as possible (loglevel = logging.ERROR) 
  ./test    # not much output         (loglevel = logging.WARNING) 
  ./test -v # a bit more verbose      (loglevel = logging.INFO)
  ./test -v # you get the idea        (loglevel = logging.DEBUG)

Individual unit tests can be run using `./testone` using the package
hierarchy in the tests directory::

  % ./testone tests.test_taxtastic.TestHelp.test01 
  --> ./taxtastic.py -h
  .
  ----------------------------------------------------------------------
  Ran 1 test in 0.219s

  OK

  % ./testone -v tests.test_taxtastic.TestHelp.test01
  test01 (tests.test_taxtastic.TestHelp) ... WARNING test_taxtastic 35 --> ./taxtastic.py -h
  INFO test_taxtastic 37 usage: taxtastic.py [-h] [-V] {create,help,check} ...

  To Be Named -- Creation, validation, and modification of reference packages
  for use with `pplacer` and related software.

  positional arguments:
    {create,help,check}
      help               Detailed help for actions using `help <action>`
      create             Create a reference package
      check              The check action is not yet implemented

  optional arguments:
    -h, --help           show this help message and exit
    -V, --version        Print the version number and exit
  ok

  ----------------------------------------------------------------------
  Ran 1 test in 0.228s

  OK

Output of unit tests are written to ./test_output.

unit tests for scripts
======================

Note that (as in the above example) scripts are tested as well. A base
class inheriting from unittest.TestBase is defined in
tests/config.py - see the docstring and existing unit tests for
instructions and examples.

