#!/bin/bash

# run from the top level directory to initialize gh-pages submodule
# from a newly-cloned git repo, eg:
#   cd bioscons
#   dev/ghpages_init.sh

git submodule init && \
git submodule update && \
cd html && \
git checkout gh-pages

