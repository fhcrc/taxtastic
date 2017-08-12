#!/bin/bash

set -e

rm -rf dist
(cd .. && python setup.py sdist --dist-dir docker/dist)
tarball=$(ls dist/taxtastic*)
name=$(basename ${tarball%.tar.gz})
tag=taxtastic:${name#taxtastic-}

cp $tarball taxtastic.tar.gz
docker build -t $tag .
