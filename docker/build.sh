#!/bin/bash

# docker tags are formatted like "taxtastic:0.9.2"

if [[ -z $1 ]]; then
    git_tag=$(git describe --tags)
    tag=taxtastic:${git_tag#v}
else
    tag=taxtastic:$1
fi

echo "building image with tag $tag"

set -e

rm -rf dist
(cd .. && python setup.py sdist --dist-dir docker/dist)
tarball=$(ls dist/taxtastic*)

cp $tarball taxtastic.tar.gz
cp ../dev/install_pplacer.sh .
docker build -t $tag .
