#!/bin/bash

if [[ -z $1 ]]; then
    echo "usage: ./build.sh TAG [<outdir>]"
    exit 1
fi

outdir=$(readlink -f ${2-.})

version=$1
singularity_version=$(singularity --version)
img=$(readlink -f $outdir/taxtastic-${version}-singularity${singularity_version}.img)

if [[ -f $img ]]; then
    echo "$img already exists"
    exit 1
fi

echo $img

rm -rf dist
(cd .. && python setup.py sdist --dist-dir singularity/dist)
tarball=$(ls dist/taxtastic*)
name=$(basename ${tarball%.tar.gz})

cp $tarball taxtastic.tar.gz
cp ../dev/install_pplacer.sh .
sudo singularity build $img Singularity
