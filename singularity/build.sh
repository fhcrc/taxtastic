#!/bin/bash

set -e

if [[ $1 == '-h' ]]; then
    echo "usage: ./build.sh [<version>] [<outdir>]"
    echo "builds :latest by default"
    exit
fi

outdir=$(readlink -f ${2-.})

if [[ -z $1 ]]; then
    version=latest
    tag=latest
else
    version=$1
    tag=release-$1
fi

img=taxtastic-${version}.img
singfile=$(mktemp Singularity-XXXXXX)
sed s"/TAG/$tag/" < Singularity > $singfile

if [[ ! -f $img ]]; then
    singularity create --size 1500 $img
    sudo singularity bootstrap $img $singfile
fi

rm $singfile

