#!/bin/bash

# build from a local docker image

set -e

if [[ -z $1 ]]; then
    echo "usage: ./build.sh TAG [<outdir>]"
    echo
    echo "choose from one of the following tags:"
    echo
    docker image list taxtastic
    exit 1
fi

outdir=$(readlink -f ${2-.})

version=$1
docker_image=taxtastic:$1
docker_name=$(echo $docker_image | tr : _)_$RANDOM
img=$(readlink -f $outdir/taxtastic-${version}.img)

if [[ -f $img ]]; then
    echo "$img already exists"
    exit 1
fi

singularity create --size 1500 $img
docker run --name $docker_name $docker_image /bin/true
docker export $docker_name | sudo singularity import $img
docker rm $docker_name
sudo singularity exec --writable $img mkdir -p /fh /app /mnt
