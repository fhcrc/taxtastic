#!/bin/bash

# see https://docs.docker.com/docker-cloud/builds/push-images/

set -e

if [[ -z "$1" ]]; then
    echo "usage ./push.sh <tag>"
    echo "example: "
    echo "$ ./push 0.6.4-12-g04409ee"
    echo "choose from tags below:"
    docker images taxtastic
    exit 1
fi

image="taxtastic:$1"

export DOCKER_ID_USER="nghoffman"
# docker login

# ensure that image exists
docker image history "$image" > /dev/null

docker tag "$image" "$DOCKER_ID_USER/$image"
docker push "$DOCKER_ID_USER/$image"

docker tag "$image" "$DOCKER_ID_USER/taxtastic:latest"
docker push "$DOCKER_ID_USER/taxtastic:latest"
