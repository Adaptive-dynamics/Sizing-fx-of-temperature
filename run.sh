#!/bin/bash

set -ex

#docker build -t docker.dragonfly.co.nz/adaptive_activity .

#docker push docker.dragonfly.co.nz/adaptive_activity

docker run --rm  -v $PWD:/work -w /work \
  docker.dragonfly.co.nz/dragonverse-17.04 ./build.sh
