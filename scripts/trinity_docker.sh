#!/usr/bin/env bash
set -euo pipefail

image="${TRINITY_DOCKER_IMAGE:-trinityrnaseq/trinityrnaseq:2.15.2}"
workspace="${TRINITY_DOCKER_WORKSPACE:-$(pwd)}"

exec docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${workspace}:${workspace}" \
  --workdir "${workspace}" \
  "${image}" \
  Trinity "$@"
