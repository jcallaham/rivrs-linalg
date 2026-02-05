#!/bin/bash
set -euo pipefail

docker rm -f rivrs-linalg-dev
docker volume rm rivrs-linalg-workspace