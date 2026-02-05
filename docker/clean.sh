#!/bin/bash
set -euo pipefail

docker rm -f rivrs-linalg
docker volume rm rivrs-linalg-workspace