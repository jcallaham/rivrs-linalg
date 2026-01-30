#!/bin/bash
set -euo pipefail

docker rm -f csrrs-dev
docker volume rm csrrs-workspace