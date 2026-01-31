#!/bin/bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}CSRRS Docker Run Script${NC}"
echo -e "${GREEN}=====================================${NC}"
echo ""

# Check if GitHub CLI is authenticated
if ! gh auth status &>/dev/null; then
	echo -e "${YELLOW}Warning: GitHub CLI not authenticated${NC}"
	echo "Git/GitHub operations may not work. Please run:"
	echo "  gh auth login"
	echo ""
	read -p "Continue anyway? (y/N) " -n 1 -r
	echo
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		exit 1
	fi
fi

# Check if container already exists
if docker ps -a --format '{{.Names}}' | grep -q '^csrrs-dev$'; then
	echo -e "${YELLOW}Container 'csrrs-dev' already exists${NC}"

	if docker ps --format '{{.Names}}' | grep -q '^csrrs-dev$'; then
		echo "Container is running. Attaching..."
		docker exec -it csrrs-dev bash
	else
		echo "Container is stopped. Starting..."
		docker start csrrs-dev
		docker exec -it csrrs-dev bash
	fi
else
	echo -e "${GREEN}Creating new container 'csrrs-dev'...${NC}"
	echo ""

	# Create named volumes if they don't exist
	docker volume create csrrs-cargo-cache 2>/dev/null || true
	docker volume create csrrs-sccache-cache 2>/dev/null || true

	# Get GitHub token from host (falls back to file-based if keyring fails)
	GH_TOKEN=$(gh auth token 2>/dev/null || echo "")

	# Run the container
	docker run -it \
		--name csrrs-dev \
		--platform linux/arm64 \
		-v csrrs-workspace:/workspace/csrrs \
		-v csrrs-cargo-cache:/root/.cargo/registry \
		-v csrrs-sccache-cache:/root/.cache/sccache \
		-e GH_TOKEN="${GH_TOKEN}" \
		-e GITHUB_TOKEN="${GH_TOKEN}" \
		-e RUSTFLAGS="-C target-cpu=native" \
		-e CARGO_BUILD_JOBS=8 \
		-e RUSTC_WRAPPER=sccache \
		--cpus=8 \
		--memory=16g \
		csrrs-dev:latest \
		/bin/bash
fi

echo ""
echo -e "${GREEN}Exited container${NC}"
echo ""
echo "To restart: docker start csrrs-dev && docker exec -it csrrs-dev bash"
echo "To remove:  docker rm csrrs-dev"
