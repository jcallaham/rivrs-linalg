#!/bin/bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Detect host architecture and set target platform
ARCH=$(uname -m)
case "$ARCH" in
	arm64|aarch64) PLATFORM="linux/arm64" ;;
	x86_64)        PLATFORM="linux/amd64" ;;
	*)             PLATFORM="linux/$ARCH"  ;;
esac

echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}rivrs-linalg Docker Run Script${NC}"
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
if docker ps -a --format '{{.Names}}' | grep -q '^rivrs-linalg$'; then
	echo -e "${YELLOW}Container 'rivrs-linalg' already exists${NC}"

	if docker ps --format '{{.Names}}' | grep -q '^rivrs-linalg$'; then
		echo "Container is running. Attaching..."
		docker attach rivrs-linalg
	else
		echo "Container is stopped. Starting and attaching..."
		docker start -ai rivrs-linalg
	fi
else
	echo -e "${GREEN}Creating new container 'rivrs-linalg'...${NC}"
	echo ""

	# Create named volumes if they don't exist
	docker volume create rivrs-linalg-workspace 2>/dev/null || true
	docker volume create rivrs-linalg-cargo-cache 2>/dev/null || true
	docker volume create rivrs-linalg-sccache-cache 2>/dev/null || true
	docker volume create rivrs-linalg-claude-config 2>/dev/null || true

	# Get GitHub token from host (falls back to file-based if keyring fails)
	GH_TOKEN=$(gh auth token 2>/dev/null || echo "")

	# Cargo registry token for publishing to crates.io
	if [ -z "${CARGO_REGISTRY_TOKEN:-}" ]; then
		echo -e "${YELLOW}Warning: CARGO_REGISTRY_TOKEN not set — 'cargo publish' will not work${NC}"
		echo "Set it in your shell environment before running this script."
		echo ""
	fi

	# Run the container
	docker run -it \
		--name rivrs-linalg \
		--platform "$PLATFORM" \
		--cap-add PERFMON \
		--security-opt seccomp=unconfined \
		-v rivrs-linalg-workspace:/workspace \
		-v rivrs-linalg-cargo-cache:/home/node/.cargo/registry \
		-v rivrs-linalg-sccache-cache:/home/node/.cache/sccache \
		-v rivrs-linalg-claude-config:/home/node/.claude \
		-e GH_TOKEN="${GH_TOKEN}" \
		-e GITHUB_TOKEN="${GH_TOKEN}" \
		-e CARGO_REGISTRY_TOKEN="${CARGO_REGISTRY_TOKEN:-}" \
		-e RUSTFLAGS="-C target-cpu=native" \
		-e CARGO_BUILD_JOBS=8 \
		-e RUSTC_WRAPPER=sccache \
		-e "TERM=${TERM:-xterm-256color}" \
		--cpus=8 \
		--memory=16g \
		rivrs-linalg:latest \
		/bin/bash
fi

echo ""
echo -e "${GREEN}Exited container${NC}"
echo ""
echo "To restart: docker start rivrs-linalg && docker exec -it rivrs-linalg bash"
echo "To remove:  docker rm rivrs-linalg"
