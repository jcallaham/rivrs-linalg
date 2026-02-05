#!/bin/bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}rivrs-linalg Docker Build Script (ARM64)${NC}"
echo -e "${GREEN}=====================================${NC}"
echo ""

# Check Docker is installed
if ! command -v docker &> /dev/null; then
	echo -e "${RED}Error: Docker is not installed${NC}"
	exit 1
fi

echo -e "${YELLOW}Docker version:${NC}"
docker --version
echo ""

# Enable BuildKit for better performance
export DOCKER_BUILDKIT=1

# Display reference material sizes
echo -e "${YELLOW}Reference material sizes:${NC}"
du -sh ../references/* 2>/dev/null || echo "Note: Some reference directories may not exist yet"
echo ""

# Extract git config from host
GIT_USER_NAME=$(git config user.name || echo "Claude Code")
GIT_USER_EMAIL=$(git config user.email || echo "claude@example.com")

echo -e "${YELLOW}Building with git config:${NC}"
echo "  Name:  $GIT_USER_NAME"
echo "  Email: $GIT_USER_EMAIL"
echo ""

# Build the image
echo -e "${GREEN}Building Docker image (ARM64)...${NC}"
echo "This will take 10-15 minutes on first build (native ARM64 + Node.js)"
echo ""
echo -e "${YELLOW}Build context: repository root${NC}"
echo -e "${YELLOW}Dockerfile: docker/Dockerfile${NC}"
echo ""

# Build from repository root so COPY references/ works
cd ..
docker build \
	--platform linux/arm64 \
	--build-arg RUST_VERSION=1.93.0 \
	--build-arg CLAUDE_CODE_VERSION=latest \
	--build-arg GIT_USER_NAME="$GIT_USER_NAME" \
	--build-arg GIT_USER_EMAIL="$GIT_USER_EMAIL" \
	-f docker/Dockerfile \
	-t rivrs-linalg-dev:latest \
	--progress=plain \
	.

echo ""
echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}Build complete!${NC}"
echo -e "${GREEN}=====================================${NC}"
echo ""

# Show image size
IMAGE_SIZE=$(docker images rivrs-linalg-dev:latest --format "{{.Size}}")
echo -e "${YELLOW}Image size:${NC} $IMAGE_SIZE"
echo ""

echo -e "${YELLOW}Next steps:${NC}"
echo ""
echo "1. Authenticate with GitHub (if not already done):"
echo "   gh auth login"
echo "   gh auth status  # Verify"
echo ""
echo "2. Choose your workflow:"
echo ""
echo "   VS Code (recommended):"
echo "   - Open project in VS Code"
echo "   - Press F1 → 'Dev Containers: Reopen in Container'"
echo "   - GitHub credentials forwarded automatically!"
echo "   - Claude Code extension will be loaded automatically"
echo ""
echo "   Command line (docker):"
echo "   - cd docker && ./run.sh"
echo ""
echo "3. Claude Code is pre-installed:"
echo "   - Run 'claude --version' inside the container to verify"
echo "   - Use the Claude Code VS Code extension for AI assistance"
echo "   - Config stored in persistent volume: rivrs-linalg-claude-config"
echo ""
echo -e "${GREEN}Happy coding!${NC}"
