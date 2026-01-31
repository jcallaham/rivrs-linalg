#!/bin/bash
# analyze-project.sh - Automated project structure analysis for sandbox devcontainer creation

set -euo pipefail

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}Sandbox Devcontainer Project Analysis${NC}"
echo -e "${GREEN}=====================================${NC}"
echo ""

# Function to detect primary language
detect_language() {
	echo -e "${BLUE}Detecting primary language...${NC}"

	# Check for language-specific files
	if [ -f "Cargo.toml" ] || [ -f "Cargo.lock" ]; then
		echo "Rust"
		return 0
	fi

	if [ -f "package.json" ]; then
		echo "Node.js/JavaScript/TypeScript"
		return 0
	fi

	if [ -f "requirements.txt" ] || [ -f "setup.py" ] || [ -f "pyproject.toml" ]; then
		echo "Python"
		return 0
	fi

	if [ -f "go.mod" ] || [ -f "go.sum" ]; then
		echo "Go"
		return 0
	fi

	if [ -f "pom.xml" ] || [ -f "build.gradle" ]; then
		echo "Java"
		return 0
	fi

	if [ -f "Gemfile" ]; then
		echo "Ruby"
		return 0
	fi

	# Count file extensions
	local rust_count=$(find . -name "*.rs" 2>/dev/null | wc -l)
	local python_count=$(find . -name "*.py" 2>/dev/null | wc -l)
	local js_count=$(find . -name "*.js" -o -name "*.ts" 2>/dev/null | wc -l)
	local go_count=$(find . -name "*.go" 2>/dev/null | wc -l)
	local java_count=$(find . -name "*.java" 2>/dev/null | wc -l)

	# Determine primary based on counts
	local max_count=0
	local max_lang="Unknown"

	if [ "$rust_count" -gt "$max_count" ]; then
		max_count=$rust_count
		max_lang="Rust"
	fi

	if [ "$python_count" -gt "$max_count" ]; then
		max_count=$python_count
		max_lang="Python"
	fi

	if [ "$js_count" -gt "$max_count" ]; then
		max_count=$js_count
		max_lang="JavaScript/TypeScript"
	fi

	if [ "$go_count" -gt "$max_count" ]; then
		max_count=$go_count
		max_lang="Go"
	fi

	if [ "$java_count" -gt "$max_count" ]; then
		max_count=$java_count
		max_lang="Java"
	fi

	echo "$max_lang"
}

# Function to find untracked files (potential static references)
find_static_references() {
	echo -e "${BLUE}Finding static reference materials (untracked files)...${NC}"

	# Check if git repository
	if ! git rev-parse --is-inside-work-tree &>/dev/null; then
		echo -e "${YELLOW}Not a git repository - cannot identify untracked files${NC}"
		return
	fi

	# Find untracked files
	local untracked=$(git status --short | grep '^??' | awk '{print $2}')

	if [ -z "$untracked" ]; then
		echo -e "${GREEN}No untracked files found${NC}"
		return
	fi

	echo ""
	echo "Untracked files/directories:"
	echo "$untracked" | while read -r item; do
		if [ -d "$item" ]; then
			local size=$(du -sh "$item" 2>/dev/null | awk '{print $1}')
			echo -e "  ${YELLOW}[DIR]${NC}  $item (Size: $size)"
		else
			local size=$(du -h "$item" 2>/dev/null | awk '{print $1}')
			echo -e "  ${GREEN}[FILE]${NC} $item (Size: $size)"
		fi
	done

	echo ""
	echo "Total size of untracked files:"
	du -sh $(echo "$untracked") 2>/dev/null | awk '{print $1}' | paste -sd+ - | bc 2>/dev/null || echo "Unable to calculate"
}

# Function to estimate total size
estimate_sizes() {
	echo -e "${BLUE}Estimating container size components...${NC}"
	echo ""

	# Base image size (Debian bookworm-slim)
	echo "  Base OS (Debian bookworm-slim): ~150MB"

	# Language runtime estimates
	local lang=$1
	case "$lang" in
		"Rust")
			echo "  Rust toolchain: ~2GB"
			echo "  sccache + tools: ~100MB"
			;;
		"Python")
			echo "  Python runtime: ~200MB"
			echo "  Common packages: ~500MB-1GB"
			;;
		"Node.js/JavaScript/TypeScript")
			echo "  Node.js runtime: ~300MB"
			echo "  node_modules: varies (check package.json)"
			;;
		"Go")
			echo "  Go runtime: ~500MB"
			echo "  Go modules: varies"
			;;
		"Java")
			echo "  OpenJDK: ~400MB"
			echo "  Maven dependencies: varies"
			;;
	esac

	# Claude Code
	echo "  Node.js 20 (for Claude Code): ~300MB"
	echo "  Claude Code npm package: ~50MB"
	echo "  Development tools: ~100MB"

	# Static references
	local untracked=$(git status --short 2>/dev/null | grep '^??' | awk '{print $2}')
	if [ -n "$untracked" ]; then
		local ref_size=$(du -sh $(echo "$untracked") 2>/dev/null | awk '{print $1}' | paste -sd+ - | bc 2>/dev/null || echo "unknown")
		echo "  Static references: $ref_size"
	fi

	echo ""
	echo "Expected total image size: 3-5GB depending on dependencies"
}

# Function to check for build tools
check_build_tools() {
	echo -e "${BLUE}Checking for build configuration files...${NC}"
	echo ""

	local found=0

	# Rust
	if [ -f "Cargo.toml" ]; then
		echo -e "  ${GREEN}✓${NC} Cargo.toml (Rust)"
		found=1
	fi

	# Python
	if [ -f "requirements.txt" ]; then
		echo -e "  ${GREEN}✓${NC} requirements.txt (Python pip)"
		found=1
	fi
	if [ -f "setup.py" ]; then
		echo -e "  ${GREEN}✓${NC} setup.py (Python setuptools)"
		found=1
	fi
	if [ -f "pyproject.toml" ]; then
		echo -e "  ${GREEN}✓${NC} pyproject.toml (Python)"
		found=1
	fi

	# Node.js
	if [ -f "package.json" ]; then
		echo -e "  ${GREEN}✓${NC} package.json (Node.js/npm)"
		found=1
	fi

	# Go
	if [ -f "go.mod" ]; then
		echo -e "  ${GREEN}✓${NC} go.mod (Go modules)"
		found=1
	fi

	# Java
	if [ -f "pom.xml" ]; then
		echo -e "  ${GREEN}✓${NC} pom.xml (Maven)"
		found=1
	fi
	if [ -f "build.gradle" ]; then
		echo -e "  ${GREEN}✓${NC} build.gradle (Gradle)"
		found=1
	fi

	# Make
	if [ -f "Makefile" ]; then
		echo -e "  ${GREEN}✓${NC} Makefile"
		found=1
	fi

	# CMake
	if [ -f "CMakeLists.txt" ]; then
		echo -e "  ${GREEN}✓${NC} CMakeLists.txt (CMake)"
		found=1
	fi

	if [ $found -eq 0 ]; then
		echo -e "  ${YELLOW}No standard build configuration files found${NC}"
	fi
}

# Function to suggest devcontainer structure
suggest_structure() {
	local lang=$1

	echo -e "${BLUE}Suggested devcontainer structure:${NC}"
	echo ""

	echo "  .devcontainer/"
	echo "    ├── devcontainer.json"
	echo "    └── docker-compose.yml"
	echo "  docker/"
	echo "    ├── Dockerfile"
	echo "    ├── build.sh"
	echo "    ├── .dockerignore"
	echo "    └── README.md"

	local untracked=$(git status --short 2>/dev/null | grep '^??' | awk '{print $2}')
	if [ -n "$untracked" ]; then
		echo ""
		echo "  Static references to include in Dockerfile:"
		echo "$untracked" | while read -r item; do
			if [ -d "$item" ]; then
				echo "    ├── $item/ (copy into container)"
			fi
		done
	fi
}

# Function to generate next steps
generate_next_steps() {
	local lang=$1

	echo -e "${BLUE}Next Steps:${NC}"
	echo ""
	echo "1. Create directory structure:"
	echo "   mkdir -p .devcontainer docker"
	echo ""
	echo "2. Choose appropriate Dockerfile template for $lang"
	echo "   See: sandbox-devcontainer/references/dockerfile-templates.md"
	echo ""
	echo "3. Create devcontainer.json"
	echo "   See: sandbox-devcontainer/references/devcontainer-configs.md"
	echo ""
	echo "4. Create docker-compose.yml with:"
	echo "   - Build args for language version"
	echo "   - Volume mounts for caches"
	echo "   - Claude config volume"
	echo ""
	echo "5. Set up GitHub authentication:"
	echo "   gh auth login"
	echo "   See: sandbox-devcontainer/references/github-auth-patterns.md"
	echo ""
	echo "6. Build container:"
	echo "   cd docker && ./build.sh"
	echo ""
	echo "7. Open in VS Code:"
	echo "   F1 → 'Dev Containers: Reopen in Container'"
}

# Main execution
main() {
	# Detect language
	LANGUAGE=$(detect_language)
	echo -e "Primary language: ${GREEN}$LANGUAGE${NC}"
	echo ""

	# Check build tools
	check_build_tools
	echo ""

	# Find static references
	find_static_references
	echo ""

	# Estimate sizes
	estimate_sizes "$LANGUAGE"
	echo ""

	# Suggest structure
	suggest_structure "$LANGUAGE"
	echo ""

	# Generate next steps
	generate_next_steps "$LANGUAGE"
	echo ""

	echo -e "${GREEN}=====================================${NC}"
	echo -e "${GREEN}Analysis Complete${NC}"
	echo -e "${GREEN}=====================================${NC}"
}

# Run main
main
