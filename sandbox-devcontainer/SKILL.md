# Sandbox Devcontainer Builder

---
name: Sandbox Devcontainer Builder
description: This skill should be used when the user asks to "create a sandbox devcontainer", "build a sandboxed container", "set up Claude Code container", "create isolated development environment", "add Claude Code to devcontainer", or mentions creating Docker containers for Claude Code to operate in safely. Provides comprehensive guidance for building language-agnostic sandboxed devcontainers with Claude Code integration, environment isolation, and static reference material embedding.
version: 1.0.0
---

## Purpose

This skill guides the creation of fully isolated Docker devcontainers that provide sandboxed environments for Claude Code to operate in. These containers embed all necessary dependencies, static reference materials, and Claude Code integration while maintaining complete isolation from the host system. The approach is language-agnostic and adaptable to any project type (Rust, Python, JavaScript, Go, Java, etc.).

## When to Use This Skill

Use this skill when:
- Creating a new devcontainer for a project with Claude Code integration
- Converting an existing project to use a sandboxed container environment
- Setting up isolated environments for AI-assisted development
- Embedding static reference materials (documentation, libraries, datasets) into containers
- Ensuring reproducible development environments across teams

## Core Principles

### 1. Complete Isolation

Containers are fully self-contained with no dependencies on host system:
- All language runtimes and tools installed in container
- Static reference materials copied into image at build time
- No mounted directories except credentials and Claude config
- Git operations via HTTPS with credential forwarding (not SSH)

### 2. Claude Code Integration

Every sandbox includes:
- Node.js 20 (required for Claude Code npm package)
- Claude Code installed globally via npm
- VS Code extension pre-configured
- Persistent configuration volume for Claude state
- Environment variables for optimal operation

### 3. Static Reference Embedding

Non-git-tracked files and reference materials embedded in image:
- Local documentation and specifications
- Third-party libraries and frameworks (for reference)
- Datasets and test fixtures
- Example code and templates
- Anything not in version control but needed for development

### 4. Clean Room Patterns

For projects requiring license isolation (like CSRRS):
- Embed GPL/proprietary code for reference only (documentation, tests)
- Document which sources can/cannot be consulted during implementation
- Provide clear guidelines in welcome scripts and documentation

## Implementation Workflow

### Step 1: Analyze Project Structure

Examine the project to determine:

1. **Primary language/runtime**: Python, Rust, Node.js, Go, Java, etc.
2. **Build tools needed**: cargo, npm, gradle, make, cmake, etc.
3. **Static references**: What files exist locally but aren't in git?
4. **Dependencies**: OS packages, language packages, system libraries
5. **Development tools**: Debuggers, formatters, linters, test runners

Use these commands to analyze:

```bash
# Find untracked files (potential static references)
git status --short | grep '^??'
find . -type f ! -path '*/.git/*' -print0 | git check-ignore --stdin -z -v

# Identify primary languages
find . -name '*.py' -o -name '*.rs' -o -name '*.go' -o -name '*.java' | head -20

# Check for build configuration files
ls -la | grep -E '(Cargo.toml|package.json|requirements.txt|go.mod|pom.xml|build.gradle)'

# Estimate static reference sizes
du -sh */ | grep -v '.git'
```

Consult `scripts/analyze-project.sh` for automated project analysis.

### Step 2: Design Container Architecture

Use multi-stage builds for efficiency:

**Stage 1 (base-builder)**: Install language runtime, tools, Claude Code
**Stage 2 (references)**: Copy static reference materials (if any)
**Stage 3 (development)**: Final environment with everything combined

See `references/architecture-patterns.md` for detailed multi-stage patterns for different languages.

### Step 3: Create Dockerfile

Build the Dockerfile following this structure:

```dockerfile
# Stage 1: Base builder with language runtime
FROM debian:bookworm-slim AS base-builder

ARG LANGUAGE_VERSION=latest
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl git build-essential \
    # ... language-specific packages
    && rm -rf /var/lib/apt/lists/*

# Install GitHub CLI (CRITICAL - must be in Dockerfile, not just devcontainer feature)
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*

# Install primary language runtime
# (See references/runtime-installation.md for examples)

# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    npm install -g npm@latest && \
    rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools suite
RUN apt-get update && apt-get install -y \
    fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
    && rm -rf /var/lib/apt/lists/*

# Configure git (CRITICAL - set credential helper to use GitHub CLI)
RUN git config --global user.name "${GIT_USER_NAME}" && \
    git config --global user.email "${GIT_USER_EMAIL}" && \
    git config --global init.defaultBranch main && \
    git config --global credential.helper "" && \
    git config --global credential.helper "!gh auth git-credential"

# Stage 2: Reference materials (if needed)
FROM scratch AS references
COPY reference-dir-1 /references/reference-dir-1
COPY reference-dir-2 /references/reference-dir-2

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project-name

# Copy reference materials from stage 2
COPY --from=references /references ./

# Copy project files
COPY .gitignore README.md ./

# Create project structure
RUN mkdir -p src tests examples .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
    echo 'echo "Project Development Environment"' >> /usr/local/bin/welcome.sh && \
    echo 'echo "Language: $(language --version)"' >> /usr/local/bin/welcome.sh && \
    echo 'echo "Node.js: $(node --version)"' >> /usr/local/bin/welcome.sh && \
    echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
    chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

Consult `references/dockerfile-templates.md` for language-specific Dockerfile templates.

### Step 4: Create devcontainer.json

Configure VS Code Dev Containers integration:

```json
{
  "name": "Project Development Container",
  "dockerComposeFile": "docker-compose.yml",
  "service": "project-dev",
  "workspaceFolder": "/workspace/project-name",
  "shutdownAction": "stopCompose",

  "features": {
    "ghcr.io/devcontainers/features/github-cli:1": {
      "version": "latest"
    }
  },

  "customizations": {
    "vscode": {
      "extensions": [
        "anthropic.claude-code"
        // Add language-specific extensions
      ],
      "settings": {
        // Language-specific settings
      }
    }
  },

  "remoteUser": "root",

  "remoteEnv": {
    "NODE_OPTIONS": "--max-old-space-size=4096",
    "CLAUDE_CONFIG_DIR": "/root/.claude"
  },

  "postCreateCommand": "language --version && claude --version && echo '✓ Ready'",

  "mounts": [
    "source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
    "source=project-claude-config,target=/root/.claude,type=volume"
  ]
}
```

See `references/devcontainer-configs.md` for complete examples per language.

### Step 5: Create docker-compose.yml

Define the container service:

```yaml
version: '3.8'

services:
  project-dev:
    build:
      context: ..
      dockerfile: docker/Dockerfile
      args:
        LANGUAGE_VERSION: "version"
        CLAUDE_CODE_VERSION: "latest"
        GIT_USER_NAME: "${GIT_USER_NAME:-Claude Code}"
        GIT_USER_EMAIL: "${GIT_USER_EMAIL:-claude@example.com}"
    platform: linux/arm64  # or linux/amd64
    image: project-dev:latest
    container_name: project-dev

    volumes:
      - language-cache:/root/.language-cache
      - claude-config:/root/.claude

    environment:
      - CLAUDE_CONFIG_DIR=/root/.claude
      - NODE_OPTIONS=--max-old-space-size=4096

    deploy:
      resources:
        limits:
          cpus: '8'
          memory: 16G

    command: sleep infinity
    stdin_open: true
    tty: true

volumes:
  language-cache:
    name: project-language-cache
  claude-config:
    name: project-claude-config
```

### Step 6: Create Build Script

Provide a convenient build script (`docker/build.sh`):

```bash
#!/bin/bash
set -euo pipefail

echo "Building Development Container..."

# Extract git config from host
GIT_USER_NAME=$(git config user.name || echo "Claude Code")
GIT_USER_EMAIL=$(git config user.email || echo "claude@example.com")

# Enable BuildKit
export DOCKER_BUILDKIT=1

cd ..
docker build \
  --platform linux/arm64 \
  --build-arg LANGUAGE_VERSION=version \
  --build-arg CLAUDE_CODE_VERSION=latest \
  --build-arg GIT_USER_NAME="$GIT_USER_NAME" \
  --build-arg GIT_USER_EMAIL="$GIT_USER_EMAIL" \
  -f docker/Dockerfile \
  -t project-dev:latest \
  --progress=plain \
  .

echo "Build complete!"
echo "Next: Open in VS Code and 'Reopen in Container'"
```

### Step 7: Create .dockerignore

Exclude unnecessary files from build context:

```
# Git internals
.git/objects/
.git/logs/

# Build artifacts
target/
dist/
build/
*.o
*.so
*.dylib

# IDE files
.vscode/
.idea/
*.swp

# Claude state
.claude/

# Docker files (avoid recursion)
.devcontainer/
docker/
Dockerfile
docker-compose.yml

# Dependencies (will be installed in container)
node_modules/
venv/
__pycache__/

# Keep static references - they must be copied!
# reference-dir-1/
# reference-dir-2/
```

### Step 8: Document the Environment

Create comprehensive documentation (`docker/README.md`):

1. **Architecture**: Platform, key features, container structure
2. **Quick Start**: Prerequisites, VS Code workflow, CLI workflow
3. **Claude Code Integration**: How to use, what it can do, config persistence
4. **Development Workflow**: Making changes, committing, container persistence
5. **GitHub Authentication**: How credentials are forwarded (HTTPS + gh CLI)
6. **Reference Materials**: What can/cannot be used (if applicable)
7. **Troubleshooting**: Common issues and solutions
8. **FAQ**: Frequently asked questions

See `examples/README-template.md` for a complete template.

## GitHub Authentication Pattern

**Critical**: Use HTTPS with GitHub CLI token forwarding, NOT SSH keys.

### Why HTTPS Over SSH

1. **VS Code Integration**: Automatic credential forwarding via `gh` CLI
2. **Simplicity**: No SSH key management in containers
3. **Security**: Tokens are read-only mounted, never embedded in images
4. **Compatibility**: Works seamlessly with Dev Containers feature

### Implementation

1. **Host Setup** (one-time):
   ```bash
   gh auth login
   gh auth status  # Verify
   ```

2. **Devcontainer Config**:
   ```json
   "features": {
     "ghcr.io/devcontainers/features/github-cli:1": {"version": "latest"}
   },
   "mounts": [
     "source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached"
   ]
   ```

3. **Inside Container**: Git and `gh` operations work automatically using host credentials

See `references/github-auth-patterns.md` for detailed authentication setup.

## Static Reference Materials

### Identifying Static References

Static references are files needed for development but not in version control:

```bash
# Find untracked files
git status --short | grep '^??'

# Find ignored files that exist
git clean -ndX

# Common patterns
# - Third-party libraries (for reference, not linking)
# - Documentation archives
# - Test datasets
# - Example code from other projects
# - Proprietary/GPL code (reference only)
```

### Embedding Strategy

**Multi-stage build with dedicated reference stage**:

```dockerfile
# Stage 2: References
FROM scratch AS references
COPY large-reference-lib/ /references/large-reference-lib/
COPY documentation/ /references/documentation/
COPY test-datasets/ /references/test-datasets/

# Stage 3: Development
FROM base-builder AS development
COPY --from=references /references ./
```

**Benefits**:
- Efficient layer caching
- Clear separation of concerns
- Easy to update references independently

### Size Considerations

Large reference materials increase image size:

```bash
# Check reference sizes before embedding
du -sh reference-dir-1 reference-dir-2

# Expected impact
# 100MB references → ~3.5GB total image
# 500MB references → ~3.9GB total image
# 1GB+ references → Consider multi-container or volume approach
```

Consult `references/reference-strategies.md` for strategies with very large references.

## Language-Specific Patterns

### Python Projects

**Runtime**: Python 3.x via apt or pyenv
**Cache volumes**: pip cache, venv
**Extensions**: Python, Pylance, Python Debugger
**Config**: PEP 8 formatting, type checking

### Rust Projects

**Runtime**: rustup with specific toolchain version
**Cache volumes**: Cargo registry, sccache
**Extensions**: rust-analyzer, LLDB
**Config**: Hard tabs, rustfmt, clippy

### Node.js/TypeScript Projects

**Runtime**: Node.js via NodeSource (same as Claude Code)
**Cache volumes**: npm cache, node_modules
**Extensions**: ESLint, Prettier, TypeScript
**Config**: Formatting, linting rules

### Go Projects

**Runtime**: Go via official install script
**Cache volumes**: Go module cache
**Extensions**: Go extension
**Config**: gofmt, golangci-lint

### Java Projects

**Runtime**: OpenJDK via apt
**Cache volumes**: Maven/Gradle cache
**Extensions**: Java Extension Pack
**Config**: Formatter, build tool integration

Detailed patterns in `references/language-runtimes.md`.

## Critical Implementation Patterns (Battle-Tested)

### 1. GitHub CLI Installation

**MUST install in Dockerfile**, not rely on devcontainer feature alone:

```dockerfile
# Install GitHub CLI
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*
```

**Why**: Ensures `gh` is available for entrypoint scripts and standalone Docker usage.

### 2. Git Credential Helper Configuration

**MUST explicitly configure** git to use GitHub CLI:

```dockerfile
RUN git config --global credential.helper "" && \
    git config --global credential.helper "!gh auth git-credential"
```

**Why**: Git doesn't automatically know to use `gh` for authentication. This makes all git operations use GitHub CLI tokens.

### 3. macOS Keychain Workaround

**Problem**: macOS stores `gh` tokens in Keychain, not in `~/.config/gh/hosts.yml` files. Containers can't access Keychain.

**Solution**: Pass token via environment variable:

**run.sh**:
```bash
GH_TOKEN=$(gh auth token 2>/dev/null || echo "")
docker run -e GH_TOKEN="${GH_TOKEN}" -e GITHUB_TOKEN="${GH_TOKEN}" ...
```

**docker-compose.yml**:
```yaml
environment:
  - GH_TOKEN=${GH_TOKEN:-}
  - GITHUB_TOKEN=${GITHUB_TOKEN:-}
```

**devcontainer.json**:
```json
{
  "remoteEnv": {
    "GH_TOKEN": "${localEnv:GH_TOKEN}",
    "GITHUB_TOKEN": "${localEnv:GITHUB_TOKEN}"
  }
}
```

**Why**: `gh` CLI automatically reads from `GH_TOKEN` environment variable, bypassing Keychain/file storage issues.

### 4. Private Repository Cloning (Runtime, Not Build-Time)

For private repos, **clone in entrypoint script** (not during docker build):

```dockerfile
# Create entrypoint that clones repo using GH_TOKEN
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'if [ ! -d "/workspace/project/.git" ]; then' >> /usr/local/bin/entrypoint.sh && \
    echo '  gh repo clone username/project /workspace/project || exit 1' >> /usr/local/bin/entrypoint.sh && \
    echo '  cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config --local safe.directory /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"' >> /usr/local/bin/entrypoint.sh && \
    echo '  git fetch --all --quiet' >> /usr/local/bin/entrypoint.sh && \
    echo 'fi' >> /usr/local/bin/entrypoint.sh && \
    echo 'cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo 'exec "$@"' >> /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
```

**Why**:
- Build-time clone fails (no auth credentials during `docker build`)
- Entrypoint runs after container starts with `GH_TOKEN` available
- Automatically sets up repository on first container start

### 5. Fetch Refspec Fix

**Critical**: `gh repo clone` sometimes fails to set proper fetch configuration. **MUST fix in entrypoint**:

```bash
git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all --quiet
```

**Symptom without fix**: `git branch -a` shows only current branch, can't checkout other branches
**Why**: GitHub CLI clone bug doesn't set fetch refspec properly

### 6. Symlink .gitignore Patterns

When symlinking reference materials, use patterns **without trailing slash**:

**Correct**:
```gitignore
faer-rs
lapack
slicot
```

**Incorrect** (won't work):
```gitignore
faer-rs/
lapack/
slicot/
```

**Why**: Git treats symlinks differently than directories. Trailing `/` only matches directories, not symlinks.

### 7. Isolated Workspace (Volume, Not Mount)

For true sandboxing, use Docker **volume**, NOT host mount:

**Correct** (isolated):
```yaml
volumes:
  - workspace:/workspace/project  # Docker volume
```

**Incorrect** (breaks sandbox):
```yaml
volumes:
  - ..:/workspace/project  # Host mount - exposes host files!
```

**run.sh correct**:
```bash
docker run -v project-workspace:/workspace/project ...
```

**run.sh incorrect**:
```bash
docker run -v "$(pwd)":/workspace/project ...  # Breaks isolation!
```

**Why**: Volumes provide complete isolation. Mounts give container/AI access to host filesystem.

### Summary of Critical Patterns

1. ✅ Install GitHub CLI in Dockerfile
2. ✅ Configure git credential helper explicitly
3. ✅ Pass GH_TOKEN via environment variable (macOS Keychain workaround)
4. ✅ Clone private repos at runtime via entrypoint
5. ✅ Fix fetch refspec after gh clone
6. ✅ Use .gitignore patterns without `/` for symlinks
7. ✅ Use Docker volumes, not mounts, for isolation

**See**: `references/github-auth-patterns.md` for complete implementation details.

## Verification Checklist

After creating the devcontainer:

**Build Verification**:
- [ ] `docker/build.sh` runs successfully
- [ ] Image builds in reasonable time (10-20 min first build)
- [ ] Final image size is acceptable (<5GB typically)

**Claude Code Functionality**:
- [ ] `claude --version` works inside container
- [ ] VS Code extension loads automatically
- [ ] Claude can read/edit project files
- [ ] Config persists across container restarts

**Language Runtime**:
- [ ] Primary language tools work (`cargo build`, `python`, etc.)
- [ ] Package managers work (`cargo`, `pip`, `npm`, etc.)
- [ ] Debugger and linter extensions work in VS Code

**GitHub Operations**:
- [ ] `gh auth status` shows authenticated
- [ ] `git clone` works via HTTPS
- [ ] `git push` works without password prompt

**Reference Materials**:
- [ ] All static references present in container
- [ ] Sizes match expectations (`du -sh reference-dirs`)
- [ ] Documentation accessible

**Development Workflow**:
- [ ] Can create/edit files in container
- [ ] Can run tests and builds
- [ ] Changes persist until container removed
- [ ] Git operations work smoothly

## Common Patterns

### Pattern: Minimal Sandbox (No References)

Simple project with just language runtime + Claude Code:

```
project/
├── .devcontainer/
│   ├── devcontainer.json
│   └── docker-compose.yml
└── docker/
    ├── Dockerfile (2 stages: base-builder → development)
    ├── build.sh
    └── .dockerignore
```

**Use case**: New projects, minimal dependencies, no local references

### Pattern: Reference-Heavy Sandbox (Like CSRRS)

Complex project with large reference materials:

```
project/
├── .devcontainer/
│   ├── devcontainer.json
│   └── docker-compose.yml
├── docker/
│   ├── Dockerfile (3 stages: base-builder → references → development)
│   ├── build.sh
│   ├── .dockerignore
│   └── README.md
├── reference-lib-1/ (500MB, not in git)
├── reference-lib-2/ (200MB, not in git)
└── documentation/ (100MB, not in git)
```

**Use case**: Projects with large local references, clean room implementations, embedded documentation

### Pattern: Multi-Language Sandbox

Project using multiple languages:

```dockerfile
# Install Python
RUN apt-get install -y python3 python3-pip

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Install Node.js (for Claude Code + project needs)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs
```

**Use case**: Full-stack projects, Python bindings for native code, polyglot systems

## Additional Resources

### Reference Files

For detailed patterns and language-specific configurations:
- **`references/architecture-patterns.md`** - Multi-stage build patterns
- **`references/dockerfile-templates.md`** - Complete Dockerfile examples per language
- **`references/devcontainer-configs.md`** - VS Code devcontainer.json examples
- **`references/language-runtimes.md`** - Runtime installation for Python, Rust, Go, Java, etc.
- **`references/github-auth-patterns.md`** - Detailed GitHub authentication setup
- **`references/reference-strategies.md`** - Handling large reference materials
- **`references/anthropic-reference.md`** - Official Anthropic devcontainer implementation

### Example Files

Complete working examples:
- **`examples/python-project/`** - Python data science sandbox
- **`examples/rust-project/`** - Rust systems programming sandbox (CSRRS)
- **`examples/nodejs-project/`** - Node.js web development sandbox
- **`examples/README-template.md`** - Docker README documentation template

### Scripts

Utility scripts for automation:
- **`scripts/analyze-project.sh`** - Automated project structure analysis
- **`scripts/generate-dockerfile.sh`** - Generate Dockerfile from project analysis
- **`scripts/test-container.sh`** - Verify container functionality

## Best Practices

### DO:
- Use multi-stage builds for efficiency
- Install Claude Code in base-builder stage
- Mount GitHub credentials, never embed them
- Use named volumes for caches (language package manager, Claude config)
- Create comprehensive README documentation
- Test the container before finalizing
- Use BuildKit for faster builds (`export DOCKER_BUILDKIT=1`)
- Set appropriate resource limits (CPU, memory)

### DON'T:
- Embed credentials or API keys in images
- Mount project directory from host (breaks isolation)
- Use SSH for git operations (use HTTPS + gh CLI)
- Skip .dockerignore (wastes time copying unnecessary files)
- Forget to clean apt lists after installs (`rm -rf /var/lib/apt/lists/*`)
- Make images unnecessarily large (combine RUN commands, clean caches)
- Use vague welcome messages (show versions and paths clearly)

## Troubleshooting

**Build fails with "platform mismatch"**:
- Remove or adjust `platform: linux/arm64` in docker-compose.yml
- Use auto-detection or specify your architecture

**Claude Code not found in container**:
- Verify Node.js installation succeeded
- Check npm install command ran (`RUN npm install -g @anthropic-ai/claude-code@latest`)
- Ensure PATH includes npm global bin (`/usr/local/bin` typically)

**Git operations require password**:
- Ensure `gh auth login` ran on host
- Verify mount in devcontainer.json: `~/.config/gh` → `/root/.config/gh`
- Check `gh auth status` inside container
- GitHub CLI feature must be enabled in devcontainer.json

**Reference materials not in container**:
- Check they're not excluded by .dockerignore
- Verify COPY commands in Dockerfile reference correct paths
- Ensure paths are relative to build context (typically project root)

**VS Code extensions not loading**:
- Rebuild container (F1 → "Dev Containers: Rebuild Container")
- Check extensions list in devcontainer.json
- Ensure Claude Code extension ID is correct: `anthropic.claude-code`

**Container consumes too much disk space**:
- Check image size: `docker images project-dev`
- Review reference material sizes: `du -sh reference-dirs`
- Consider cleaning intermediate layers, combining RUN commands
- For >5GB images, evaluate if all references are necessary

## Summary

Creating a sandboxed devcontainer involves:

1. **Analyze project**: Determine language, tools, static references
2. **Design architecture**: Plan multi-stage Dockerfile
3. **Create Dockerfile**: Base-builder → references → development
4. **Configure devcontainer.json**: Extensions, mounts, environment
5. **Define docker-compose.yml**: Service, volumes, resources
6. **Write build script**: Convenient build automation
7. **Add .dockerignore**: Exclude unnecessary files
8. **Document**: Comprehensive README for users
9. **Verify**: Test all functionality (build, Claude Code, git, language tools)
10. **Iterate**: Refine based on usage

The result is a fully isolated, reproducible development environment where Claude Code can operate safely with complete access to project files and references.
