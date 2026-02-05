# rivrs-linalg Docker Development Environment

Complete, isolated Docker development environment for rivrs-linalg with all reference materials and GitHub integration.

## Architecture

**Platform**: Native ARM64 for Apple Silicon (M1/M2/M3 Macs)

**Key Features**:
- Rust 1.93.0 with rust-analyzer, clippy, rustfmt
- Claude Code AI assistant pre-installed
- Node.js 20 + npm for Claude Code and JavaScript/TypeScript development
- Full development tool suite (fzf, zsh, vim, jq, etc.)
- All reference materials included (730MB: faer-rs, lapack, SLICOT-Reference, slicot)
- GitHub authentication via GitHub CLI token forwarding
- sccache for fast compilation
- Persistent Cargo and Claude Code config caches across container restarts
- 8 CPUs, 16GB RAM by default

## Quick Start

### Prerequisites

1. **Install Docker Desktop** (with ARM64 support for Apple Silicon)
2. **Authenticate with GitHub CLI** (one-time setup):
   ```bash
   gh auth login
   ```
   Follow the prompts to authenticate via web browser or token.

### Option 1: VS Code Dev Containers (Recommended)

1. Install "Dev Containers" extension in VS Code
2. Ensure you're logged in to GitHub: `gh auth status`
3. Open project folder in VS Code
4. Press `F1` → "Dev Containers: Reopen in Container"
5. Wait for build (8-12 minutes first time, native ARM64 speed)
6. Start coding! rust-analyzer and all extensions are pre-configured

**Authentication**: VS Code automatically forwards your GitHub CLI credentials and Git credential helper. No manual setup needed!

### Option 2: Docker Compose (Command Line)

```bash
# Build the image
cd docker
./build.sh

# Uncomment the gh config mount in .devcontainer/docker-compose.yml:
# - ${HOME}/.config/gh:/root/.config/gh:ro

# Start container in background
docker-compose up -d

# Enter container
docker exec -it rivrs-linalg-dev bash

# Inside container
rustc --version  # Verify Rust works
gh auth status   # Verify GitHub CLI access
git pull         # Test Git operations

# Stop container (preserves state)
docker-compose down
```

**Note**: Docker Compose doesn't have VS Code's automatic credential forwarding, so you need to manually mount the GitHub CLI config.

### Option 3: Pure Docker (Advanced)

```bash
# Build the image
cd docker
./build.sh

# Run container
./run.sh

# Container automatically starts with interactive shell
```

## Container Structure

```
/workspace/rivrs-linalg/              # Project root
├── faer-rs/                   # Rust linear algebra reference
├── lapack/                    # BSD-licensed LAPACK (consult freely)
├── SLICOT-Reference/          # SLICOT docs and tests ONLY
├── slicot/                    # SLICOT docs and tests ONLY
├── src/                       # Your Rust code
├── tests/                     # Tests
├── examples/                  # Examples
├── .rustfmt.toml              # Hard tabs, 80-char width
└── .claude/                   # Claude Code state

/root/.cargo/                  # Rust installation
/root/.cache/sccache/          # Compilation cache
/root/.claude/                 # Claude Code config (persisted)
```

## Claude Code Integration

**Claude Code** is pre-installed and ready to use for AI-assisted development.

### How to Use Claude Code

**In VS Code (Recommended)**:
1. Reopen the project in the container (F1 → "Dev Containers: Reopen in Container")
2. The Claude Code extension is automatically loaded
3. Click the Claude icon in the sidebar to start interacting
4. Claude can read/edit Rust files, run tests, and assist with rivrs-linalg development

**From Command Line**:
```bash
# Inside container
claude --version  # Verify installation

# Claude Code CLI is available for terminal-based AI assistance
# Note: You'll need to authenticate via the VS Code extension first
```

### What Claude Code Can Do

- **Read and analyze** Rust code in the rivrs-linalg project
- **Edit files** based on your instructions
- **Run commands** (cargo build, cargo test, etc.)
- **Explain algorithms** from reference materials (following clean room rules)
- **Suggest implementations** based on academic papers and LAPACK patterns
- **Debug issues** by analyzing error messages and code

### Clean Room Compliance

Claude Code is configured to follow rivrs-linalg clean room implementation rules:
- Will NOT read SLICOT Fortran source code (slicot/src/*.f files)
- Will use SLICOT documentation, test cases, and LAPACK source as references
- Will help implement algorithms from academic papers and textbooks
- Will document which references were used for each implementation

### Persistent Configuration

Claude Code configuration is stored in a persistent Docker volume:
- Volume name: `rivrs-linalg-claude-config`
- Mount point: `/root/.claude`
- Survives container restarts and rebuilds
- Plugins and settings are preserved

To reset Claude Code config:
```bash
docker volume rm rivrs-linalg-claude-config
# Next container start will create a fresh volume
```

## Development Workflow

### 1. Make Changes in Container

All development happens inside the container:

```bash
# Inside container
cd /workspace/rivrs-linalg

# Create/edit code
vim src/lib.rs

# Build with caching
cargo build  # sccache speeds this up

# Run tests
cargo test

# Format and lint
cargo fmt
cargo clippy
```

### 2. Commit and Push from Container

```bash
# Inside container
git add .
git commit -m "feat: add Sylvester equation solver"
git push origin 001-sylvester-solver

# SSH key forwarded from host - no password needed
```

### 3. Pull Changes on Host (if needed)

```bash
# On host machine
cd /Users/jared/Dropbox/projects/rivrs-linalg
git pull origin 001-sylvester-solver
```

### 4. Container Persistence

Container state persists until you remove it:

```bash
# Stop container (state saved)
docker-compose down

# Restart container (state restored)
docker-compose up -d
docker exec -it rivrs-linalg-dev bash

# Remove container (state lost)
docker rm rivrs-linalg-dev
```

**Important**: Uncommitted code in the container is lost if you remove the container. Always push to GitHub!

## GitHub Authentication

**Credentials are NEVER embedded in the Docker image**. They're forwarded at runtime:

### How It Works

- **VS Code**: Automatically forwards Git credentials and mounts `~/.config/gh` directory
- **Docker Compose**: Manually mount `~/.config/gh` (uncomment line in docker-compose.yml)
- **Protocol**: HTTPS with GitHub CLI tokens (not SSH keys)
- **Security**: Tokens are read-only mounted, never copied into image layers

### Setup (One-Time)

```bash
# Authenticate with GitHub CLI on your host
gh auth login

# Verify
gh auth status  # Should show logged in to github.com
```

**Inside Container**: Git and GitHub CLI operations work automatically using your host credentials.

## Reference Materials

### What You Can Use

✅ **LAPACK source** (`lapack/`): BSD-licensed, consult freely for numerical algorithms

✅ **SLICOT documentation** (`SLICOT-Reference/doc/`, `slicot/doc/`): Understand what algorithms do

✅ **SLICOT test cases** (`slicot/examples/`): Validate your implementations

✅ **faer-rs source** (`faer-rs/`): Learn modern Rust linear algebra patterns

### What You CANNOT Use

❌ **SLICOT Fortran source** (`slicot/src/*.f`): GPL-licensed, reading it contaminates clean room

**Rationale**: rivrs-linalg must remain permissively licensed (MIT/Apache-2.0). Reading GPL code during implementation creates copyright issues. Use academic papers, textbooks, and BSD-licensed LAPACK instead.

## Performance Notes

### ARM64 Optimization

The container is built natively for ARM64 (Apple Silicon):

- **RUSTFLAGS**: `-C target-cpu=native` (ARM-specific optimization)
- **SIMD**: Uses ARM NEON or portable implementations (not x86 SSE/AVX)
- **Build times**: 8-12 min first build, 1-2 min subsequent (native speed)
- **Runtime**: Excellent performance on M1/M2/M3 Macs

**Trade-off**: Won't match x86_64 SIMD characteristics in production. For x86 testing, build with:
```bash
docker build --platform linux/amd64 ...  # Slow, uses emulation
```

### Compilation Caching

**sccache** caches compilation artifacts:

```bash
# Check cache statistics
sccache --show-stats

# Clear cache if needed
sccache --stop-server
rm -rf /root/.cache/sccache/*
```

**First build** of faer-rs: 3-8 minutes (native ARM64)
**Second build**: <30 seconds (cache hit)

### Cargo Configuration

Pre-configured for performance:

```toml
[build]
jobs = 8                # Parallel compilation
incremental = true      # Incremental builds

[profile.dev]
opt-level = 3           # Match faer-rs (fast debug builds)
```

## Troubleshooting

### "Permission denied" or authentication errors with GitHub

**Cause**: GitHub CLI not authenticated on host

**Fix**:
```bash
# On host, check authentication status
gh auth status

# If not logged in, authenticate
gh auth login

# For VS Code users: restart container
# Press F1 → "Dev Containers: Rebuild Container"

# For docker-compose users: restart container
docker-compose restart
```

### "Cannot connect to Docker daemon"

**Cause**: Docker Desktop not running

**Fix**: Start Docker Desktop and wait for it to fully initialize

### Image build fails with "COPY failed"

**Cause**: Reference materials not present in project root

**Fix**: Ensure `faer-rs/`, `lapack/`, `SLICOT-Reference/`, `slicot/` exist:
```bash
ls -l faer-rs lapack SLICOT-Reference slicot
```

### Container won't start: "platform mismatch"

**Cause**: Running on non-ARM64 system

**Fix**: Remove `platform: linux/arm64` from `docker-compose.yml` (will auto-detect)

### rust-analyzer not working in VS Code

**Cause**: Extensions not installed or container rebuild needed

**Fix**:
1. Press `F1` → "Dev Containers: Rebuild Container"
2. Wait for rebuild
3. rust-analyzer should activate (check bottom-right status bar)

### sccache not caching

**Cause**: Volume not mounted or sccache server crashed

**Fix**:
```bash
# Check volume is mounted
docker inspect rivrs-linalg-dev | grep sccache

# Restart sccache
sccache --stop-server
sccache --start-server

# Check it's working
sccache --show-stats
```

### Out of disk space

**Cause**: Docker images/volumes consuming space

**Fix**:
```bash
# See disk usage
docker system df

# Clean up (careful: removes unused images)
docker system prune -a

# Remove only rivrs-linalg volumes
docker volume rm rivrs-linalg-cargo-cache rivrs-linalg-sccache-cache
```

## Clean Room Implementation Reminder

**This Docker environment includes reference materials, but you must still follow clean room rules**:

1. **Never read** `slicot/src/*.f` files (GPL source code)
2. **Do read**:
   - `slicot/doc/` (what algorithms do)
   - `slicot/examples/` (test cases)
   - `lapack/` (BSD-licensed, freely consultable)
   - `faer-rs/` (Rust patterns)
3. **Implement from**: academic papers, textbooks, LAPACK source
4. **Document**: which papers/books you used for each implementation

## Backup Strategy

Container state is ephemeral. **Always push commits to GitHub**:

```bash
# Good: Changes safe on GitHub
git add .
git commit -m "progress checkpoint"
git push origin 001-sylvester-solver

# Bad: Changes only in container
# (lost if container removed)
```

## Resource Limits

Default limits (adjust in `docker-compose.yml`):

- **CPUs**: 8 cores
- **Memory**: 16GB
- **Swap**: Unlimited (uses host swap)

For lower-end machines:

```yaml
deploy:
  resources:
    limits:
      cpus: '4'      # Reduce to 4 cores
      memory: 8G     # Reduce to 8GB
```

## Advanced Usage

### Multiple Concurrent Builds

```bash
# Terminal 1: Run tests
docker exec -it rivrs-linalg-dev bash -c "cd /workspace/rivrs-linalg && cargo test"

# Terminal 2: Check docs
docker exec -it rivrs-linalg-dev bash -c "cd /workspace/rivrs-linalg && cargo doc --open"

# Terminal 3: Format code
docker exec -it rivrs-linalg-dev bash -c "cd /workspace/rivrs-linalg && cargo fmt"
```

### Custom Docker Build

```bash
# Build with different Rust version
cd docker
docker build \
  --platform linux/arm64 \
  --build-arg RUST_VERSION=1.83.0 \
  -f Dockerfile \
  -t rivrs-linalg-dev:custom \
  ..
```

### Export Container State

```bash
# Commit container to new image
docker commit rivrs-linalg-dev rivrs-linalg-dev:checkpoint

# Export as tarball
docker save rivrs-linalg-dev:checkpoint | gzip > rivrs-linalg-dev-checkpoint.tar.gz
```

## FAQ

**Q: Do I need SSH keys for GitHub?**
A: No! This setup uses HTTPS + GitHub CLI tokens, which is simpler and more secure for dev containers.

**Q: How do I authenticate with GitHub?**
A: Run `gh auth login` on your host machine once. VS Code automatically forwards the credentials.

**Q: Do I need to commit and push every change?**
A: No, work freely in the container. Only push when you want to save work permanently or share with others.

**Q: Can I edit files on my host Mac and see changes in container?**
A: No, this is a fully isolated environment. Edit inside the container or push/pull via git.

**Q: What happens if my Mac restarts?**
A: Container stops. Run `docker-compose up -d` to restart it. Uncommitted work is preserved.

**Q: Can I use this on Linux or Intel Macs?**
A: Yes, remove `platform: linux/arm64` from `docker-compose.yml`. On Intel, it will use x86_64 instead.

**Q: How do I update Rust version?**
A: Change `RUST_VERSION` in `docker-compose.yml` and rebuild: `./build.sh`

**Q: Is the 3.3GB image size too large?**
A: It includes complete Rust toolchain + 710MB reference materials. Trade-off for total isolation.

**Q: Why not mount the project directory from host?**
A: Isolation. Prevents accidental dependency on host tools/libraries. Forces clean build.

**Q: Can I use Claude Code for Rust development?**
A: Yes! Claude Code is pre-installed and the VS Code extension is automatically loaded. Use it for AI-assisted coding, debugging, and algorithm implementation.

**Q: Why is Node.js installed in a Rust container?**
A: Claude Code is distributed as an npm package, so Node.js is required. This also enables JavaScript/TypeScript development if needed for Python bindings or web interfaces.

**Q: Does Claude Code work offline?**
A: No, Claude Code requires internet access to call the Anthropic API. However, rust-analyzer and other local tools work offline.

**Q: How much larger is the image with Claude Code?**
A: Approximately 500MB larger (3.3GB → 3.8GB) due to Node.js ecosystem and development tools. Trade-off for integrated AI assistance.

---

For issues or questions, see main project README or Claude Code documentation.
