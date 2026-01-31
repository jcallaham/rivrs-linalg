# [Project Name] Docker Development Environment

Complete, isolated Docker development environment for [Project Name] with Claude Code AI assistance and all dependencies.

## Architecture

**Platform**: [linux/arm64 | linux/amd64 | multi-platform]

**Key Features**:
- [Primary Language] [version] with [key tools]
- Claude Code AI assistant pre-installed
- Node.js 20 + npm for Claude Code and [optional: project needs]
- Full development tool suite (fzf, zsh, vim, jq, etc.)
- [Optional: Static reference materials included (size: XMB)]
- GitHub authentication via GitHub CLI token forwarding
- [Build cache tool] for fast compilation
- Persistent caches across container restarts
- [X] CPUs, [X]GB RAM by default

## Quick Start

### Prerequisites

1. **Install Docker Desktop** ([with ARM64/AMD64 support])
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
5. Wait for build ([X]-[Y] minutes first time)
6. Start coding! Claude Code and all extensions are pre-configured

**Authentication**: VS Code automatically forwards your GitHub CLI credentials. No manual setup needed!

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
docker exec -it [container-name] bash

# Inside container
[language-command] --version  # Verify language works
gh auth status                 # Verify GitHub CLI access
git pull                       # Test Git operations

# Stop container (preserves state)
docker-compose down
```

**Note**: Docker Compose doesn't have VS Code's automatic credential forwarding, so you need to manually mount the GitHub CLI config.

## Container Structure

```
/workspace/[project-name]/     # Project root
├── [reference-dir-1]/         # [Description of reference materials]
├── [reference-dir-2]/         # [Description of reference materials]
├── src/                       # Source code
├── tests/                     # Tests
├── [other-dirs]/              # [Description]
└── .claude/                   # Claude Code state

/root/.[language-cache]/       # Language package cache
/root/.claude/                 # Claude Code config (persisted)
```

## Claude Code Integration

**Claude Code** is pre-installed and ready to use for AI-assisted development.

### How to Use Claude Code

**In VS Code (Recommended)**:
1. Reopen the project in the container (F1 → "Dev Containers: Reopen in Container")
2. The Claude Code extension is automatically loaded
3. Click the Claude icon in the sidebar to start interacting
4. Claude can read/edit files, run tests, and assist with development

**From Command Line**:
```bash
# Inside container
claude --version  # Verify installation

# Claude Code CLI is available for terminal-based AI assistance
# Note: You'll need to authenticate via the VS Code extension first
```

### What Claude Code Can Do

- **Read and analyze** code in the project
- **Edit files** based on your instructions
- **Run commands** ([build-command], [test-command], etc.)
- **Explain** complex code and algorithms
- **Suggest implementations** based on best practices
- **Debug issues** by analyzing error messages and code

[Optional section for clean room projects:]
### Clean Room Compliance

Claude Code is configured to follow [project] clean room implementation rules:
- Will NOT read [restricted files/directories]
- Will use [allowed references] as permitted sources
- Will help implement algorithms from [academic papers, textbooks, etc.]
- Will document which references were used for each implementation

### Persistent Configuration

Claude Code configuration is stored in a persistent Docker volume:
- Volume name: `[project]-claude-config`
- Mount point: `/root/.claude`
- Survives container restarts and rebuilds
- Plugins and settings are preserved

To reset Claude Code config:
```bash
docker volume rm [project]-claude-config
# Next container start will create a fresh volume
```

## Development Workflow

### 1. Make Changes in Container

All development happens inside the container:

```bash
# Inside container
cd /workspace/[project-name]

# Create/edit code
vim src/[file]

# Build
[build-command]

# Run tests
[test-command]

# Format and lint
[format-command]
[lint-command]
```

### 2. Commit and Push from Container

```bash
# Inside container
git add .
git commit -m "feat: add new feature"
git push origin [branch-name]

# GitHub CLI credentials forwarded from host - no password needed
```

### 3. Pull Changes on Host (if needed)

```bash
# On host machine
cd /path/to/project
git pull origin [branch-name]
```

### 4. Container Persistence

Container state persists until you remove it:

```bash
# Stop container (state saved)
docker-compose down

# Restart container (state restored)
docker-compose up -d
docker exec -it [container-name] bash

# Remove container (state lost)
docker rm [container-name]
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

[Optional section for projects with reference materials:]
## Reference Materials

### What You Can Use

✅ **[Reference 1]** (`[path]/`): [Description, license, usage]

✅ **[Reference 2]** (`[path]/`): [Description, license, usage]

### What You CANNOT Use

❌ **[Restricted materials]** (`[path]/`): [Why restricted, licensing issues]

**Rationale**: [Explanation of licensing strategy or clean room requirements]

## Performance Notes

### [Architecture] Optimization

The container is built [natively for ARM64 | for AMD64 | multi-platform]:

- **[Compiler flags]**: [Optimization settings]
- **Build times**: [X]-[Y] min first build, [X]-[Y] min subsequent
- **Runtime**: [Performance characteristics]

[If cross-compiling:]
**Trade-off**: [Explanation of architecture-specific considerations]

### [Build Cache] Caching

**[Cache tool]** caches compilation artifacts:

```bash
# Check cache statistics
[cache-stats-command]

# Clear cache if needed
[cache-clear-command]
```

**First build**: [X]-[Y] minutes
**Second build**: <[X] seconds (cache hit)

### [Language] Configuration

Pre-configured for performance:

```[toml|yaml|json]
[configuration example]
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

**Cause**: [Reference materials | dependencies] not present in project root

**Fix**: Ensure required files/directories exist:
```bash
ls -l [expected-files-or-dirs]
```

### Container won't start: "platform mismatch"

**Cause**: Running on non-[specified] system

**Fix**: Remove `platform: linux/[arch]` from `docker-compose.yml` (will auto-detect)

### [Language-specific extension] not working in VS Code

**Cause**: Extensions not installed or container rebuild needed

**Fix**:
1. Press `F1` → "Dev Containers: Rebuild Container"
2. Wait for rebuild
3. Extension should activate (check bottom-right status bar)

### [Build cache] not caching

**Cause**: Volume not mounted or cache server crashed

**Fix**:
```bash
# Check volume is mounted
docker inspect [container-name] | grep [cache-name]

# Restart cache server
[cache-restart-command]

# Check it's working
[cache-stats-command]
```

### Out of disk space

**Cause**: Docker images/volumes consuming space

**Fix**:
```bash
# See disk usage
docker system df

# Clean up (careful: removes unused images)
docker system prune -a

# Remove only [project] volumes
docker volume rm [project]-[cache-1] [project]-[cache-2]
```

[Optional section for clean room projects:]
## Clean Room Implementation Reminder

**This Docker environment includes reference materials, but you must still follow clean room rules**:

1. **Never read** `[restricted-paths]` (GPL source code or proprietary code)
2. **Do read**:
   - `[allowed-path-1]` ([usage description])
   - `[allowed-path-2]` ([usage description])
3. **Implement from**: academic papers, textbooks, [BSD-licensed references]
4. **Document**: which papers/books you used for each implementation

## Backup Strategy

Container state is ephemeral. **Always push commits to GitHub**:

```bash
# Good: Changes safe on GitHub
git add .
git commit -m "progress checkpoint"
git push origin [branch-name]

# Bad: Changes only in container
# (lost if container removed)
```

## Resource Limits

Default limits (adjust in `docker-compose.yml`):

- **CPUs**: [X] cores
- **Memory**: [X]GB
- **Swap**: Unlimited (uses host swap)

For lower-end machines:

```yaml
deploy:
  resources:
    limits:
      cpus: '[reduced-value]'
      memory: [reduced-value]G
```

## Advanced Usage

### Multiple Concurrent Commands

```bash
# Terminal 1: Run tests
docker exec -it [container] bash -c "cd /workspace/[project] && [test-command]"

# Terminal 2: Check docs
docker exec -it [container] bash -c "cd /workspace/[project] && [docs-command]"

# Terminal 3: Format code
docker exec -it [container] bash -c "cd /workspace/[project] && [format-command]"
```

### Custom Docker Build

```bash
# Build with different [language] version
cd docker
docker build \
  --platform linux/[arch] \
  --build-arg [LANGUAGE]_VERSION=[version] \
  -f Dockerfile \
  -t [project-dev]:custom \
  ..
```

### Export Container State

```bash
# Commit container to new image
docker commit [container-name] [project-dev]:checkpoint

# Export as tarball
docker save [project-dev]:checkpoint | gzip > [project]-checkpoint.tar.gz
```

## FAQ

**Q: Do I need SSH keys for GitHub?**
A: No! This setup uses HTTPS + GitHub CLI tokens, which is simpler and more secure for dev containers.

**Q: How do I authenticate with GitHub?**
A: Run `gh auth login` on your host machine once. VS Code automatically forwards the credentials.

**Q: Do I need to commit and push every change?**
A: No, work freely in the container. Only push when you want to save work permanently or share with others.

**Q: Can I edit files on my host and see changes in container?**
A: No, this is a fully isolated environment. Edit inside the container or push/pull via git.

**Q: What happens if my computer restarts?**
A: Container stops. Run `docker-compose up -d` to restart it. Uncommitted work is preserved.

**Q: Can I use this on [different platform]?**
A: [Yes/No, platform-specific instructions]

**Q: How do I update [language] version?**
A: Change `[LANGUAGE]_VERSION` in `docker-compose.yml` and rebuild: `./build.sh`

**Q: Is the [X]GB image size too large?**
A: [Explanation of what contributes to size and trade-offs]

**Q: Can I use Claude Code for [language] development?**
A: Yes! Claude Code is pre-installed and the VS Code extension is automatically loaded. Use it for AI-assisted coding, debugging, and algorithm implementation.

**Q: Why is Node.js installed in a [language] container?**
A: Claude Code is distributed as an npm package, so Node.js is required. [Optional: This also enables JavaScript/TypeScript development if needed.]

**Q: Does Claude Code work offline?**
A: No, Claude Code requires internet access to call the Anthropic API. However, [language-analyzer] and other local tools work offline.

---

For issues or questions, see main project README or Claude Code documentation.
