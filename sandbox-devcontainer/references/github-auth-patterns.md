# GitHub Authentication Patterns for Devcontainers

Comprehensive guide to setting up GitHub authentication in sandboxed devcontainers using HTTPS and GitHub CLI.

## Why HTTPS Over SSH

### Advantages of HTTPS + GitHub CLI

1. **Automatic Credential Forwarding**: VS Code Dev Containers automatically forward GitHub CLI credentials
2. **No Key Management**: No need to generate, copy, or manage SSH keys in containers
3. **Simpler Setup**: One-time `gh auth login` on host, works everywhere
4. **Read-Only Mounting**: Credentials are mounted read-only, never embedded in images
5. **Better Security**: Tokens can be scoped and revoked easily
6. **Native Integration**: GitHub CLI is designed for this workflow

### SSH Limitations in Containers

1. **Manual Key Generation**: Must generate keys inside container or copy from host
2. **Image Contamination**: Risk of accidentally embedding private keys in images
3. **Permission Issues**: SSH keys need exact permissions (0600) which can break in Docker
4. **No Auto-Forwarding**: VS Code doesn't automatically forward SSH agents
5. **Rotation Complexity**: Rotating keys requires updating multiple containers

## Critical Discovery: macOS Keychain Issue

**Problem**: On macOS, `gh auth login` stores tokens in the macOS Keychain by default, not in `~/.config/gh/hosts.yml`. Containers cannot access the Keychain, so authentication fails even though credentials are mounted.

**Solution**: Use `GH_TOKEN` environment variable to pass the token into the container.

## Complete Setup Guide

### Step 1: Install GitHub CLI on Host (One-Time)

**macOS**:
```bash
brew install gh
```

**Linux (Debian/Ubuntu)**:
```bash
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | \
  sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg

echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] \
  https://cli.github.com/packages stable main" | \
  sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null

sudo apt update && sudo apt install gh
```

**Windows**:
```powershell
winget install --id GitHub.cli
```

### Step 2: Authenticate GitHub CLI (One-Time)

```bash
# Start authentication flow
gh auth login

# Follow prompts:
# 1. Select "GitHub.com"
# 2. Select "HTTPS" protocol
# 3. Select "Login with a web browser" (recommended)
# 4. Copy one-time code shown
# 5. Press Enter to open browser
# 6. Paste code and authorize
```

**Verify authentication**:
```bash
gh auth status

# Expected output:
# github.com
#   ✓ Logged in to github.com as username (/Users/you/.config/gh/hosts.yml)
#   ✓ Git operations for github.com configured to use https protocol.
#   ✓ Token: *******************
```

**Check credential location**:
```bash
ls -la ~/.config/gh/

# Should see:
# config.yml
# hosts.yml
```

### Step 3: Configure Dockerfile

**CRITICAL**: Install GitHub CLI in the Dockerfile itself, not just via devcontainer feature:

```dockerfile
# Stage 1: Base builder
FROM debian:bookworm-slim AS base-builder

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl git build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install GitHub CLI (REQUIRED - not optional)
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*

# Configure git to use GitHub CLI as credential helper
RUN git config --global user.name "${GIT_USER_NAME}" && \
    git config --global user.email "${GIT_USER_EMAIL}" && \
    git config --global init.defaultBranch main && \
    git config --global credential.helper "" && \
    git config --global credential.helper "!gh auth git-credential"

# ... rest of Dockerfile
```

**Key points**:
- Install `gh` CLI in the Dockerfile (don't rely only on devcontainer feature)
- Configure git credential helper explicitly: `git config --global credential.helper "!gh auth git-credential"`
- This makes `git` commands use `gh` for authentication automatically

### Step 4: Pass Token via Environment Variable (macOS Keychain Workaround)

**For standalone Docker scripts** (docker/run.sh):

```bash
#!/bin/bash
set -euo pipefail

# Get GitHub token from host (works with macOS Keychain)
GH_TOKEN=$(gh auth token 2>/dev/null || echo "")

# Run the container
docker run -it \
    --name project-dev \
    -v project-workspace:/workspace/project \
    -e GH_TOKEN="${GH_TOKEN}" \
    -e GITHUB_TOKEN="${GH_TOKEN}" \
    # ... other flags
    project-dev:latest
```

**For docker-compose.yml**:

```yaml
environment:
  - GH_TOKEN=${GH_TOKEN:-}
  - GITHUB_TOKEN=${GITHUB_TOKEN:-}
  # ... other vars
```

**For devcontainer.json**:

```json
{
  "remoteEnv": {
    "GH_TOKEN": "${localEnv:GH_TOKEN}",
    "GITHUB_TOKEN": "${localEnv:GITHUB_TOKEN}"
  }
}
```

Then on your host, export the token:
```bash
# Add to ~/.zshrc or ~/.bash_profile
export GH_TOKEN=$(gh auth token)
```

**Why this works**:
- `gh auth token` extracts token from macOS Keychain
- Token is passed as environment variable (never stored in files)
- `gh` CLI inside container reads from `GH_TOKEN` env var
- Works whether credentials are in Keychain or file-based

### Step 5: Configure devcontainer.json

Add GitHub CLI feature (optional, since we install in Dockerfile):

```json
{
	"name": "Project Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "project-dev",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached"
	],

	// ... rest of config
}
```

**Key components**:
- **feature**: Installs `gh` CLI in container
- **mount**: Forwards host credentials to container (read-only)
- **path**: Uses environment variable `${localEnv:HOME}` for cross-platform compatibility

### Step 5: Configure Git in Dockerfile

Set up git to use HTTPS protocol:

```dockerfile
# Configure git for HTTPS workflow
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main
```

**Note**: These are just display values for commits. Authentication uses GitHub CLI tokens.

### Step 6: Test Inside Container

After opening project in container:

```bash
# Verify GitHub CLI is authenticated
gh auth status

# Expected: ✓ Logged in to github.com as username

# Test git clone via HTTPS
git clone https://github.com/username/repo.git

# Should work without password prompt

# Test push
cd existing-repo
git add .
git commit -m "test commit"
git push origin branch-name

# Should work without password prompt
```

## How It Works

### Credential Flow

1. **Host Authentication**: `gh auth login` stores token in `~/.config/gh/hosts.yml`
2. **Container Mount**: VS Code mounts `~/.config/gh` into container at `/root/.config/gh`
3. **GitHub CLI Integration**: `gh` CLI inside container reads mounted credentials
4. **Git Configuration**: Git uses `gh` as credential helper automatically
5. **HTTPS Operations**: All `git clone/push/pull` use HTTPS with token authentication

### Token Storage

**On Host**:
```
~/.config/gh/
├── config.yml          # General configuration
└── hosts.yml           # Credentials (contains OAuth token)
```

**In Container**:
```
/root/.config/gh/       # Mounted from host (read-only)
├── config.yml
└── hosts.yml
```

**Security**: Token is never copied into image layers, only mounted at runtime.

### Git Credential Helper

GitHub CLI automatically configures git to use `gh` as credential helper:

```bash
# Inside container, check git config
git config --global credential.helper

# Output:
# /usr/bin/gh auth git-credential
```

When git needs credentials for `https://github.com`:
1. Git calls `gh auth git-credential`
2. GitHub CLI reads token from `/root/.config/gh/hosts.yml`
3. Returns token to git
4. Git uses token for HTTPS authentication

## Advanced Patterns

### Custom Git Configuration

If you need custom git settings:

```dockerfile
# In Dockerfile
RUN git config --global core.editor vim && \
	git config --global pull.rebase false && \
	git config --global push.default current
```

### Enterprise GitHub

For GitHub Enterprise Server:

```bash
# On host, authenticate to enterprise
gh auth login --hostname github.enterprise.com

# In devcontainer.json, mount works the same
# GitHub CLI supports multiple hosts in hosts.yml
```

### Token Scopes

Control what the token can do:

```bash
# Re-authenticate with specific scopes
gh auth login --scopes "repo,read:org,workflow"

# Check current scopes
gh auth status
```

### Multiple GitHub Accounts

GitHub CLI supports multiple accounts:

```bash
# Authenticate additional account
gh auth login --hostname github.com

# Switch between accounts
gh auth switch

# Or specify per-command
gh repo list --user other-username
```

### Token Refresh

Tokens are long-lived (1 year default), but can be refreshed:

```bash
# Check token expiration
gh auth status

# Refresh if needed
gh auth refresh
```

## Private Repository Pattern (Runtime Clone)

For **private repositories**, clone at container startup (not during build):

**Entrypoint script in Dockerfile**:

```dockerfile
# Create entrypoint script that clones repository on first run
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'set -e' >> /usr/local/bin/entrypoint.sh && \
    echo '' >> /usr/local/bin/entrypoint.sh && \
    echo '# Check if repository already exists (container restart)' >> /usr/local/bin/entrypoint.sh && \
    echo 'if [ ! -d "/workspace/project/.git" ]; then' >> /usr/local/bin/entrypoint.sh && \
    echo '  echo "Cloning repository from GitHub..."' >> /usr/local/bin/entrypoint.sh && \
    echo '  mkdir -p /workspace' >> /usr/local/bin/entrypoint.sh && \
    echo '  cd /workspace' >> /usr/local/bin/entrypoint.sh && \
    echo '  ' >> /usr/local/bin/entrypoint.sh && \
    echo '  # Use gh CLI to clone (uses GH_TOKEN env var)' >> /usr/local/bin/entrypoint.sh && \
    echo '  if ! gh repo clone username/repo project; then' >> /usr/local/bin/entrypoint.sh && \
    echo '    echo "ERROR: Failed to clone. Ensure gh auth token is set."' >> /usr/local/bin/entrypoint.sh && \
    echo '    exit 1' >> /usr/local/bin/entrypoint.sh && \
    echo '  fi' >> /usr/local/bin/entrypoint.sh && \
    echo '  cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  ' >> /usr/local/bin/entrypoint.sh && \
    echo '  # CRITICAL: Fix fetch refspec (gh repo clone sometimes misses this)' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config --local safe.directory /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"' >> /usr/local/bin/entrypoint.sh && \
    echo '  git fetch --all --quiet' >> /usr/local/bin/entrypoint.sh && \
    echo '  echo "Repository cloned successfully!"' >> /usr/local/bin/entrypoint.sh && \
    echo 'else' >> /usr/local/bin/entrypoint.sh && \
    echo '  echo "Repository already exists, skipping clone."' >> /usr/local/bin/entrypoint.sh && \
    echo 'fi' >> /usr/local/bin/entrypoint.sh && \
    echo '' >> /usr/local/bin/entrypoint.sh && \
    echo 'cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo 'exec "$@"' >> /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
```

**Why entrypoint approach**:
- Build-time clone fails (no credentials available during `docker build`)
- Entrypoint runs after container starts with `GH_TOKEN` available
- Uses `gh repo clone` which respects `GH_TOKEN` env var
- **Critical fix**: Sets fetch refspec (gh clone sometimes omits this)

**Fetch refspec issue**:
Without the fetch refspec configuration, `git fetch` won't download branches:
```bash
# Symptom: git branch -a shows nothing except current branch
# Fix: git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
```

## Troubleshooting

### macOS Keychain: Credentials mounted but not working

**Symptom**: `gh auth status` shows "token is invalid" even though host shows authenticated

**Cause**: macOS stores token in Keychain, not in `~/.config/gh/hosts.yml` file

**Solution**: Use `GH_TOKEN` environment variable approach (see Step 4 above)

### Git operations ask for password despite gh being authenticated

**Symptom**: `gh auth status` works but `git push` asks for username/password

**Cause**: Git doesn't know to use GitHub CLI for credentials

**Solution**: Configure git credential helper in Dockerfile:
```dockerfile
RUN git config --global credential.helper "!gh auth git-credential"
```

### Can't checkout branches after cloning

**Symptom**: `git branch -a` shows only current branch, `git checkout other-branch` fails

**Cause**: Fetch refspec not configured (gh repo clone bug)

**Solution**: Add to entrypoint after cloning:
```bash
git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all --quiet
```

### "Permission denied" when pushing

**Symptom**: `git push` fails with permission error

**Causes**:
1. Not authenticated on host
2. Mount not configured
3. Token expired
4. Insufficient token scopes

**Solutions**:
```bash
# On host, verify authentication
gh auth status

# If not authenticated
gh auth login

# If authenticated but still failing, refresh token
gh auth refresh --scopes "repo,workflow"

# Restart container
# In VS Code: F1 → "Dev Containers: Rebuild Container"
```

### GitHub CLI not found in container

**Symptom**: `gh: command not found`

**Cause**: GitHub CLI feature not configured

**Solution**:
Add to devcontainer.json:
```json
"features": {
	"ghcr.io/devcontainers/features/github-cli:1": {
		"version": "latest"
	}
}
```

Rebuild container.

### Credentials not accessible in container

**Symptom**: `gh auth status` shows not logged in

**Causes**:
1. Mount not configured correctly
2. Path doesn't exist on host
3. Docker file sharing not enabled

**Solutions**:
```bash
# On host, verify credentials exist
ls -la ~/.config/gh/hosts.yml

# Check mount in devcontainer.json
"mounts": [
	"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached"
]

# On macOS, ensure Docker Desktop has access
# Preferences → Resources → File Sharing
# Add: /Users/yourname/.config/gh

# Rebuild container
```

### Clone works but push fails

**Symptom**: `git clone` succeeds but `git push` requires password

**Cause**: Repository uses SSH URL instead of HTTPS

**Solution**:
```bash
# Check remote URL
git remote -v

# If using git@github.com: format, change to HTTPS
git remote set-url origin https://github.com/username/repo.git

# Verify
git remote -v

# Now push should work
git push
```

### Works on host but not in container

**Symptom**: `gh` works on host, fails in container

**Causes**:
1. Container not rebuilt after adding feature
2. Mount path mismatch
3. User mismatch (root vs non-root)

**Solutions**:
```bash
# Rebuild container completely
# F1 → "Dev Containers: Rebuild Container"

# Check mount inside container
ls -la /root/.config/gh/

# If using non-root user, adjust mount target
"mounts": [
	"source=${localEnv:HOME}/.config/gh,target=/home/developer/.config/gh,type=bind,consistency=cached"
]
```

### Docker Compose vs VS Code Dev Containers

**Issue**: Manual docker-compose doesn't auto-mount credentials

**Solution**:

For docker-compose CLI usage (not VS Code):

1. Uncomment manual mount in docker-compose.yml:
```yaml
volumes:
  # - ${HOME}/.config/gh:/root/.config/gh:ro
```

2. Rebuild:
```bash
docker-compose down
docker-compose up -d
```

**Note**: VS Code Dev Containers handle this automatically via devcontainer.json mounts.

## Best Practices

### DO:
- Use `gh auth login` on host before container setup
- Mount credentials read-only (`:ro` flag for manual mounts)
- Use HTTPS URLs for git repositories
- Test `gh auth status` after container creation
- Keep `gh` CLI updated on host (`gh upgrade`)
- Use scoped tokens (minimum necessary permissions)

### DON'T:
- Copy `~/.config/gh` into Dockerfile (embeds secrets in image)
- Use SSH keys for GitHub in containers
- Commit `.config/gh/hosts.yml` to version control
- Share credentials between users (each dev authenticates individually)
- Use personal access tokens directly (let `gh` manage tokens)

## Docker Compose Configuration

For completeness, docker-compose.yml configuration:

```yaml
version: '3.8'

services:
  project-dev:
    build:
      context: ..
      dockerfile: docker/Dockerfile

    volumes:
      # VS Code Dev Containers handles this via devcontainer.json mounts
      # For docker-compose CLI usage, uncomment:
      # - ${HOME}/.config/gh:/root/.config/gh:ro

    # ... rest of config
```

**Note**: The mount is commented out because VS Code Dev Containers manages it. Uncomment only for standalone docker-compose usage.

## Security Considerations

### Token Storage

**Safe**:
- Tokens stored in `~/.config/gh/hosts.yml` on host
- Mounted read-only into container
- Never copied into image layers
- Never committed to git

**Unsafe**:
- Copying `hosts.yml` into Dockerfile
- Setting token as environment variable in Dockerfile
- Embedding token in git configuration

### Token Lifecycle

1. **Generation**: `gh auth login` creates OAuth token
2. **Storage**: Encrypted in `hosts.yml`
3. **Usage**: Mounted into containers at runtime
4. **Rotation**: `gh auth refresh` or re-login
5. **Revocation**: `gh auth logout` or GitHub settings

### Sharing Containers

When sharing Docker images:
- Credentials are NOT included (only in mounts)
- Other users authenticate with their own `gh auth login`
- Each developer uses their own credentials
- No risk of credential leakage in images

## Migration from SSH

If currently using SSH, migrate to HTTPS:

### Step 1: Update Remote URLs

```bash
# Show current remotes
git remote -v

# Change from SSH to HTTPS
git remote set-url origin https://github.com/username/repo.git

# Verify
git remote -v
```

### Step 2: Remove SSH Configuration

```bash
# Remove SSH config from git
git config --global --unset core.sshCommand

# Remove GitHub SSH host config
# Edit ~/.ssh/config, remove github.com sections
```

### Step 3: Authenticate with GitHub CLI

```bash
gh auth login
# Select HTTPS protocol
```

### Step 4: Test

```bash
git pull
git push

# Should work without SSH keys
```

## Additional Patterns

### Symlinked Reference Directories

When symlinking reference materials into the repository, use .gitignore patterns **without trailing slash**:

**Correct** (.gitignore):
```
faer-rs
lapack
slicot
SLICOT-Reference
```

**Incorrect**:
```
faer-rs/
lapack/
slicot/
SLICOT-Reference/
```

**Why**: Git treats symlinks differently than directories. Patterns with trailing `/` match directories only, not symlinks.

### Workspace Volume (Not Mount)

For sandboxed environments, use Docker **volume** (not mount) for isolation:

**Correct** (docker-compose.yml):
```yaml
volumes:
  - workspace:/workspace/project  # Docker volume (isolated)
```

**Incorrect** (breaks sandbox):
```yaml
volumes:
  - ..:/workspace/project  # Mount (exposes host files to container)
```

**run.sh**:
```bash
docker run \
    -v project-workspace:/workspace/project \  # Volume
    # NOT: -v $(pwd):/workspace/project      # Mount
```

**Why**: Volumes provide complete isolation. Mounts give AI agents access to host filesystem.

## Summary

Complete GitHub authentication pattern for sandboxed devcontainers:

1. **Install `gh` CLI in Dockerfile** (not just devcontainer feature):
   ```dockerfile
   RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | ... && apt install gh
   ```

2. **Configure git credential helper** (critical):
   ```dockerfile
   RUN git config --global credential.helper "!gh auth git-credential"
   ```

3. **Pass token via environment variable** (macOS Keychain workaround):
   ```bash
   # run.sh
   GH_TOKEN=$(gh auth token)
   docker run -e GH_TOKEN="${GH_TOKEN}" ...
   ```

4. **Clone at runtime for private repos** (entrypoint script):
   ```dockerfile
   RUN echo 'gh repo clone user/repo project' > /entrypoint.sh
   ENTRYPOINT ["/entrypoint.sh"]
   ```

5. **Fix fetch refspec after clone** (critical):
   ```bash
   git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
   git fetch --all
   ```

6. **Use Docker volumes, not mounts** (isolation):
   ```yaml
   volumes:
     - workspace:/workspace/project  # Volume, NOT ../:/workspace/project
   ```

This complete pattern is:
- **Secure**: Tokens via env var, never in images, workspace isolated from host
- **Works on macOS**: Handles Keychain vs file-based credentials
- **Handles private repos**: Runtime clone with authentication
- **Fetches all branches**: Fixes gh clone fetch refspec bug
- **Production-ready**: Battle-tested on CSRRS project
