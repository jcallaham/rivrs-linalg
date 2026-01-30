# Lessons Learned: Sandbox Devcontainer Implementation

Battle-tested solutions to critical issues discovered during CSRRS implementation (2026-01-30).

## Issue 1: GitHub CLI Not Found

**Problem**: Container starts but `gh: command not found`

**Root Cause**: Relying only on VS Code devcontainer feature to install GitHub CLI. This doesn't work for:
- Standalone Docker usage (docker run, docker-compose)
- Entrypoint scripts that run before features are installed
- Docker builds that need gh during setup

**Solution**: Install GitHub CLI directly in Dockerfile:

```dockerfile
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*
```

**Impact**: Critical - without this, private repo cloning and git authentication fail completely.

---

## Issue 2: macOS Keychain Token Storage

**Problem**: `gh auth status` shows authenticated on host, but inside container shows "token is invalid" even with credentials mounted

**Root Cause**: On macOS, `gh auth login` stores tokens in macOS Keychain by default, not in `~/.config/gh/hosts.yml`. Containers can't access the macOS Keychain, so the mounted `hosts.yml` file doesn't contain the actual token.

**What `hosts.yml` looks like** (missing token):
```yaml
github.com:
    git_protocol: https
    users:
        username:
    user: username
```

**Solution**: Pass token via environment variable:

```bash
# run.sh
GH_TOKEN=$(gh auth token 2>/dev/null || echo "")
docker run -e GH_TOKEN="${GH_TOKEN}" -e GITHUB_TOKEN="${GH_TOKEN}" ...
```

```yaml
# docker-compose.yml
environment:
  - GH_TOKEN=${GH_TOKEN:-}
  - GITHUB_TOKEN=${GITHUB_TOKEN:-}
```

```json
// devcontainer.json
{
  "remoteEnv": {
    "GH_TOKEN": "${localEnv:GH_TOKEN}",
    "GITHUB_TOKEN": "${localEnv:GITHUB_TOKEN}"
  }
}
```

**Why it works**: `gh` CLI automatically reads from `GH_TOKEN` environment variable if present, bypassing file/Keychain storage completely.

**Impact**: High - affects all macOS users. Linux users may be fine with file-based storage.

---

## Issue 3: Git Operations Ask for Password

**Problem**: `gh auth status` shows authenticated, but `git push` asks for username/password

**Root Cause**: Git doesn't automatically know to use GitHub CLI for credential management. Even with `gh` authenticated, git uses its own credential helper (or none).

**Solution**: Explicitly configure git to use GitHub CLI as credential helper:

```dockerfile
RUN git config --global credential.helper "" && \
    git config --global credential.helper "!gh auth git-credential"
```

**Verification**:
```bash
# Inside container
git config --global credential.helper
# Output: !gh auth git-credential

git push
# Should work without password prompt
```

**Impact**: Critical - without this, all git operations fail even with gh authenticated.

---

## Issue 4: Private Repository Clone Fails During Build

**Problem**: `git clone https://github.com/user/private-repo.git` in Dockerfile fails with authentication error

**Root Cause**: Docker build runs before container starts, so:
- No `GH_TOKEN` environment variable available yet
- No mounted credentials available yet
- Can't authenticate to clone private repos

**Solution**: Clone at runtime via entrypoint script:

```dockerfile
# Create entrypoint script
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'if [ ! -d "/workspace/project/.git" ]; then' >> /usr/local/bin/entrypoint.sh && \
    echo '  gh repo clone user/project /workspace/project || exit 1' >> /usr/local/bin/entrypoint.sh && \
    echo '  cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  # Critical fixes after clone' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config --local safe.directory /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"' >> /usr/local/bin/entrypoint.sh && \
    echo '  git fetch --all --quiet' >> /usr/local/bin/entrypoint.sh && \
    echo 'fi' >> /usr/local/bin/entrypoint.sh && \
    echo 'cd /workspace/project && exec "$@"' >> /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
```

**Workflow**:
1. Image builds with entrypoint script (no clone yet)
2. Container starts, entrypoint runs with `GH_TOKEN` available
3. First start: clones repository
4. Subsequent starts: skips clone (already exists)

**Impact**: Critical for private repositories. Public repos can still use build-time clone if desired.

---

## Issue 5: Can't Checkout Other Branches

**Problem**: After cloning, `git branch -a` only shows current branch. `git checkout feature-branch` fails with "pathspec did not match any file(s) known to git"

**Root Cause**: `gh repo clone` sometimes doesn't set the fetch refspec properly. Without it:
```bash
git config --get remote.origin.fetch
# Returns nothing (should return "+refs/heads/*:refs/remotes/origin/*")
```

This means `git fetch` doesn't know to fetch other branches.

**Solution**: Set fetch refspec explicitly after cloning:

```bash
git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all --quiet
```

**Add to entrypoint** right after `gh repo clone`:

```bash
if ! gh repo clone user/repo project; then exit 1; fi
cd project
git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all --quiet  # Now fetches all branches
```

**Verification**:
```bash
git branch -a
# Should show:
#   * main
#   remotes/origin/feature-branch
#   remotes/origin/another-branch

git checkout feature-branch
# Should work now
```

**Impact**: High - users can't work on feature branches without this fix.

---

## Issue 6: Symlinked Reference Directories Show as Untracked

**Problem**: Reference materials are symlinked into repository (e.g., `ln -s /opt/references/docs ./docs`). Despite being in `.gitignore`, `git status` shows them as untracked files.

**Root Cause**: `.gitignore` patterns with trailing slashes (`docs/`) match directories only, not symlinks. Git treats symlinks as files, not directories.

**Incorrect .gitignore**:
```gitignore
docs/
references/
lib/
```

**Correct .gitignore**:
```gitignore
docs
references
lib
```

**Solution**: Use patterns without trailing slash to match both directories and symlinks.

**Verification**:
```bash
git status
# Should not show symlinked directories as untracked
```

**Impact**: Low (cosmetic), but annoying. Clutters git status output.

---

## Issue 7: Host Files Accessible from Container (Breaks Sandbox)

**Problem**: When using `-v $(pwd):/workspace/project`, AI agents and users inside container can access all host files in the project directory, breaking the sandbox model.

**Root Cause**: Mounting host directory into container creates a direct link. Any changes in container affect host, and vice versa. This is a "development container" pattern, not a "sandbox container" pattern.

**Incorrect (breaks sandbox)**:
```bash
docker run -v "$(pwd)":/workspace/project ...
```

```yaml
# docker-compose.yml
volumes:
  - ..:/workspace/project
```

**Correct (isolated sandbox)**:
```bash
docker run -v project-workspace:/workspace/project ...  # Docker volume
```

```yaml
# docker-compose.yml
volumes:
  - workspace:/workspace/project  # Named volume

volumes:
  workspace:
    name: project-workspace
```

**Why**: Docker volumes are isolated storage managed by Docker. Data exists only inside Docker's storage, not on host filesystem. Perfect for sandboxing.

**Workflow with volumes**:
1. Work in isolated container
2. Commit changes to git inside container
3. Push to GitHub from container
4. Pull changes on host if needed
5. Repository exists in two places (host + container volume), synced via git

**Impact**: Critical for security/sandboxing. Mounts completely break the isolation model.

---

## Summary: Complete Working Pattern

Based on all lessons learned, here's the complete pattern:

### 1. Dockerfile

```dockerfile
FROM debian:bookworm-slim AS base-builder

# Install GitHub CLI (CRITICAL)
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*

# Configure git credential helper (CRITICAL)
RUN git config --global user.name "Claude Code" && \
    git config --global user.email "claude@example.com" && \
    git config --global init.defaultBranch main && \
    git config --global credential.helper "" && \
    git config --global credential.helper "!gh auth git-credential"

# Create entrypoint for private repo cloning (CRITICAL for private repos)
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'if [ ! -d "/workspace/project/.git" ]; then' >> /usr/local/bin/entrypoint.sh && \
    echo '  gh repo clone username/project /workspace/project || exit 1' >> /usr/local/bin/entrypoint.sh && \
    echo '  cd /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config --local safe.directory /workspace/project' >> /usr/local/bin/entrypoint.sh && \
    echo '  git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"' >> /usr/local/bin/entrypoint.sh && \
    echo '  git fetch --all --quiet' >> /usr/local/bin/entrypoint.sh && \
    echo 'fi' >> /usr/local/bin/entrypoint.sh && \
    echo 'cd /workspace/project && exec "$@"' >> /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
```

### 2. run.sh

```bash
#!/bin/bash
GH_TOKEN=$(gh auth token 2>/dev/null || echo "")

docker run -it \
    --name project-dev \
    -v project-workspace:/workspace/project \
    -e GH_TOKEN="${GH_TOKEN}" \
    -e GITHUB_TOKEN="${GH_TOKEN}" \
    project-dev:latest
```

### 3. docker-compose.yml

```yaml
version: '3.8'
services:
  project-dev:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    image: project-dev:latest
    volumes:
      - workspace:/workspace/project
    environment:
      - GH_TOKEN=${GH_TOKEN:-}
      - GITHUB_TOKEN=${GITHUB_TOKEN:-}

volumes:
  workspace:
    name: project-workspace
```

### 4. devcontainer.json

```json
{
  "name": "Project Dev",
  "dockerComposeFile": "docker-compose.yml",
  "service": "project-dev",
  "workspaceFolder": "/workspace/project",
  "remoteEnv": {
    "GH_TOKEN": "${localEnv:GH_TOKEN}",
    "GITHUB_TOKEN": "${localEnv:GITHUB_TOKEN}"
  }
}
```

### 5. .gitignore (for symlinks)

```gitignore
# Reference directories (symlinked) - NO trailing slash
references
docs
lib
```

## Testing Checklist

After implementing all fixes:

```bash
# Inside container
gh --version                    # Should show version (Issue 1)
gh auth status                  # Should show authenticated (Issue 2)
git config credential.helper    # Should show "!gh auth git-credential" (Issue 3)
ls /workspace/project/.git      # Should exist (Issue 4)
git branch -a                   # Should show all remote branches (Issue 5)
git status                      # Symlinks should not show as untracked (Issue 6)
ls /host/files                  # Should fail - no host access (Issue 7)
git push                        # Should work without password
git checkout feature-branch     # Should work
```

## Platform-Specific Notes

### macOS
- **Must** use `GH_TOKEN` environment variable (Keychain issue)
- Test with `gh auth token` to ensure token extractable
- For long-term sessions, add `export GH_TOKEN=$(gh auth token)` to `~/.zshrc`

### Linux
- May work with mounted `~/.config/gh` if file-based storage
- `GH_TOKEN` approach still recommended for consistency

### Windows (WSL)
- Same as Linux
- Ensure Docker Desktop has WSL integration enabled
- `GH_TOKEN` approach recommended

## Timeline

All issues discovered and fixed during CSRRS sandbox devcontainer implementation (2026-01-30):

1. ✅ Issue 1: GitHub CLI installation (fixed first)
2. ✅ Issue 2: macOS Keychain (discovered during auth testing)
3. ✅ Issue 3: Git credential helper (discovered when git push failed)
4. ✅ Issue 4: Private repo clone (discovered during build attempts)
5. ✅ Issue 5: Fetch refspec (discovered after successful clone)
6. ✅ Issue 6: Symlink .gitignore (discovered after adding references)
7. ✅ Issue 7: Volume vs mount (architectural decision for sandboxing)

Total time to discover and fix all issues: ~2 hours of iterative troubleshooting.

## References

- **CSRRS Implementation**: See `/docker/Dockerfile`, `/docker/run.sh`, `/.devcontainer/` for working examples
- **GitHub Auth Patterns**: See `references/github-auth-patterns.md` for complete authentication documentation
- **SKILL.md**: See "Critical Implementation Patterns" section for quick reference

---

**Status**: All patterns tested and verified working on macOS ARM64 with private GitHub repository (jcallaham/csrrs).
