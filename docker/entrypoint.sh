#!/bin/bash
set -e

# Run as root for initial setup
if [ "$(id -u)" = "0" ]; then
  # Ensure critical directories are owned by node user
  chown -R node:node /workspace 2>/dev/null || true
  chown -R node:node /home/node/.claude 2>/dev/null || true
  chown -R node:node /home/node/.cargo 2>/dev/null || true
  chown -R node:node /home/node/.config 2>/dev/null || true

  # Create Claude debug directory with proper permissions
  mkdir -p /home/node/.claude/debug 2>/dev/null || true
  chown -R node:node /home/node/.claude 2>/dev/null || true

  # Check if repository already exists (container restart)
  if [ ! -d "/workspace/rivrs-linalg/.git" ]; then
    echo "Cloning repository from GitHub..."
    cd /workspace

    # Authenticate gh CLI for node user if GH_TOKEN is available
    if [ -n "${GH_TOKEN}" ]; then
      su - node -c "echo '${GH_TOKEN}' | gh auth login --with-token" 2>/dev/null || true
    fi

    # Clone repository as node user
    if ! su - node -c "cd /workspace && gh repo clone jcallaham/rivrs-linalg rivrs-linalg" 2>/tmp/clone_error.log; then
      echo "ERROR: Failed to clone repository"
      cat /tmp/clone_error.log
      echo ""
      echo "GitHub authentication options:"
      echo "  1. On host: gh auth login && gh auth token"
      echo "  2. Run container with: docker run -e GH_TOKEN=\$(gh auth token) ..."
      echo "  3. Mount gh config: -v ~/.config/gh:/home/node/.config/gh:ro"
      exit 1
    fi

    # Configure git
    su - node -c 'cd /workspace/rivrs-linalg && git config --local safe.directory /workspace/rivrs-linalg'
    su - node -c 'cd /workspace/rivrs-linalg && git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"'
    su - node -c 'cd /workspace/rivrs-linalg && git fetch --all --quiet'
    su - node -c 'cd /workspace/rivrs-linalg && git lfs pull'

    # Symlink reference materials
    ln -sf /opt/references /workspace/rivrs-linalg/references
    chown -h node:node /workspace/rivrs-linalg/references

    echo "Repository cloned successfully!"
  else
    echo "Repository already exists, skipping clone."

    # Ensure symlink exists
    if [ ! -e "/workspace/rivrs-linalg/references" ]; then
      ln -sf /opt/references /workspace/rivrs-linalg/references
      chown -h node:node /workspace/rivrs-linalg/references
    fi

    # Re-authenticate gh CLI if GH_TOKEN is available and gh not authenticated
    if [ -n "${GH_TOKEN}" ]; then
      if ! su - node -c "gh auth status" >/dev/null 2>&1; then
        echo "Re-authenticating GitHub CLI..."
        su - node -c "echo '${GH_TOKEN}' | gh auth login --with-token" 2>/dev/null || true
      fi
    fi
  fi

  # Extract SuiteSparse test matrices from reference archive (if available and not already extracted)
  SUITESPARSE_ARCHIVE="/opt/references/ssids/suitesparse.tar.gz"
  SUITESPARSE_DEST="/workspace/rivrs-linalg/sparse/test-data/suitesparse"
  if [ -f "$SUITESPARSE_ARCHIVE" ] && [ ! -f "$SUITESPARSE_DEST/.extracted" ]; then
    echo "Extracting SuiteSparse test matrices..."
    su - node -c "mkdir -p '$SUITESPARSE_DEST'"
    su - node -c "tar xzf '$SUITESPARSE_ARCHIVE' -C '$SUITESPARSE_DEST'" 2>/dev/null && \
      su - node -c "touch '$SUITESPARSE_DEST/.extracted'" && \
      echo "SuiteSparse matrices extracted ($(du -sh "$SUITESPARSE_DEST" | awk '{print $1}'))" || \
      echo "WARNING: Failed to extract SuiteSparse matrices (non-fatal)"
  fi

  # Ensure ownership is correct for volumes (may have been created with wrong permissions)
  chown -R node:node /home/node/.claude 2>/dev/null || true
  chown -R node:node /home/node/.cargo/registry 2>/dev/null || true
  chown -R node:node /home/node/.cache/sccache 2>/dev/null || true

  # Switch to node user and execute command in workspace
  cd /workspace/rivrs-linalg
  exec gosu node "$@"
fi

# If already node user, just execute
cd /workspace/rivrs-linalg
exec "$@"
