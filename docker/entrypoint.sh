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
  fi

  # Switch to node user and execute command in workspace
  cd /workspace/rivrs-linalg
  exec gosu node "$@"
fi

# If already node user, just execute
cd /workspace/rivrs-linalg
exec "$@"
