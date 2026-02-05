# Stage 1: Base builder with Rust toolchain
FROM debian:bookworm-slim AS base-builder

ARG RUST_VERSION=1.93.0
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
	ca-certificates \
	git \
	build-essential \
	gcc \
	g++ \
	gfortran \
	cmake \
	openssh-client \
	pkg-config \
	libssl-dev \
	&& rm -rf /var/lib/apt/lists/*

# Install GitHub CLI
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
	&& chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
	&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
	&& apt-get update \
	&& apt-get install -y gh \
	&& rm -rf /var/lib/apt/lists/*

# Install Rust (auto-detects ARM64 architecture)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- \
	-y \
	--default-toolchain ${RUST_VERSION} \
	--profile default \
	--component rustfmt,clippy,rust-analyzer

ENV PATH="/root/.cargo/bin:${PATH}"

# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install additional development tools (full suite)
RUN apt-get update && apt-get install -y \
	fzf \
	zsh \
	man-db \
	vim \
	nano \
	unzip \
	gnupg2 \
	dnsutils \
	jq \
	&& rm -rf /var/lib/apt/lists/*

# Install additional Cargo tools (without sccache wrapper initially)
RUN cargo install sccache cargo-watch cargo-expand

# Now set sccache as the wrapper
ENV RUSTC_WRAPPER=sccache

# Configure git (HTTPS-based workflow)
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main && \
	git config --global credential.helper "" && \
	git config --global credential.helper "!gh auth git-credential"

# Stage 2: Reference materials (separate for caching)
FROM scratch AS references
# Copy the entire references directory
COPY references /references

# Stage 3: Development environment
FROM base-builder AS development

# Copy reference materials to /opt/references (permanent location)
COPY --from=references /references/ /opt/references/

# Stage 4: Python environment

# Install UV
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"
ENV UV_LINK_MODE=copy

# Install GitHub spec-kit
RUN uv tool install specify-cli --from git+https://github.com/github/spec-kit.git

# Create entrypoint script that clones repository on first run
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
	echo 'set -e' >> /usr/local/bin/entrypoint.sh && \
	echo '' >> /usr/local/bin/entrypoint.sh && \
	echo '# Check if repository already exists (container restart)' >> /usr/local/bin/entrypoint.sh && \
	echo 'if [ ! -d "/workspace/rivrs-linalg/.git" ]; then' >> /usr/local/bin/entrypoint.sh && \
	echo '  echo "Cloning repository from GitHub..."' >> /usr/local/bin/entrypoint.sh && \
	echo '  mkdir -p /workspace' >> /usr/local/bin/entrypoint.sh && \
	echo '  cd /workspace' >> /usr/local/bin/entrypoint.sh && \
	echo '  ' >> /usr/local/bin/entrypoint.sh && \
	echo '  # Use gh CLI to clone (uses mounted credentials)' >> /usr/local/bin/entrypoint.sh && \
	echo '  if ! gh repo clone jcallaham/rivrs-linalg rivrs-linalg 2>/tmp/clone_error.log; then' >> /usr/local/bin/entrypoint.sh && \
	echo '    echo "ERROR: Failed to clone repository. Check GitHub authentication:"' >> /usr/local/bin/entrypoint.sh && \
	echo '    cat /tmp/clone_error.log' >> /usr/local/bin/entrypoint.sh && \
	echo '    echo ""' >> /usr/local/bin/entrypoint.sh && \
	echo '    echo "Please ensure gh CLI is authenticated on your host machine:"' >> /usr/local/bin/entrypoint.sh && \
	echo '    echo "  gh auth login"' >> /usr/local/bin/entrypoint.sh && \
	echo '    exit 1' >> /usr/local/bin/entrypoint.sh && \
	echo '  fi' >> /usr/local/bin/entrypoint.sh && \
	echo '  ' >> /usr/local/bin/entrypoint.sh && \
	echo '  cd /workspace/rivrs-linalg' >> /usr/local/bin/entrypoint.sh && \
	echo '  git config --local safe.directory /workspace/rivrs-linalg' >> /usr/local/bin/entrypoint.sh && \
	echo '  ' >> /usr/local/bin/entrypoint.sh && \
	echo '  # Fix fetch refspec (gh repo clone sometimes misses this)' >> /usr/local/bin/entrypoint.sh && \
	echo '  git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"' >> /usr/local/bin/entrypoint.sh && \
	echo '  git fetch --all --quiet' >> /usr/local/bin/entrypoint.sh && \
	echo '' >> /usr/local/bin/entrypoint.sh && \
	echo '  # Symlink reference materials into repository' >> /usr/local/bin/entrypoint.sh && \
	echo '  ln -s /opt/references/faer-rs /workspace/rivrs-linalg/faer-rs' >> /usr/local/bin/entrypoint.sh && \
	echo '  ln -s /opt/references/lapack /workspace/rivrs-linalg/lapack' >> /usr/local/bin/entrypoint.sh && \
	echo '  ln -s /opt/references/SLICOT-Reference /workspace/rivrs-linalg/SLICOT-Reference' >> /usr/local/bin/entrypoint.sh && \
	echo '  ln -s /opt/references/slicot /workspace/rivrs-linalg/slicot' >> /usr/local/bin/entrypoint.sh && \
	echo '  ' >> /usr/local/bin/entrypoint.sh && \
	echo '  echo "Repository cloned successfully!"' >> /usr/local/bin/entrypoint.sh && \
	echo 'else' >> /usr/local/bin/entrypoint.sh && \
	echo '  echo "Repository already exists, skipping clone."' >> /usr/local/bin/entrypoint.sh && \
	echo 'fi' >> /usr/local/bin/entrypoint.sh && \
	echo '' >> /usr/local/bin/entrypoint.sh && \
	echo 'cd /workspace/rivrs-linalg' >> /usr/local/bin/entrypoint.sh && \
	echo 'exec "$@"' >> /usr/local/bin/entrypoint.sh && \
	chmod +x /usr/local/bin/entrypoint.sh

WORKDIR /workspace/rivrs-linalg

# Configure Cargo for performance
RUN mkdir -p /root/.cargo && \
	echo '[build]' > /root/.cargo/config.toml && \
	echo 'jobs = 8' >> /root/.cargo/config.toml && \
	echo 'incremental = true' >> /root/.cargo/config.toml && \
	echo '' >> /root/.cargo/config.toml && \
	echo '[profile.dev]' >> /root/.cargo/config.toml && \
	echo 'opt-level = 3' >> /root/.cargo/config.toml

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "======================================"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "rivrs-linalg Monorepo Development Environment (ARM64)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "======================================"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust:        $(rustc --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Cargo:       $(cargo --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Platform:    $(uname -m)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
