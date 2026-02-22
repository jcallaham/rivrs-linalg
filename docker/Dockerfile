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
	linux-perf \
	gcc \
	g++ \
	gfortran \
	cmake \
	openssh-client \
	pkg-config \
	libssl-dev \
	# Additional development tools
	sudo \
	gosu \
	fzf \
	zsh \
	man-db \
	vim \
	nano \
	unzip \
	gnupg2 \
	dnsutils \
	jq \
	procps \
	less \
	&& rm -rf /var/lib/apt/lists/*

# Install GitHub CLI
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
	&& chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
	&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
	&& apt-get update \
	&& apt-get install -y gh \
	&& rm -rf /var/lib/apt/lists/*

# Create node user (uid 1000, same as Anthropic reference)
RUN groupadd --gid 1000 node \
	&& useradd --uid 1000 --gid node --shell /bin/bash --create-home node

# Configure sudo for node user (passwordless)
RUN echo "node ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/node \
	&& chmod 0440 /etc/sudoers.d/node

# Install Rust for node user
USER node
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- \
	-y \
	--default-toolchain ${RUST_VERSION} \
	--profile default \
	--component rustfmt,clippy,rust-analyzer

ENV PATH="/home/node/.cargo/bin:${PATH}"

# Install additional Cargo tools (without sccache wrapper initially)
RUN cargo install sccache cargo-watch cargo-expand

# Now set sccache as the wrapper
ENV RUSTC_WRAPPER=sccache

# Install Claude Code
RUN curl -fsSL https://claude.ai/install.sh | bash

# Configure git (HTTPS-based workflow)
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main && \
	git config --global credential.helper "" && \
	git config --global credential.helper "!gh auth git-credential"

USER root

# Install Git LFS
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
	&& apt-get install -y git-lfs && git lfs install --system

# Stage 2: Reference materials (separate for caching)

FROM scratch AS references
# Copy the entire references directory
COPY references /references

# =============================================================================
# Stage: SPRAL Build (optional - for validation and benchmarking only)
#
# Builds the SPRAL sparse linear algebra library (BSD-3) from source as a
# reference implementation for SSIDS and other sparse solvers.
# See: https://github.com/ralna/spral
#
# To create a leaner development image without SPRAL, remove:
#   - This entire stage (FROM debian:bookworm-slim AS spral-builder ... meson install)
#   - The "SPRAL installation" block in the development stage below
# =============================================================================
FROM debian:bookworm-slim AS spral-builder

# Install build dependencies for SPRAL
RUN apt-get update && apt-get install -y \
	build-essential \
	gfortran \
	meson \
	ninja-build \
	python3 \
	pkg-config \
	git \
	libopenblas-dev \
	libmetis-dev \
	libhwloc-dev \
	&& rm -rf /var/lib/apt/lists/*

# Clone SPRAL source
RUN git clone --depth=1 https://github.com/ralna/spral.git /opt/spral-src

# Build and install SPRAL with examples and tests enabled
RUN meson setup /opt/spral-src/builddir /opt/spral-src \
	-Dexamples=true \
	-Dtests=true \
	-Dlibblas=openblas \
	-Dliblapack=openblas \
	--prefix=/opt/spral \
	&& meson compile -C /opt/spral-src/builddir \
	&& meson install -C /opt/spral-src/builddir

# Stage 3: Development environment
FROM base-builder AS development

# Copy reference materials to /opt/references (permanent location)
COPY --from=references /references/ /opt/references/

# =============================================================================
# SPRAL installation (optional - for validation and benchmarking only)
# Remove this block and the spral-builder stage above for a leaner image
# =============================================================================

# Install SPRAL runtime dependencies
RUN apt-get update && apt-get install -y \
	libopenblas-dev \
	libmetis-dev \
	libhwloc-dev \
	&& rm -rf /var/lib/apt/lists/*

# Copy SPRAL build artifacts from spral-builder stage
COPY --from=spral-builder /opt/spral /opt/spral

# Create workspace and config directories with proper ownership
RUN mkdir -p /workspace /opt/references /home/node/.claude /home/node/.cargo \
	&& chown -R node:node /workspace /opt/references /home/node/.claude /home/node/.cargo

# Stage 4: Python environment
USER node

# Install UV
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/home/node/.local/bin:${PATH}"
ENV UV_LINK_MODE=copy

# Install GitHub spec-kit
RUN uv tool install specify-cli --from git+https://github.com/github/spec-kit.git

USER root

# Configure Cargo for performance (for node user)
RUN echo '[build]' > /home/node/.cargo/config.toml && \
	echo 'jobs = 8' >> /home/node/.cargo/config.toml && \
	echo 'incremental = true' >> /home/node/.cargo/config.toml && \
	echo '' >> /home/node/.cargo/config.toml && \
	echo '[profile.dev]' >> /home/node/.cargo/config.toml && \
	echo 'opt-level = 3' >> /home/node/.cargo/config.toml && \
	chown -R node:node /home/node/.cargo


# Create entrypoint script that clones repository on first run

# Copy entrypoint script
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "======================================"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "rivrs-linalg Monorepo Development Environment ($(uname -m))"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "======================================"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust:        $(rustc --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Cargo:       $(cargo --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Platform:    $(uname -m)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
