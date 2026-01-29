# Stage 1: Base builder with Rust toolchain
FROM debian:bookworm-slim AS base-builder

ARG RUST_VERSION=1.93.0
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
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
	git config --global init.defaultBranch main

# Stage 2: Reference materials (separate for caching)
FROM scratch AS references
COPY faer-rs /references/faer-rs
COPY lapack /references/lapack
COPY SLICOT-Reference /references/SLICOT-Reference
COPY slicot /references/slicot

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/csrrs

# Copy reference materials from stage 2
COPY --from=references /references/faer-rs ./faer-rs
COPY --from=references /references/lapack ./lapack
COPY --from=references /references/SLICOT-Reference ./SLICOT-Reference
COPY --from=references /references/slicot ./slicot

# Copy project metadata files
COPY .gitignore CLAUDE.md README.md ./

# Copy Cargo.toml and rustfmt.toml from faer-rs for reference
COPY faer-rs/Cargo.toml ./faer-Cargo.toml.example
COPY faer-rs/rustfmt.toml ./.rustfmt.toml

# Create project structure
RUN mkdir -p src tests examples benches docs .claude

# Initialize git repository
RUN git init && \
	git config --local safe.directory /workspace/csrrs

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
	echo 'echo "CSRRS Development Environment (ARM64)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "======================================"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust:        $(rustc --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Cargo:       $(cargo --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Platform:    $(uname -m)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Reference materials:"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "  - faer-rs/           (Rust linear algebra)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "  - lapack/            (BSD-licensed reference)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "  - SLICOT-Reference/  (docs and tests only)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "  - slicot/            (docs and tests only)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Clean room reminder: NEVER read slicot/src/*.f files"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
