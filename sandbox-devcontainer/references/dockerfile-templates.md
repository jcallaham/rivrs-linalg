# Dockerfile Templates by Language

This reference provides complete Dockerfile templates for different language ecosystems, all including Claude Code integration.

## Python Data Science Template

```dockerfile
# Stage 1: Base builder with Python and Claude Code
FROM debian:bookworm-slim AS base-builder

ARG PYTHON_VERSION=3.11
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
	git \
	build-essential \
	python3 \
	python3-pip \
	python3-venv \
	libssl-dev \
	libffi-dev \
	&& rm -rf /var/lib/apt/lists/*

# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY datasets/ /references/datasets/
# COPY docs/ /references/docs/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Copy project files
COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt

# Create project structure
RUN mkdir -p src tests notebooks .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Python Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Python:      $(python3 --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Pip:         $(pip3 --version | cut -d\" \" -f2)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Rust Systems Programming Template

```dockerfile
# Stage 1: Base builder with Rust and Claude Code
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
	pkg-config \
	libssl-dev \
	&& rm -rf /var/lib/apt/lists/*

# Install Rust
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

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Install Cargo tools
RUN cargo install sccache cargo-watch cargo-expand

# Configure sccache
ENV RUSTC_WRAPPER=sccache

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY external-lib/ /references/external-lib/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Copy project files
COPY Cargo.toml ./
COPY .rustfmt.toml ./

# Create project structure
RUN mkdir -p src tests examples benches .claude

# Configure Cargo for performance
RUN mkdir -p /root/.cargo && \
	echo '[build]' > /root/.cargo/config.toml && \
	echo 'jobs = 8' >> /root/.cargo/config.toml && \
	echo 'incremental = true' >> /root/.cargo/config.toml

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust:        $(rustc --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Cargo:       $(cargo --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Node.js/TypeScript Template

```dockerfile
# Stage 1: Base builder with Node.js (serves both project and Claude Code)
FROM debian:bookworm-slim AS base-builder

ARG NODE_VERSION=20
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
	git \
	build-essential \
	&& rm -rf /var/lib/apt/lists/*

# Install Node.js (same for project and Claude Code)
RUN curl -fsSL https://deb.nodesource.com/setup_${NODE_VERSION}.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY libs/ /references/libs/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Copy package files
COPY package.json package-lock.json ./
RUN npm ci

# Create project structure
RUN mkdir -p src tests dist .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "npm:         $(npm --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Go Template

```dockerfile
# Stage 1: Base builder with Go and Claude Code
FROM debian:bookworm-slim AS base-builder

ARG GO_VERSION=1.22.0
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
	git \
	build-essential \
	&& rm -rf /var/lib/apt/lists/*

# Install Go
RUN curl -LO https://go.dev/dl/go${GO_VERSION}.linux-arm64.tar.gz && \
	tar -C /usr/local -xzf go${GO_VERSION}.linux-arm64.tar.gz && \
	rm go${GO_VERSION}.linux-arm64.tar.gz

ENV PATH="/usr/local/go/bin:/root/go/bin:${PATH}"
ENV GOPATH="/root/go"

# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Install Go tools
RUN go install golang.org/x/tools/gopls@latest && \
	go install github.com/go-delve/delve/cmd/dlv@latest && \
	go install honnef.co/go/tools/cmd/staticcheck@latest

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY vendor/ /references/vendor/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Copy Go module files
COPY go.mod go.sum ./
RUN go mod download

# Create project structure
RUN mkdir -p cmd pkg internal .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Go Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Go:          $(go version | cut -d\" \" -f3)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Java/Maven Template

```dockerfile
# Stage 1: Base builder with Java and Claude Code
FROM debian:bookworm-slim AS base-builder

ARG JAVA_VERSION=17
ARG GIT_USER_NAME="Claude Code"
ARG GIT_USER_EMAIL="claude@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
	curl \
	git \
	build-essential \
	openjdk-${JAVA_VERSION}-jdk \
	maven \
	&& rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME="/usr/lib/jvm/java-${JAVA_VERSION}-openjdk-arm64"
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY lib/ /references/lib/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Copy Maven files
COPY pom.xml ./
RUN mvn dependency:go-offline

# Create project structure
RUN mkdir -p src/main/java src/test/java .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Java Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Java:        $(java -version 2>&1 | head -n1 | cut -d'\''\"'\'' -f2)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Maven:       $(mvn -version | head -n1 | cut -d\" \" -f3)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Multi-Language Polyglot Template

For projects using multiple languages (e.g., Python backend + React frontend, Rust library + Python bindings):

```dockerfile
# Stage 1: Base builder with multiple runtimes
FROM debian:bookworm-slim AS base-builder

ARG PYTHON_VERSION=3.11
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
	cmake \
	pkg-config \
	libssl-dev \
	python3 \
	python3-pip \
	python3-venv \
	&& rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- \
	-y \
	--default-toolchain ${RUST_VERSION} \
	--profile default

ENV PATH="/root/.cargo/bin:${PATH}"

# Install Node.js 20 (for Claude Code AND frontend)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# Install development tools
RUN apt-get update && apt-get install -y \
	fzf zsh man-db vim nano unzip gnupg2 dnsutils jq \
	&& rm -rf /var/lib/apt/lists/*

# Configure git
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main

# Stage 2: Reference materials (if any)
FROM scratch AS references
# COPY references/ /references/

# Stage 3: Development environment
FROM base-builder AS development

WORKDIR /workspace/project

# Copy reference materials
# COPY --from=references /references ./

# Create project structure
RUN mkdir -p backend frontend bindings tests .claude

# Create welcome script
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Polyglot Development Environment"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "====================================="' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Python:      $(python3 --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Rust:        $(rustc --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	echo 'echo ""' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh

CMD ["/bin/bash", "-c", "welcome.sh && exec /bin/bash"]
```

## Common Patterns Across Templates

### Node.js Installation (Required for Claude Code)

All templates include:
```dockerfile
# Install Node.js 20 for Claude Code
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
	apt-get install -y nodejs && \
	npm install -g npm@latest && \
	rm -rf /var/lib/apt/lists/*

# Install Claude Code
ARG CLAUDE_CODE_VERSION=latest
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}
```

### Development Tools Suite

Standard toolkit included in all templates:
```dockerfile
RUN apt-get update && apt-get install -y \
	fzf \       # Fuzzy finder
	zsh \       # Advanced shell
	man-db \    # Manual pages
	vim \       # Editor
	nano \      # Simpler editor
	unzip \     # Archive extraction
	gnupg2 \    # Encryption
	dnsutils \  # Network debugging
	jq \        # JSON processing
	&& rm -rf /var/lib/apt/lists/*
```

### Git Configuration

Essential for all containers:
```dockerfile
RUN git config --global user.name "${GIT_USER_NAME}" && \
	git config --global user.email "${GIT_USER_EMAIL}" && \
	git config --global init.defaultBranch main
```

### Multi-Stage Pattern

All templates use 3-stage builds:
1. **base-builder**: Language runtime + tools + Claude Code
2. **references**: Static reference materials (if any)
3. **development**: Final environment with project files

### Welcome Script

Shows installed versions for verification:
```dockerfile
RUN echo '#!/bin/bash' > /usr/local/bin/welcome.sh && \
	echo 'echo "Language:    $(language --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Node.js:     $(node --version)"' >> /usr/local/bin/welcome.sh && \
	echo 'echo "Claude Code: $(claude --version 2>/dev/null || echo '\''not installed'\'')"' >> /usr/local/bin/welcome.sh && \
	chmod +x /usr/local/bin/welcome.sh
```

## Customization Points

When adapting these templates:

1. **Architecture**: Change `linux/arm64` to `linux/amd64` for Intel/AMD
2. **Versions**: Update ARG values for language/tool versions
3. **System packages**: Add language-specific libraries (e.g., `libpq-dev` for PostgreSQL)
4. **Package managers**: Use project-specific install commands
5. **Cache volumes**: Configure in docker-compose.yml for package caches
6. **Reference stage**: Uncomment and customize COPY commands for static materials
7. **Project structure**: Adjust `mkdir -p` to match project layout
8. **Welcome script**: Customize to show relevant version information

## Platform Considerations

### ARM64 (Apple Silicon)
- Native performance
- Use `linux/arm64` platform
- Fast builds (8-15 minutes)

### AMD64 (Intel/AMD)
- Standard x86_64
- Use `linux/amd64` platform
- Standard build times

### Multi-Platform
```dockerfile
# Auto-detect or specify both
--platform linux/arm64,linux/amd64
```

## Size Optimization

To keep images smaller:

1. **Combine RUN commands**: Reduces layers
   ```dockerfile
   RUN apt-get update && apt-get install -y \
       package1 package2 package3 \
       && rm -rf /var/lib/apt/lists/*
   ```

2. **Clean caches**: Remove after install
   ```dockerfile
   && rm -rf /var/lib/apt/lists/*  # APT
   && rm -rf ~/.cache/pip          # Pip
   && rm -rf ~/.npm                # npm
   ```

3. **Use specific base**: Debian bookworm-slim is minimal
4. **Multi-stage builds**: Only copy what's needed to final stage
5. **Minimize layers**: Combine related operations

Expected sizes:
- Minimal (Python/Node): 2-3GB
- Moderate (Rust/Go): 3-4GB
- Heavy (Java + refs): 4-5GB
- Polyglot + refs: 4-6GB
