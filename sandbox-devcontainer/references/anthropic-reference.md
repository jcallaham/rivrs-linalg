# Anthropic Official Claude Code Devcontainer Reference

This document captures key learnings from Anthropic's official Claude Code devcontainer implementation available at:
https://github.com/anthropics/claude-code/tree/main/.devcontainer

## Official Implementation Overview

Anthropic provides a reference devcontainer specifically designed for Claude Code development and testing. Understanding this implementation provides valuable insights for building sandboxed containers.

## Key Components from Anthropic's Implementation

### 1. Base Image Choice

**Anthropic uses**: `mcr.microsoft.com/devcontainers/typescript-node:1-20-bookworm`

**Why**:
- Pre-built by Microsoft for devcontainer usage
- Node.js 20 pre-installed (required for Claude Code)
- Debian Bookworm base
- Optimized for VS Code Dev Containers

**For sandbox devcontainers**:
- Language-specific projects should use language-specific base OR
- Build custom from `debian:bookworm-slim` + install language + Node.js
- Anthropic's choice is specific to Node.js-centric development

### 2. Firewall and Network Security

**Anthropic implements**: Sophisticated firewall with IP allowlisting

```bash
# init-firewall.sh (from Anthropic implementation)
# - Blocks all outbound traffic by default
# - Allows specific IPs for:
#   - GitHub (for code access)
#   - npm registry (for package installs)
#   - Anthropic API (for Claude Code functionality)
# - Uses iptables rules
# - Requires NET_ADMIN capability
```

**Key script**: `.devcontainer/init-firewall.sh`

**For sandbox devcontainers**:
- **Optional but recommended for security-critical projects**
- Prevents Claude Code from accessing arbitrary network resources
- Requires `--cap-add=NET_ADMIN` in docker-compose.yml
- Adds complexity to setup and maintenance
- IP addresses need periodic updates (GitHub/npm IPs can change)

**CSRRS decision**: Skipped firewall for simplicity. Can be added later if needed.

### 3. Dockerfile Structure

**Anthropic's Dockerfile** (`.devcontainer/Dockerfile`):

```dockerfile
FROM mcr.microsoft.com/devcontainers/typescript-node:1-20-bookworm

# Install additional tools
RUN apt-get update && apt-get install -y \
    iptables \
    iputils-ping \
    dnsutils \
    && rm -rf /var/lib/apt/lists/*

# Copy and set up firewall script
COPY init-firewall.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/init-firewall.sh
```

**Key characteristics**:
- Minimal additional packages
- Firewall script embedded
- No language runtime installation (already in base)

**For sandbox devcontainers**:
- Multi-stage builds more efficient for complex projects
- Language-specific: Install runtime in Dockerfile
- Embed reference materials in separate stage
- More comprehensive dev tool suite

### 4. devcontainer.json Configuration

**Anthropic's configuration**:

```json
{
  "name": "Claude Code Development",
  "build": {
    "dockerfile": "Dockerfile",
    "context": ".."
  },

  "features": {
    "ghcr.io/devcontainers/features/github-cli:1": {}
  },

  "customizations": {
    "vscode": {
      "extensions": [
        "anthropic.claude-code"
      ]
    }
  },

  "capAdd": ["NET_ADMIN"],  // For firewall

  "postCreateCommand": "/usr/local/bin/init-firewall.sh",

  "mounts": [
    "source=claude-bash-history,target=/home/node/.bash_history,type=volume",
    "source=claude-config,target=/home/node/.claude,type=volume"
  ],

  "remoteUser": "node"
}
```

**Key features**:
- `capAdd: NET_ADMIN` - Required for firewall setup
- `postCreateCommand` - Runs firewall script on container start
- Named volumes for bash history and Claude config
- `remoteUser: node` - Non-root user for security

**For sandbox devcontainers**:
- Skip `capAdd` if no firewall
- Skip `postCreateCommand` firewall script
- Keep Claude config volume mount
- Optionally use root user (simpler permissions) or create custom user

### 5. Volume Strategy

**Anthropic uses**:
- `claude-bash-history` - Persist bash command history
- `claude-config` - Persist Claude Code configuration

**For sandbox devcontainers**:
- Always include Claude config volume
- Add language-specific cache volumes (Cargo cache, pip cache, npm cache, etc.)
- Bash history volume optional but nice to have

### 6. User Configuration

**Anthropic uses**: `remoteUser: node` (non-root)

**Benefits**:
- Better security (limited permissions)
- Closer to production environment
- Follows principle of least privilege

**Trade-offs**:
- More complex permissions
- May need `sudo` for some operations
- Additional configuration for tools

**For sandbox devcontainers**:
- **Development/personal use**: `remoteUser: root` is simpler
- **Team/production**: Create custom non-root user
- CSRRS uses root for simplicity

## Key Differences: Anthropic vs Language-Specific Sandboxes

| Aspect | Anthropic Implementation | Language-Specific Sandbox (e.g., CSRRS) |
|--------|-------------------------|------------------------------------------|
| **Base Image** | Microsoft devcontainer (Node.js) | Custom build from debian:bookworm-slim |
| **Primary Language** | Node.js/TypeScript | Any (Rust, Python, Go, etc.) |
| **Firewall** | Yes, with IP allowlisting | Optional (skipped for simplicity) |
| **User** | `node` (non-root) | `root` (simpler permissions) |
| **Reference Materials** | None (code-only) | Large embedded references (700MB+) |
| **Build Stages** | Single stage | Multi-stage (base → references → dev) |
| **Cache Volumes** | Bash history, Claude config | Language cache + Claude config |
| **Network Security** | Restricted via firewall | Unrestricted (trust environment) |
| **Complexity** | Higher (security-focused) | Lower (development-focused) |

## Lessons Applied to CSRRS

1. **Claude Code Volume**: Adopted persistent volume for Claude config
2. **GitHub CLI**: Used same feature-based installation pattern
3. **Environment Variables**: Set `CLAUDE_CONFIG_DIR` and `NODE_OPTIONS`
4. **VS Code Extension**: Pre-configured in devcontainer.json
5. **Credential Forwarding**: Used same GitHub CLI mount pattern

## Firewall Implementation (Optional Enhancement)

If you want to add Anthropic-style firewall security to your sandbox:

### 1. Create `init-firewall.sh`

```bash
#!/bin/bash
# Based on Anthropic's implementation

# Allow loopback
iptables -A OUTPUT -o lo -j ACCEPT

# Allow established connections
iptables -A OUTPUT -m state --state ESTABLISHED,RELATED -j ACCEPT

# Allow DNS
iptables -A OUTPUT -p udp --dport 53 -j ACCEPT

# Allow GitHub IPs (example - update regularly)
# Get current IPs: curl https://api.github.com/meta | jq .git
iptables -A OUTPUT -d 140.82.112.0/20 -j ACCEPT
iptables -A OUTPUT -d 143.55.64.0/20 -j ACCEPT

# Allow npm registry
iptables -A OUTPUT -d 104.16.0.0/12 -j ACCEPT  # Cloudflare (npm CDN)

# Allow Anthropic API
iptables -A OUTPUT -d 160.79.104.0/23 -j ACCEPT

# Drop everything else
iptables -A OUTPUT -j DROP

echo "Firewall configured"
```

### 2. Update devcontainer.json

```json
{
  "capAdd": ["NET_ADMIN"],
  "postCreateCommand": "/usr/local/bin/init-firewall.sh"
}
```

### 3. Update docker-compose.yml

```yaml
services:
  project-dev:
    cap_add:
      - NET_ADMIN
```

### 4. Copy script in Dockerfile

```dockerfile
COPY docker/init-firewall.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/init-firewall.sh
```

### Firewall Considerations

**Pros**:
- Restricts Claude Code network access to known services
- Prevents accidental data leakage
- Defense in depth security

**Cons**:
- Adds complexity to setup
- Requires NET_ADMIN capability
- IP addresses can change (GitHub, npm, etc.)
- May break legitimate use cases (accessing project-specific APIs)
- Harder to debug network issues

**When to use**:
- Processing sensitive data
- Regulatory compliance requirements
- Untrusted code execution
- Production-like environments

**When to skip**:
- Personal development
- Trusted code only
- Need flexible network access
- Simplicity preferred

## GitHub IP Allowlisting

If implementing firewall, GitHub IPs can be fetched:

```bash
# Get current GitHub IPs
curl -s https://api.github.com/meta | jq -r '.git[]' | grep -E '^[0-9]'

# Returns ranges like:
# 192.30.252.0/22
# 185.199.108.0/22
# 140.82.112.0/20
# etc.
```

**Important**: IPs change periodically. Firewall script needs regular updates or dynamic fetching.

## npm Registry IPs

npm uses Cloudflare CDN, making IP allowlisting challenging:

```bash
# Cloudflare IP ranges (large)
curl -s https://www.cloudflare.com/ips-v4 | head -10

# Alternative: Allow specific domains via DNS resolution
# More complex to implement with iptables
```

## Non-Root User Pattern

Anthropic uses `node` user. To implement in custom container:

```dockerfile
# Create user
RUN useradd -m -s /bin/bash developer && \
    usermod -aG sudo developer && \
    echo "developer ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Switch to user
USER developer

# Adjust paths in devcontainer.json
# /root/.claude → /home/developer/.claude
```

**Trade-offs**:
- Better security posture
- More complex permission management
- Need sudo for system operations
- Language tools may need user-specific config

## Bash History Persistence

Nice quality-of-life feature from Anthropic:

```json
"mounts": [
  "source=project-bash-history,target=/root/.bash_history,type=volume"
]
```

**Benefits**:
- Command history persists across container rebuilds
- Improves developer experience
- Minimal overhead

**Implementation**: Just add the mount, no Dockerfile changes needed.

## Summary of Anthropic Best Practices

✅ **Adopt for all sandboxes**:
- Claude Code VS Code extension pre-configured
- Persistent Claude config volume
- GitHub CLI feature integration
- GitHub credential forwarding via mount
- Environment variables for Claude Code

⚠️ **Consider based on requirements**:
- Firewall with IP allowlisting (security vs complexity)
- Non-root user (security vs simplicity)
- Bash history persistence (nice to have)

❌ **Skip for language-specific sandboxes**:
- Microsoft devcontainer base image (use language-specific base)
- Node.js-only focus (support multiple languages)

## References

- **Official Anthropic Devcontainer**: https://github.com/anthropics/claude-code/tree/main/.devcontainer
- **Dockerfile**: https://github.com/anthropics/claude-code/blob/main/.devcontainer/Dockerfile
- **devcontainer.json**: https://github.com/anthropics/claude-code/blob/main/.devcontainer/devcontainer.json
- **Firewall Script**: https://github.com/anthropics/claude-code/blob/main/.devcontainer/init-firewall.sh
- **GitHub Meta API**: https://api.github.com/meta (for current IPs)
- **Dev Containers Specification**: https://containers.dev/

## Evolution from Anthropic to CSRRS

The CSRRS implementation evolved Anthropic's patterns for Rust development:

1. **Started with**: Anthropic's security-focused, Node.js-centric container
2. **Adapted for**: Rust systems programming with large reference materials
3. **Added**: Multi-stage builds, language runtime, comprehensive tooling
4. **Simplified**: Removed firewall, used root user
5. **Enhanced**: Better documentation, build scripts, README templates
6. **Generalized**: Created language-agnostic patterns for any project

This skill codifies that evolution for reuse across any language ecosystem.
