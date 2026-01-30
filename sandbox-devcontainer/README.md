# Sandbox Devcontainer Builder Skill

A comprehensive skill for creating isolated Docker devcontainers with Claude Code integration for any programming language or framework.

## Quick Start

Invoke this skill when you need to create a sandboxed development environment:

```
/sandbox-devcontainer Create a sandbox image for this project
```

Claude will:
1. Analyze your project structure and detect the primary language
2. Identify static reference materials (untracked files) to embed
3. Generate appropriate Dockerfile, devcontainer.json, and docker-compose.yml
4. Create comprehensive documentation
5. Set up Claude Code integration

## What This Skill Provides

### Core Guidance (SKILL.md)

- Complete workflow for creating sandboxed devcontainers
- Language-agnostic approach supporting Rust, Python, Node.js, Go, Java, etc.
- Claude Code integration patterns
- Static reference material embedding
- GitHub authentication via HTTPS + CLI
- Clean room implementation guidelines (for license-sensitive projects)

### Reference Files

Detailed documentation for specific aspects:

- **`dockerfile-templates.md`** - Complete Dockerfiles for Python, Rust, Node.js, Go, Java, multi-language
- **`devcontainer-configs.md`** - VS Code devcontainer.json examples per language
- **`github-auth-patterns.md`** - HTTPS authentication with GitHub CLI (not SSH)
- **`anthropic-reference.md`** - Learnings from Anthropic's official implementation
- **`architecture-patterns.md`** - Multi-stage build strategies (TODO: add if needed)
- **`language-runtimes.md`** - Runtime installation guides (TODO: add if needed)
- **`reference-strategies.md`** - Handling large reference materials (TODO: add if needed)

### Examples

Real-world implementations:

- **`csrrs-implementation.md`** - Complete Rust scientific computing example with 710MB references
- **`README-template.md`** - Template for docker/README.md documentation
- (Add more examples as you create sandboxes for other projects)

### Scripts

Automation tools:

- **`analyze-project.sh`** - Automated project analysis (language detection, static references, size estimation)
- (Add more scripts like generate-dockerfile.sh, test-container.sh as needed)

## File Organization

```
sandbox-devcontainer/
├── SKILL.md                      # Main skill content (always loaded)
├── README.md                     # This file
├── references/                   # Detailed patterns (loaded as needed)
│   ├── dockerfile-templates.md
│   ├── devcontainer-configs.md
│   ├── github-auth-patterns.md
│   └── anthropic-reference.md
├── examples/                     # Working examples
│   ├── csrrs-implementation.md
│   └── README-template.md
└── scripts/                      # Utility scripts
    └── analyze-project.sh
```

## Typical Usage

### Scenario 1: New Project

```
User: "Create a sandbox devcontainer for my Python data science project"

Claude (using this skill):
1. Runs analyze-project.sh to detect Python
2. Identifies any untracked datasets/docs
3. Generates Dockerfile based on Python template
4. Creates devcontainer.json with Python extensions
5. Sets up docker-compose.yml with pip cache volumes
6. Creates build.sh and .dockerignore
7. Generates comprehensive README.md
```

### Scenario 2: Existing Project with References

```
User: "Add Claude Code to my Rust project with embedded SLICOT references"

Claude (using this skill):
1. Detects Rust via Cargo.toml
2. Finds slicot/ directory (710MB untracked)
3. Creates multi-stage Dockerfile (base → references → dev)
4. Configures devcontainer.json with Rust extensions
5. Documents clean room rules in README
6. Sets up GitHub HTTPS authentication
```

### Scenario 3: Multi-Language Project

```
User: "Create a sandbox for my full-stack app (Python backend + React frontend)"

Claude (using this skill):
1. Detects both Python and Node.js
2. Uses polyglot Dockerfile template
3. Installs Python + Node.js (Node.js serves both frontend and Claude Code)
4. Configures extensions for both languages
5. Sets up port forwarding for web app
6. Creates appropriate cache volumes for pip and npm
```

## Key Principles Codified

1. **Complete Isolation**: No host dependencies, all runtimes embedded
2. **Claude Code Integration**: Always included for AI assistance
3. **HTTPS + GitHub CLI**: Simpler than SSH, auto-forwarded by VS Code
4. **Static Reference Embedding**: Copy untracked materials into image
5. **Multi-Stage Builds**: Efficient caching and layer separation
6. **Persistent Volumes**: Cache packages and Claude config
7. **Comprehensive Documentation**: Always generate docker/README.md

## When NOT to Use This Skill

- Simple projects that don't need containerization
- Projects that should mount host directories (not isolated)
- When you just want to add Claude Code to existing devcontainer (use simpler approach)

## Related Skills

This skill complements:

- **Plugin development skills**: For creating Claude Code plugins that run IN these sandboxes
- **Frontend/backend skills**: Provide domain expertise for code WRITTEN in these sandboxes
- **Testing/CI skills**: For validating code in these sandboxes

## Evolution

This skill was created by codifying the CSRRS devcontainer implementation process:

1. Started with Anthropic's official Node.js devcontainer
2. Adapted for Rust scientific computing with large references
3. Generalized to support any language/framework
4. Captured as reusable skill for future projects

## Contributing

As you create sandboxes for new projects:

1. Add language-specific Dockerfile to `references/dockerfile-templates.md`
2. Add devcontainer config to `references/devcontainer-configs.md`
3. Document implementation in `examples/[project]-implementation.md`
4. Update this README with new use cases

## Version History

- **v1.0.0** (2024-01-29): Initial skill creation based on CSRRS implementation
  - Python, Rust, Node.js, Go, Java templates
  - GitHub HTTPS authentication patterns
  - Multi-stage build guidance
  - CSRRS real-world example
