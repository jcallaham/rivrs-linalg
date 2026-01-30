# VS Code Devcontainer Configuration Examples

Complete devcontainer.json configurations for different language ecosystems, all with Claude Code integration.

## Python Data Science Configuration

```json
{
	"name": "Python Data Science Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "python-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"ms-python.python",
				"ms-python.vscode-pylance",
				"ms-python.debugpy",
				"ms-toolsai.jupyter",
				"ms-toolsai.vscode-jupyter-cell-tags",
				"ms-toolsai.jupyter-renderers",
				"charliermarsh.ruff"
			],
			"settings": {
				"editor.tabSize": 4,
				"editor.insertSpaces": true,
				"editor.formatOnSave": true,
				"editor.rulers": [88],
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"python.defaultInterpreterPath": "/usr/bin/python3",
				"python.formatting.provider": "black",
				"python.linting.enabled": true,
				"python.linting.pylintEnabled": false,
				"python.linting.flake8Enabled": true,
				"python.testing.pytestEnabled": true,
				"[python]": {
					"editor.defaultFormatter": "charliermarsh.ruff"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"PYTHONPATH": "/workspace/project/src",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "python3 --version && pip3 --version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=python-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Rust Systems Programming Configuration

```json
{
	"name": "Rust Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "rust-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"rust-lang.rust-analyzer",
				"vadimcn.vscode-lldb",
				"serayuzgur.crates",
				"tamasfe.even-better-toml",
				"dustypomerleau.rust-syntax"
			],
			"settings": {
				"editor.tabSize": 8,
				"editor.insertSpaces": false,
				"editor.detectIndentation": false,
				"editor.rulers": [80, 100],
				"editor.formatOnSave": true,
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"rust-analyzer.checkOnSave.command": "clippy",
				"rust-analyzer.cargo.features": "all",
				"rust-analyzer.inlayHints.parameterHints.enable": false,
				"[rust]": {
					"editor.defaultFormatter": "rust-lang.rust-analyzer"
				},
				"[toml]": {
					"editor.defaultFormatter": "tamasfe.even-better-toml"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"RUSTFLAGS": "-C target-cpu=native",
		"CARGO_BUILD_JOBS": "8",
		"RUSTC_WRAPPER": "sccache",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "rustc --version && cargo --version && gh --version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=rust-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Node.js/TypeScript Configuration

```json
{
	"name": "Node.js Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "node-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"dbaeumer.vscode-eslint",
				"esbenp.prettier-vscode",
				"ms-vscode.vscode-typescript-next",
				"orta.vscode-jest",
				"christo pher.vscode-test-adapter",
				"firsttris.vscode-jest-runner"
			],
			"settings": {
				"editor.tabSize": 2,
				"editor.insertSpaces": true,
				"editor.formatOnSave": true,
				"editor.rulers": [80, 120],
				"editor.codeActionsOnSave": {
					"source.fixAll.eslint": "explicit"
				},
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"eslint.validate": [
					"javascript",
					"javascriptreact",
					"typescript",
					"typescriptreact"
				],
				"[javascript]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[typescript]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[json]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"NODE_ENV": "development",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "node --version && npm --version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=node-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Go Configuration

```json
{
	"name": "Go Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "go-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"golang.go"
			],
			"settings": {
				"editor.tabSize": 4,
				"editor.insertSpaces": false,
				"editor.formatOnSave": true,
				"editor.rulers": [80, 120],
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"go.toolsManagement.autoUpdate": true,
				"go.useLanguageServer": true,
				"go.lintTool": "golangci-lint",
				"go.lintOnSave": "workspace",
				"go.formatTool": "gofmt",
				"go.testFlags": ["-v", "-race"],
				"[go]": {
					"editor.defaultFormatter": "golang.go"
				},
				"[go.mod]": {
					"editor.defaultFormatter": "golang.go"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"GOPATH": "/root/go",
		"GO111MODULE": "on",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "go version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=go-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Java/Maven Configuration

```json
{
	"name": "Java Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "java-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"vscjava.vscode-java-pack",
				"vscjava.vscode-maven",
				"vscjava.vscode-java-test",
				"vscjava.vscode-java-debug",
				"redhat.java"
			],
			"settings": {
				"editor.tabSize": 4,
				"editor.insertSpaces": true,
				"editor.formatOnSave": true,
				"editor.rulers": [80, 120],
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"java.configuration.runtimes": [
					{
						"name": "JavaSE-17",
						"path": "/usr/lib/jvm/java-17-openjdk-arm64"
					}
				],
				"java.format.settings.url": ".vscode/java-formatter.xml",
				"java.saveActions.organizeImports": true,
				"[java]": {
					"editor.defaultFormatter": "redhat.java"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"JAVA_HOME": "/usr/lib/jvm/java-17-openjdk-arm64",
		"MAVEN_OPTS": "-Xmx2048m",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "java -version && mvn --version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=java-claude-config,target=/root/.claude,type=volume"
	]
}
```

## React/Frontend Configuration

```json
{
	"name": "React Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "react-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"dbaeumer.vscode-eslint",
				"esbenp.prettier-vscode",
				"dsznajder.es7-react-js-snippets",
				"styled-components.vscode-styled-components",
				"bradlc.vscode-tailwindcss",
				"ms-vscode.vscode-typescript-next"
			],
			"settings": {
				"editor.tabSize": 2,
				"editor.insertSpaces": true,
				"editor.formatOnSave": true,
				"editor.rulers": [80, 120],
				"editor.codeActionsOnSave": {
					"source.fixAll.eslint": "explicit"
				},
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"emmet.includeLanguages": {
					"javascript": "javascriptreact",
					"typescript": "typescriptreact"
				},
				"[javascript]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[javascriptreact]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[typescript]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[typescriptreact]": {
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"NODE_ENV": "development",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "node --version && npm --version && claude --version && echo '✓ Development environment ready'",

	"forwardPorts": [3000, 5173, 8080],
	"portsAttributes": {
		"3000": {
			"label": "Application",
			"onAutoForward": "notify"
		}
	},

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=react-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Multi-Language Polyglot Configuration

```json
{
	"name": "Polyglot Development Container",
	"dockerComposeFile": "docker-compose.yml",
	"service": "polyglot-dev",
	"workspaceFolder": "/workspace/project",
	"shutdownAction": "stopCompose",

	"features": {
		"ghcr.io/devcontainers/features/github-cli:1": {
			"version": "latest"
		}
	},

	"customizations": {
		"vscode": {
			"extensions": [
				"anthropic.claude-code",
				"ms-python.python",
				"ms-python.vscode-pylance",
				"rust-lang.rust-analyzer",
				"vadimcn.vscode-lldb",
				"dbaeumer.vscode-eslint",
				"esbenp.prettier-vscode"
			],
			"settings": {
				"editor.formatOnSave": true,
				"files.insertFinalNewline": true,
				"files.trimTrailingWhitespace": true,
				"[python]": {
					"editor.tabSize": 4,
					"editor.insertSpaces": true
				},
				"[rust]": {
					"editor.tabSize": 8,
					"editor.insertSpaces": false,
					"editor.defaultFormatter": "rust-lang.rust-analyzer"
				},
				"[javascript]": {
					"editor.tabSize": 2,
					"editor.insertSpaces": true,
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				},
				"[typescript]": {
					"editor.tabSize": 2,
					"editor.insertSpaces": true,
					"editor.defaultFormatter": "esbenp.prettier-vscode"
				}
			}
		}
	},

	"remoteUser": "root",

	"remoteEnv": {
		"PYTHONPATH": "/workspace/project/backend",
		"RUSTFLAGS": "-C target-cpu=native",
		"NODE_ENV": "development",
		"NODE_OPTIONS": "--max-old-space-size=4096",
		"CLAUDE_CONFIG_DIR": "/root/.claude"
	},

	"postCreateCommand": "python3 --version && rustc --version && node --version && claude --version && echo '✓ Development environment ready'",

	"mounts": [
		"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
		"source=polyglot-claude-config,target=/root/.claude,type=volume"
	]
}
```

## Common Configuration Patterns

### Required for All Configurations

**Claude Code Extension**:
```json
"extensions": [
	"anthropic.claude-code",
	// ... other extensions
]
```

**Environment Variables**:
```json
"remoteEnv": {
	"NODE_OPTIONS": "--max-old-space-size=4096",
	"CLAUDE_CONFIG_DIR": "/root/.claude"
}
```

**GitHub CLI Feature**:
```json
"features": {
	"ghcr.io/devcontainers/features/github-cli:1": {
		"version": "latest"
	}
}
```

**Persistent Mounts**:
```json
"mounts": [
	"source=${localEnv:HOME}/.config/gh,target=/root/.config/gh,type=bind,consistency=cached",
	"source=project-claude-config,target=/root/.claude,type=volume"
]
```

### Editor Settings Best Practices

**Universal Settings**:
```json
"editor.formatOnSave": true,
"files.insertFinalNewline": true,
"files.trimTrailingWhitespace": true,
"editor.rulers": [80, 120]
```

**Language-Specific Overrides**:
```json
"[language]": {
	"editor.tabSize": 4,
	"editor.insertSpaces": true,
	"editor.defaultFormatter": "extension.id"
}
```

### Port Forwarding (for Web Apps)

```json
"forwardPorts": [3000, 8080, 5173],
"portsAttributes": {
	"3000": {
		"label": "Frontend",
		"onAutoForward": "notify"
	},
	"8080": {
		"label": "Backend API",
		"onAutoForward": "silent"
	}
}
```

### Post-Create Commands

**Verification Script**:
```json
"postCreateCommand": "language --version && claude --version && echo '✓ Ready'"
```

**Setup Script**:
```json
"postCreateCommand": "bash .devcontainer/setup.sh"
```

**Chain Commands**:
```json
"postCreateCommand": "npm install && npm run build && echo '✓ Ready'"
```

## Customization Guide

When creating your own configuration:

1. **Service Name**: Match `"service"` to docker-compose.yml
2. **Workspace Folder**: Match Dockerfile's WORKDIR
3. **Extensions**: Add language-specific extensions
4. **Settings**: Configure formatters, linters, etc.
5. **Environment Variables**: Add language-specific vars
6. **Post-Create**: Verify installations
7. **Mounts**: Always include gh config + Claude config
8. **Ports**: Forward if running web services

## Volume Naming Convention

Use project-specific volume names:
```json
"source=projectname-claude-config,target=/root/.claude,type=volume"
```

This prevents volume conflicts between projects.

## Remote User Consideration

These configurations use `"remoteUser": "root"` for simplicity. For production or shared environments, consider creating a non-root user:

```dockerfile
RUN useradd -m -s /bin/bash developer
USER developer
```

Then in devcontainer.json:
```json
"remoteUser": "developer"
```

Adjust paths accordingly (`/home/developer/.claude` instead of `/root/.claude`).

## Platform-Specific Settings

### ARM64 (Apple Silicon)
No special settings needed - VS Code detects platform automatically.

### AMD64 (Intel/AMD)
No special settings needed - platform specified in docker-compose.yml.

### Multi-Platform
VS Code uses the image built for your architecture.

## Troubleshooting

**Extensions not loading**:
- Rebuild container: F1 → "Dev Containers: Rebuild Container"
- Check extension IDs are correct
- Ensure extensions support remote containers

**Post-create command fails**:
- Check syntax (bash commands, proper escaping)
- Verify commands exist in container
- Test commands manually after container starts

**Mounts not working**:
- Check paths exist on host (`~/.config/gh`)
- Verify volume names are unique
- Ensure Docker Desktop has file sharing permissions

**GitHub auth not working**:
- Run `gh auth login` on host
- Verify `gh auth status` shows authenticated
- Check mount path is correct
- Restart container after authenticating
