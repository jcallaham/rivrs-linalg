# Control Systems Algorithms

This directory contains control systems algorithm implementations for rivrs-linalg.

## Overview

Implementations of control systems algorithms similar to those in SLICOT, using modern Rust with `faer` for high-performance linear algebra. All algorithms are implemented from academic papers and permissively-licensed reference code (LAPACK, SLICOT-Reference).

## Current Algorithms

### Sylvester Equation Solvers

| Equation | Function | Algorithm |
|---|---|---|
| Continuous: `AX + XB = C` | `solve_continuous` | Bartels-Stewart (Schur-Schur) |
| Discrete: `AXB + X = C` | `solve_discrete` | Modified Bartels-Stewart |

Features:
- Eigenvalue separation estimation for singularity detection
- Overflow prevention via scaling
- Blocked Level-3 BLAS triangular solver (continuous, N > 64)
- Pre-computed Schur form APIs

## Quick Start

```bash
cd control
cargo build
cargo test
cargo bench
```

## Examples

```rust
use rivrs_control::sylvester::solve_continuous;
use faer::mat;

let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
let c = mat![[4.0, 5.0], [6.0, 12.0f64]];

let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
let x = &result.solution * (1.0 / result.scale);
```

See `examples/` directory for more.

## Development

See [CLAUDE.md](CLAUDE.md) for detailed development guidance.

### Docker Development Environment

Use the centralized Docker environment at the repository root:

```bash
cd ../docker
./build.sh
./run.sh
# Once inside container:
cd control/
```

This provides a unified environment with:
- Rust toolchain
- Claude Code
- Reference materials (LAPACK, SLICOT-Reference, faer-rs)
- GitHub CLI integration
- Access to all monorepo projects

## Documentation

- [CLAUDE.md](CLAUDE.md) - Development guidance and clean room methodology
- [specs/](specs/) - Feature specifications and planning documents
- Parent [NOTICE](../NOTICE) - Academic sources and licensing attribution

## Roadmap

**Near-term:**
- Lyapunov equation solvers
- Continuous and discrete Riccati solvers

**Future:**
- Generalized Sylvester/Lyapunov equations
- State-space canonical forms
- Model reduction and balancing

## License

Part of rivrs-linalg, licensed under Apache-2.0. See [../LICENSE](../LICENSE).
