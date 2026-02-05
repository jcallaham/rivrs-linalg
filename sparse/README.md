# Sparse Linear Algebra Solvers

This directory contains sparse linear algebra solver implementations for rivrs-linalg.

## Overview

Implementations of sparse solvers inspired by SPRAL, using modern Rust with `faer` for high-performance linear algebra. All algorithms are implemented from academic papers and permissively-licensed reference code (SPRAL, LAPACK).

## Planned Algorithms

### SSIDS — Sparse Symmetric Indefinite Direct Solver

| Component | Description | Status |
|---|---|---|
| Symbolic analysis | Ordering, elimination tree, symbolic factorization | Planned |
| Numeric factorization | LDL^T with APTP pivoting | Planned |
| Triangular solve | Forward/backward substitution | Planned |

Target features:
- A Posteriori Threshold Pivoting (APTP) for robust indefinite factorization
- Simplicial factorization (Phase 1), supernodal optimization (Phase 2)
- Integration with faer's sparse infrastructure (CSC, elimination trees, AMD ordering)

## Quick Start

```bash
cd sparse
cargo build
cargo test
```

## Development

See [CLAUDE.md](CLAUDE.md) for detailed development guidance.

### Docker Development Environment

Use the centralized Docker environment at the repository root:

```bash
cd ../docker
./build.sh
./run.sh
# Once inside container:
cd sparse/
```

## Documentation

- [CLAUDE.md](CLAUDE.md) - Development guidance and clean room methodology
- [specs/](specs/) - Feature specifications and planning documents
- Parent [NOTICE](../NOTICE) - Academic sources and licensing attribution

## Roadmap

**Near-term:**
- Simplicial LDL^T factorization with APTP
- Symbolic analysis (ordering, elimination tree)

**Future:**
- Supernodal optimization
- GPU offload for large frontal matrices
- Iterative refinement
- Python bindings via PyO3

## License

Part of rivrs-linalg, licensed under Apache-2.0. See [../LICENSE](../LICENSE).
