# rivrs-linalg — Numerical Linear Algebra for Rivrs

A scientific computing library providing numerical linear algebra implementations
for the [Rivrs](https://github.com/jcallaham/rivrs) symbolic-numeric framework.

**Current Status**: Phase 1 - Modular development structure
**License**: Apache-2.0

## Project Structure

This repository is organized as a monorepo with isolated project directories for different algorithm domains. This Phase 1 structure enables independent development of distinct numerical algorithms before eventual integration into a unified Cargo workspace.

```
rivrs-linalg/
├── control/          # Control systems algorithms (Sylvester, Lyapunov, Riccati)
├── sparse/           # Sparse direct solvers
├── references/       # Shared reference implementations (LAPACK, SPRAL, faer-rs)
├── NOTICE            # Attribution and licensing information
└── LICENSE           # Apache-2.0 license
```

Each project directory (`control/`, `sparse/`) is a complete, standalone Rust project with:
- Own `Cargo.toml` (immediately buildable with `cargo build`)
- Own `CLAUDE.md` (domain-specific development guidance)
- Own `README.md` (algorithm descriptions and status)
- Own `Dockerfile` and `.devcontainer/` (isolated development environment)
- Own `.claude/` and `.specify/` (coding agent and spec-kit configuration)
- Own `specs/` (feature specifications and research)

## Current Projects

### [control/](control/) - Control Systems Algorithms

**Status**: Active development
**Algorithms**: Sylvester equations (continuous/discrete)

Clean room implementations of control systems algorithms similar to SLICOT, using modern Rust with `faer` for high-performance linear algebra.

**Completed:**
- ✅ Continuous Sylvester solver (AX + XB = C) - Bartels-Stewart algorithm
- ✅ Discrete Sylvester solver (AXB + X = C)
- ✅ Blocked Level-3 BLAS variant for continuous case
- ✅ Overflow prevention and condition estimation

**Roadmap:**
- 📋 Lyapunov equation solvers
- 📋 Algebraic Riccati equation solvers

See [control/README.md](control/README.md) for details.

## Development Workflow

### Working on a Specific Domain

```bash
# Control systems work
cd control/
cargo build
cargo test
cargo bench

# Sparse solvers work
cd sparse/
cargo build
cargo test
```

Each directory is isolated and can be developed independently.

### Docker Development Environments

Each domain has its own Docker environment with domain-specific reference materials:

```bash
cd control/docker/
./build.sh
./run.sh

# Or for sparse solvers
cd sparse/docker/
./build.sh
./run.sh
```

### Shared References

Reference implementations are stored in `references/` and are accessible to all projects:
- **faer-rs/**: Pure Rust linear algebra library
- **lapack/**: Reference LAPACK (BSD-licensed)
- **SLICOT-Reference/**: BSD-3-licensed SLICOT implementation
- **slicot/**: GPL SLICOT (documentation and test cases only)
- **spral/**: SPRAL sparse library (to be added, BSD-3-licensed)

## Clean Room Implementation

**Critical**: All implementations are clean room to maintain Apache-2.0 licensing.

### Permitted Sources
- ✅ Academic papers and textbooks
- ✅ Permissively-licensed reference code (LAPACK, SLICOT-Reference, SPRAL)
- ✅ Documentation and test cases from any library

### Prohibited Sources
- ❌ GPL SLICOT Fortran source code (slicot/src/*.f)
- ❌ Proprietary HSL library source code
- ❌ Any GPL-licensed implementation source code

See [NOTICE](NOTICE) for complete attribution and licensing information.

## Migration Path: Phase 1 → Workspace

**Phase 1 (Current)**: Isolated projects for independent experimentation
- Each domain is a complete standalone project
- Enables rapid iteration and easy pivoting
- Simple Docker isolation per domain
- Can abandon or restructure without affecting others

**Phase 2 (Future)**: Cargo workspace with shared utilities
```
rivrs-linalg/
├── Cargo.toml          # Workspace manifest
├── rivrs-core/         # Shared utilities, traits, errors
├── rivrs-control/      # Control systems (renamed from control/)
├── rivrs-sparse/       # Sparse solvers (renamed from sparse/)
└── references/         # Shared references
```

Migration will be straightforward:
1. Create `rivrs-core/` with shared code
2. Move `control/` → `rivrs-control/`, update dependencies
3. Move `sparse/` → `rivrs-sparse/`, update dependencies
4. Add workspace `Cargo.toml`
5. Git history preserved via file moves

**Phase 3 (Integration)**: Rivrs dependency

Eventually, rivrs-linalg will be used by the main Rivrs symbolic-numeric framework via path dependency or published crates.

## Contributing

See individual project CLAUDE.md files for domain-specific development guidance:
- [control/CLAUDE.md](control/CLAUDE.md) - Control systems algorithms
- [sparse/CLAUDE.md](sparse/CLAUDE.md) - Sparse solvers

**Key Requirements:**
- Maintain clean room status (cite academic sources, avoid GPL code)
- Write comprehensive tests
- Document numerical stability characteristics
- Benchmark against reference implementations

## License

Licensed under Apache License, Version 2.0 ([LICENSE](LICENSE))

## Attribution

See [NOTICE](NOTICE) for complete attribution to academic sources and reference implementations:
- LAPACK (BSD-3-Clause)
- SLICOT-Reference (BSD-3-Clause)
- SPRAL (BSD-3-Clause, planned)
- Academic papers and textbooks

Each implementation cites the specific sources used.
