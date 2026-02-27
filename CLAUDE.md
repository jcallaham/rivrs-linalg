# CLAUDE.md - rivrs-linalg

This file provides guidance to Claude Code when working at the top level of the rivrs-linalg repository.

## Repository Overview

This repository contains numerical linear algebra implementations for the Rivrs symbolic-numeric framework. Each algorithm domain is a standalone Rust crate (`rivrs-sparse`, `rivrs-control`, etc.). Thin facade crates (`rivrs-linalg`, `rivrs`) re-export domain crates so users can depend at whatever granularity they prefer.

**License**: Apache-2.0

## Crate Structure

| Crate | Location | Role |
|-------|----------|------|
| `rivrs-sparse` | `sparse/` | Real implementation — sparse symmetric indefinite solver (SSIDS) |
| `rivrs-control` | `control/` | Real implementation — control systems algorithms (not yet published) |
| `rivrs-linalg` | `crates/rivrs-linalg/` | Facade — re-exports domain crates as `rivrs_linalg::sparse`, etc. |
| `rivrs` | `crates/rivrs/` | Umbrella — re-exports `rivrs-linalg` as `rivrs::linalg` |

The facade/umbrella crates contain minimal code (re-exports + documentation). All real implementation lives in the domain crates.

## When to Use This File vs Crate-Specific CLAUDE.md

**Use this file when:**
- Working at repository root level
- Creating or modifying facade/umbrella crates
- Managing shared infrastructure (references/, CI/CD)
- Planning repository-wide changes

**Use crate-specific CLAUDE.md when:**
- Implementing algorithms in `sparse/` -> use `sparse/CLAUDE.md`
- Implementing algorithms in `control/` -> use `control/CLAUDE.md`

## Domain Crates

### sparse/ — rivrs-sparse

**Status**: Feature-complete. Preparing for 0.1.0 crates.io release.
**Focus**: SPRAL-inspired SSIDS (Sparse Symmetric Indefinite Direct Solver)

See [sparse/CLAUDE.md](sparse/CLAUDE.md) for complete guidance.

### control/ — rivrs-control

**Status**: Sylvester solvers complete. Lyapunov/Riccati planned. Not yet published.
**Focus**: SLICOT-inspired algorithms using faer

See [control/CLAUDE.md](control/CLAUDE.md) for complete guidance.

## Release Plan

1. **Publish** `rivrs-sparse` 0.1.0 (the existing `sparse/` project, as-is).
2. **Publish** `rivrs-linalg` 0.1.0 as a facade re-exporting `rivrs-sparse`.
3. **Publish** `rivrs` 0.1.0 as an umbrella re-exporting `rivrs-linalg`.
4. **Add domain crates** over time and wire them into the facades.
5. **Optionally consolidate** into a monolithic `rivrs-linalg` later if facade indirection becomes undesirable.

## Clean Room Implementation - Repository-Wide Policy

**Critical**: All implementations must maintain clean room status to preserve Apache-2.0 licensing.

**NEVER do this:**
- Read or copy GPL-licensed source code during algorithm implementation
- Read proprietary/restrictive source code (HSL libraries)
- Implement algorithms by translating GPL code line-by-line

**ALWAYS do this:**
- Implement from academic papers and textbooks
- Consult permissively-licensed reference code (LAPACK, SLICOT-Reference, SPRAL)
- Use documentation and test cases from any library (non-copyrightable facts)
- Document all academic sources used for each implementation
- Cite reference implementations consulted

### Crate-Specific Guidelines

**control/** (SLICOT-based):
- SLICOT-Reference (BSD-3): Consult freely
- SLICOT documentation (slicot/doc/): Consult freely
- SLICOT test cases (slicot/examples/): Use for validation
- SLICOT GPL source (slicot/src/*.f): NEVER read during implementation
- LAPACK source (BSD-3): Consult freely

**sparse/** (SPRAL-based):
- SPRAL (BSD-3): Consult freely - primary reference
- HSL documentation: Consult for algorithm descriptions
- HSL source code: NEVER read (proprietary/restrictive)
- Academic papers: Primary source for algorithms

## Development Workflow

### Working on a Domain Crate

```bash
# cd into the crate directory
cd sparse/    # or cd control/

# Standard Rust workflow
cargo build
cargo test
cargo bench
```

### Managing Shared References

Reference implementations in `references/` are accessible to all crates:
- `faer-rs/` - Rust linear algebra library
- `lapack/` - Reference LAPACK (BSD)
- `SLICOT-Reference/` - BSD-3 SLICOT
- `slicot/` - GPL SLICOT (docs/tests only, NOT source)
- `spral/` - SPRAL sparse library (BSD-3)

## Git Workflow

### Commit Practices
- Commit frequently after logical units
- Commit after tests pass or milestones are verified
- Use clear, descriptive commit messages
- Include Co-Authored-By for AI assistance when appropriate

### Branch Strategy
- `main` branch for stable development
- Feature branches for major changes

## Documentation Standards

### Algorithm Documentation
Every algorithm implementation MUST include:
1. Mathematical description of the operation
2. **Academic sources cited** (papers, textbooks, standards)
3. **Reference implementations consulted** (LAPACK routines, SLICOT routines, SPRAL functions)
4. Numerical stability characteristics
5. Complexity analysis (time, space)
6. Usage examples
7. Cross-references to equivalent routines in other libraries

### Clean Room Audit Trail
Each implementation should be traceable to non-GPL sources:
- Cite specific papers, textbook sections, or permissively-licensed code
- Document decision points where algorithm choices were made
- Note any deviations from reference implementations

## License and Attribution

**License**: Apache-2.0 (see [LICENSE](LICENSE))

**Attribution**: See NOTICE for complete attribution to:
- Academic sources (papers, textbooks)
- LAPACK (BSD-3-Clause)
- SLICOT-Reference (BSD-3-Clause)
- SPRAL (BSD-3-Clause)
- MatrixEquations.jl (MIT, for validation)

Each crate's implementation cites specific sources used.
