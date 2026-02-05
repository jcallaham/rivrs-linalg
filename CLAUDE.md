# CLAUDE.md - rivrs-linalg Monorepo

This file provides guidance to Claude Code when working at the top level of the rivrs-linalg monorepo.

## Repository Overview

rivrs-linalg is a monorepo containing numerical linear algebra implementations for the Rivrs symbolic-numeric framework. The repository is organized into isolated project directories for different algorithm domains.

**Current Structure**: Phase 1 - Isolated projects for modular development
**Future Structure**: Phase 2 - Cargo workspace with shared utilities
**License**: Apache-2.0

## Directory Structure

```
rivrs-linalg/
├── control/          # Control systems algorithms (active development)
│   ├── CLAUDE.md     # Domain-specific guidance for control work
│   ├── Cargo.toml    # Standalone Rust project
│   ├── src/          # Implementation
│   ├── specs/        # Feature specifications
│   └── docker/       # Isolated development environment
├── references/       # Shared reference implementations
│   ├── faer-rs/      # Rust linear algebra library
│   ├── lapack/       # Reference LAPACK (BSD)
│   ├── SLICOT-Reference/  # BSD-3 SLICOT
│   ├── slicot/       # GPL SLICOT (docs/tests only)
│   └── spral/        # To be added (BSD-3)
├── NOTICE            # Attribution and licensing
├── LICENSE           # Apache-2.0 license
└── CLAUDE.md         # This file (monorepo-level guidance)
```

## When to Use This File vs Domain-Specific CLAUDE.md

**Use this file when:**
- Working at repository root level
- Setting up new project directories
- Managing shared infrastructure (references/, CI/CD)
- Planning repository-wide changes
- Understanding the overall structure

**Use domain-specific CLAUDE.md when:**
- Implementing algorithms in `control/` → use `control/CLAUDE.md`
- Domain-specific development questions

Each domain directory is self-contained with its own development guidance.

## Project Domains

### control/ - Control Systems Algorithms

**Status**: Active development
**Focus**: SLICOT-inspired algorithms using faer
**Current**: Sylvester equation solvers (continuous, discrete)
**Roadmap**: Lyapunov, Riccati, state-space analysis

See [control/CLAUDE.md](control/CLAUDE.md) for complete guidance.

### sparse/ - Sparse Solvers

**Status**: Phase 0 - Scaffolding and literature review
**Focus**: SPRAL-inspired SSIDS (Sparse Symmetric Indefinite Direct Solver)
**Current**: Project skeleton, APTP algorithm research
**Roadmap**: Symbolic analysis, simplicial LDL^T, supernodal optimization

See [sparse/CLAUDE.md](sparse/CLAUDE.md) for complete guidance.

## Clean Room Implementation - Repository-Wide Policy

**Critical**: All implementations must maintain clean room status to preserve Apache-2.0 licensing.

### Universal Rules (Apply to All Domains)

**NEVER do this:**
- ❌ Read or copy GPL-licensed source code during algorithm implementation
- ❌ Read proprietary/restrictive source code (HSL libraries)
- ❌ Implement algorithms by translating GPL code line-by-line

**ALWAYS do this:**
- ✅ Implement from academic papers and textbooks
- ✅ Consult permissively-licensed reference code (LAPACK, SLICOT-Reference, SPRAL)
- ✅ Use documentation and test cases from any library (non-copyrightable facts)
- ✅ Document all academic sources used for each implementation
- ✅ Cite reference implementations consulted

### Domain-Specific Guidelines

**control/** (SLICOT-based):
- ✅ SLICOT-Reference (BSD-3): Consult freely
- ✅ SLICOT documentation (slicot/doc/): Consult freely
- ✅ SLICOT test cases (slicot/examples/): Use for validation
- ❌ SLICOT GPL source (slicot/src/*.f): NEVER read during implementation
- ✅ LAPACK source (BSD-3): Consult freely

**sparse/** (SPRAL-based):
- ✅ SPRAL (BSD-3): Consult freely - primary reference
- ✅ HSL documentation: Consult for algorithm descriptions
- ❌ HSL source code: NEVER read (proprietary/restrictive)
- ✅ Academic papers: Primary source for algorithms

## Development Workflow

### Working on a Specific Domain

```bash
# Always cd into the domain directory
cd control/    # or cd sparse/

# Standard Rust workflow
cargo build
cargo test
cargo bench

# Docker development (isolated environment)
cd docker/
./build.sh
./run.sh
```

### Adding New Domains

When adding a new algorithm domain (e.g., `dense/`, `iterative/`):

1. Create new directory: `mkdir new-domain/`
2. Copy template structure from existing domain
3. Create domain-specific `CLAUDE.md`
4. Set up `Cargo.toml` as standalone project
5. Create `README.md` with domain overview
6. Set up Docker and `.devcontainer/`
7. Create `.claude/` for spec-kit
8. Update this top-level `CLAUDE.md` and `README.md`

### Managing Shared References

Reference implementations in `references/` are accessible to all domains:

**Current:**
- `faer-rs/` - Rust linear algebra library
- `lapack/` - Reference LAPACK (BSD)
- `SLICOT-Reference/` - BSD-3 SLICOT
- `slicot/` - GPL SLICOT (docs/tests only, NOT source)

**To Add:**
- `spral/` - SPRAL sparse library (BSD-3)

To add new reference:
```bash
cd references/
git clone <reference-repo-url>
# Update domain CLAUDE.md files to mention new reference
```

## Migration to Workspace (Phase 2)

When domains stabilize, migrate to Cargo workspace:

### Pre-Migration Checklist
- [ ] Core algorithms in each domain are implemented and tested
- [ ] APIs are relatively stable
- [ ] Identify common code that should be extracted to `rivrs-core/`
- [ ] All implementations have proper attribution

### Migration Steps
1. Create workspace `Cargo.toml` at root
2. Create `rivrs-core/` with shared utilities, traits, errors
3. Rename domains:
   - `control/` → `rivrs-control/`
   - `sparse/` → `rivrs-sparse/`
4. Update each domain's `Cargo.toml` to depend on `rivrs-core`
5. Move domain `Cargo.toml` and `src/` to workspace member structure
6. Test workspace builds: `cargo build --workspace`
7. Update documentation

### Workspace Structure (Future)
```
rivrs-linalg/
├── Cargo.toml              # Workspace manifest
├── rivrs-core/
│   ├── Cargo.toml
│   └── src/                # Shared traits, errors, utilities
├── rivrs-control/
│   ├── Cargo.toml          # depends on rivrs-core
│   ├── CLAUDE.md
│   └── src/
├── rivrs-sparse/
│   ├── Cargo.toml          # depends on rivrs-core
│   ├── CLAUDE.md
│   └── src/
└── references/
```

## Git Workflow

### Commit Practices (Repository-Wide)
- Commit frequently after logical units
- Commit after tests pass or milestones are verified
- Use clear, descriptive commit messages
- Include Co-Authored-By for AI assistance when appropriate

### Branch Strategy (Current: Simple)
- `main` branch for all development (Phase 1)
- Feature branches for major changes
- No complex branching needed yet (single developer, early stage)

### Branch Strategy (Future: Workspace)
- `main` branch for stable workspace
- Domain-specific branches: `control/*`, `sparse/*`
- Feature branches: `feature/algorithm-name`

## Testing and CI/CD

### Current Testing
Each domain runs its own tests:
```bash
cd control/ && cargo test
cd sparse/ && cargo test
```

### Future CI/CD (Workspace)
```yaml
# GitHub Actions workflow
- Run tests for all workspace members
- Run benchmarks
- Check formatting (rustfmt)
- Check lints (clippy)
- Verify clean room compliance (check for cited sources)
```

## Documentation Standards (Repository-Wide)

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

## Repository Maintenance

### Regular Tasks
- Update reference implementations periodically
- Keep Rust toolchain up to date (check MSRV)
- Update dependencies (cargo update)
- Review and update documentation
- Verify clean room compliance

### Adding New Contributors
New contributors must:
1. Read this CLAUDE.md
2. Read domain-specific CLAUDE.md for their work area
3. Read NOTICE file to understand attribution requirements
4. Understand clean room methodology
5. Know which sources they can/cannot consult

## License and Attribution

**License**: Apache-2.0 (see [LICENSE](LICENSE))

**Attribution**: See NOTICE in each project for complete attribution to:
- Academic sources (papers, textbooks)
- LAPACK (BSD-3-Clause)
- SLICOT-Reference (BSD-3-Clause)
- SPRAL (BSD-3-Clause, planned)
- MatrixEquations.jl (MIT, for validation)

Each domain's implementation cites specific sources used.

## Questions or Issues?

For domain-specific questions, consult the appropriate CLAUDE.md:
- Control systems: [control/CLAUDE.md](control/CLAUDE.md)

For repository-wide questions or structural issues, refer to this file.
