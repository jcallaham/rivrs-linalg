# Sparse Symmetric Indefinite Direct Solvers (SSIDS)

This directory contains sparse symmetric indefinite direct solver implementations for rivrs-linalg.

**Status**: Phase 1 - Early scaffolding (no algorithms implemented yet)

## Overview

Clean room implementation of sparse symmetric indefinite direct solvers based on SPRAL (Sparse Parallel Robust Algorithms Library). Implements multifrontal LDLT factorization for solving sparse symmetric indefinite linear systems Ax = b.

**Primary Reference**: SPRAL (BSD-3-Clause) - https://github.com/ralgr/spral

## Planned Algorithms

### Multifrontal LDLT Factorization

Solve sparse symmetric indefinite systems:
```
A x = b
```

**Components:**
1. **Symbolic Analysis**: Matrix ordering (AMD, nested dissection), elimination tree
2. **Numeric Factorization**: Multifrontal LDLT with threshold partial pivoting
3. **Solve Phase**: Forward/backward substitution using factorization

**Related Solvers**: HSL MA27, MA57, MA97 (proprietary) - this is a permissively-licensed alternative

## Current Status

⚠️ **No algorithms implemented yet.** This is a scaffold for future development.

**Completed:**
- ✅ Project structure
- ✅ Build configuration
- ✅ Error type definitions

**Next Steps:**
1. Add SPRAL to parent `references/` directory
2. Implement sparse matrix validation
3. Build elimination tree construction
4. Implement AMD ordering
5. Symbolic factorization
6. Numeric LDLT factorization

## Quick Start

```bash
cd ssids
cargo build
cargo test
```

Currently only has placeholder tests to verify the build works.

## Development

See [CLAUDE.md](CLAUDE.md) for detailed development guidance, including:
- Clean room implementation methodology (SPRAL is permissively licensed, HSL is not)
- Algorithm references and academic papers
- Implementation roadmap
- Testing strategy

### Docker Development Environment

(To be created - will provide isolated environment with SPRAL reference materials)

## Academic References

### Core Algorithm
- Duff, I.S. & Reid, J.K. (1983). "The Multifrontal Solution of Indefinite Sparse
  Symmetric Linear Systems". *ACM TOMS*, 9(3):302-325.

### Parallel Algorithms
- Hogg, J.D., Reid, J.K. & Scott, J.A. (2010). "Design of a Multicore Sparse
  Cholesky Factorization Using DAGs". *SIAM SISC*, 32(6):3627-3649.

### Ordering
- Amestoy, P.R., Davis, T.A. & Duff, I.S. (2004). "Algorithm 837: AMD, an
  Approximate Minimum Degree Ordering Algorithm". *ACM TOMS*, 30(3):381-388.

See [CLAUDE.md](CLAUDE.md) for complete reference list.

## License

Part of rivrs-linalg, licensed under Apache-2.0. See [../LICENSE](../LICENSE).

Based on algorithms from SPRAL (BSD-3-Clause) and academic publications.
