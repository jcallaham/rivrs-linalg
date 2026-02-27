# rivrs-linalg — Numerical Linear Algebra for Rivrs

A scientific computing library providing numerical linear algebra implementations
for the [Rivrs](https://github.com/jcallaham/rivrs) symbolic-numeric framework.

**License**: Apache-2.0

## Crates

| Crate | Description | Status |
|-------|-------------|--------|
| [`rivrs-sparse`](sparse/) | Sparse symmetric indefinite direct solver (SSIDS) | Feature-complete, preparing for 0.1.0 release |
| [`rivrs-linalg`](crates/rivrs-linalg/) | Facade re-exporting domain crates as modules | Planned (0.1.0 alongside rivrs-sparse) |
| [`rivrs`](crates/rivrs/) | Umbrella re-exporting rivrs-linalg and future siblings | Planned (0.1.0 alongside rivrs-sparse) |

Users can depend on whichever level they prefer:
- `rivrs-sparse` for just the sparse solver
- `rivrs-linalg` for `rivrs_linalg::sparse`, `rivrs_linalg::control`, etc.
- `rivrs` for `rivrs::linalg`, `rivrs::ode`, `rivrs::optimize`, etc. (future)

## Repository Structure

```
rivrs-linalg/                 # Repository root
├── sparse/                   # rivrs-sparse crate (standalone Rust project)
├── control/                  # rivrs-control crate (standalone, not yet published)
├── crates/
│   ├── rivrs-linalg/         # Facade crate
│   └── rivrs/                # Umbrella crate
├── references/               # Shared reference implementations (LAPACK, SPRAL, faer-rs)
├── NOTICE                    # Attribution and licensing
└── LICENSE                   # Apache-2.0
```

### rivrs-sparse

Feature-complete sparse symmetric indefinite direct solver (SSIDS). Parallel
multifrontal LDL^T with APTP pivoting. 65/65 SuiteSparse matrices passing,
competitive with SPRAL and MUMPS. See [sparse/README.md](sparse/README.md) for
API, benchmarks, and usage.

### rivrs-control (not yet published)

Clean room implementations of control systems algorithms (SLICOT-inspired), using
`faer` for dense linear algebra. Sylvester equation solvers (continuous/discrete)
complete. Lyapunov and Riccati solvers planned.
See [control/README.md](control/README.md) for details.

## Development

Each domain crate is a standalone Rust project:

```bash
cd sparse/   # or control/
cargo build
cargo test
```

See the crate-level CLAUDE.md for detailed development guidance:
- [sparse/CLAUDE.md](sparse/CLAUDE.md)
- [control/CLAUDE.md](control/CLAUDE.md)

## Release Plan

1. **Publish** `rivrs-sparse` 0.1.0 on crates.io (the existing `sparse/` project, as-is).
2. **Publish** `rivrs-linalg` 0.1.0 as a facade that re-exports `rivrs-sparse` as `rivrs_linalg::sparse`.
3. **Publish** `rivrs` 0.1.0 as an umbrella that re-exports `rivrs-linalg` as `rivrs::linalg`.
4. **Add domain crates** over time (`rivrs-control`, then `rivrs-ode`, `rivrs-optimize`, etc.) and wire them into the facades.
5. **Optionally consolidate** into a monolithic `rivrs-linalg` later if the facade indirection becomes undesirable.

## Clean Room Implementation

All implementations are clean room to maintain Apache-2.0 licensing.

- **Permitted**: Academic papers/textbooks, permissively-licensed reference code
  (LAPACK, SLICOT-Reference, SPRAL), documentation and test cases from any library
- **Prohibited**: GPL source code (SLICOT `src/*.f`), proprietary HSL source code

See [NOTICE](NOTICE) for complete attribution.

## License

Apache-2.0 — see [LICENSE](LICENSE) and [NOTICE](NOTICE).
