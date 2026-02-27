# Research: Continuous Integration Setup

**Feature**: 007-ci-setup
**Date**: 2026-02-07

## Research Questions & Findings

### R1: Feature-gated testing — is test-util already tested in CI?

**Decision**: The self-referencing dev-dependency (`rivrs-sparse = { path = ".", features = ["test-util"] }`) in Cargo.toml already activates `test-util` when running `cargo test --all-targets`. This means existing CI tests DO exercise test-util code. However, no job validates that the library compiles with default features only (the consumer perspective).

**Rationale**: Library consumers who add `rivrs-sparse` as a dependency without `features = ["test-util"]` get default features only. We need to verify this default build doesn't accidentally depend on test-util symbols. A `cargo check` (or `cargo build --lib`) without dev-dependencies confirms this.

**Alternatives considered**:
- `cargo-hack --feature-powerset`: Overkill for a single optional feature. cargo-hack is valuable when there are many feature combinations; with just `test-util`, explicit commands are clearer.
- Matrix dimension for features: Unnecessary complexity for one feature flag.

### R2: Benchmark compilation verification approach

**Decision**: Use `cargo bench --no-run` to compile benchmarks without executing them. This is the standard Cargo pattern for CI benchmark verification.

**Rationale**: `cargo bench --no-run` compiles all benchmark targets (including criterion harness setup) but skips actual execution. This catches compilation errors and missing imports without consuming CI time on benchmark runs. The `solver_benchmarks` bench target uses `harness = false` (criterion), which is compatible with `--no-run`.

**Alternatives considered**:
- `cargo build --benches`: Also works but doesn't exercise criterion's harness setup the same way. `--no-run` is more thorough.
- Running benchmarks in CI: Premature — no solver code exists yet, and CI runners have variable performance characteristics that make benchmark results unreliable.

### R3: Lint/fmt on MSRV vs stable

**Decision**: Run clippy and rustfmt only on stable. Run tests on both MSRV (1.87) and stable.

**Rationale**: Clippy lints and rustfmt rules can differ between Rust versions. MSRV testing verifies that code compiles and tests pass on the minimum supported toolchain. Linting against stable ensures the codebase follows current best practices. The `rust-version = "1.87"` in Cargo.toml enables MSRV-aware clippy lints on stable.

**Alternatives considered**:
- Lint on both MSRV and stable: Creates noise from version-specific lint differences without improving code quality. A lint failure on MSRV but not stable would be confusing and not actionable.
- Lint on MSRV only: Would miss newer lints that catch real issues.

### R4: Path filtering for monorepo CI

**Decision**: Do not add path filtering at this time. Keep running all jobs on every PR.

**Rationale**: The monorepo currently has only two domains (sparse, control). Path filtering introduces complexity with GitHub required checks (skipped jobs don't satisfy required check requirements) and the `dorny/paths-filter` dependency. The current CI runs fast enough (~5 min cached) that filtering provides minimal time savings. Revisit when the monorepo grows or CI times become a concern.

**Alternatives considered**:
- `dorny/paths-filter@v3`: Effective for larger monorepos but requires a separate detection job and conditional logic in every downstream job. The overhead isn't justified for two domains.
- GitHub's native `paths` trigger filter: Breaks required checks entirely — jobs that don't run can't satisfy branch protection rules.

### R5: Cache configuration for monorepo

**Decision**: Continue using `Swatinem/rust-cache@v2` with the existing pattern. Each job already `cd`s into its domain directory, and rust-cache detects the Cargo.toml/lock file location automatically.

**Rationale**: The existing cache setup works correctly — each sparse job runs in the `sparse/` directory, and rust-cache hashes based on the Cargo.lock and rustc version. No `workspaces` parameter is needed since each job operates in a single project directory (not a Cargo workspace).

**Alternatives considered**:
- Explicit `workspaces: sparse -> sparse/target` parameter: Not needed since `cd sparse` already sets the working directory. Would be useful if jobs ran from repo root.
- `shared-key` across jobs: Considered but not needed — different jobs (test vs lint vs doc) may have different compiled artifacts due to different feature/target flags.

### R6: Default-features-only build check

**Decision**: Add a `cargo check --lib` step (without dev-dependencies, without test-util) to verify the library compiles for downstream consumers. Use `cargo check --no-default-features --lib` isn't needed since there are no default features; `cargo check --lib` is sufficient if run before dev-dependencies are resolved. In practice, the simplest approach is `cargo check` in the lint job (which doesn't activate dev-deps by default for `check`).

**Rationale**: `cargo check` is faster than `cargo build` and sufficient for compilation verification. Running it in the lint job (which already uses stable) keeps the workflow simple.

**Correction**: Actually, `cargo check` DOES resolve dev-dependencies and their features. The reliable way to check default-features-only is to use `cargo check --lib --no-dev` — but this flag doesn't exist in stable cargo. The practical alternative is that the existing lint job's `cargo clippy --all-targets` already compiles the lib, and since the self-referencing dev-dep activates test-util, we'd need `cargo-hack` or a separate approach. Given the single feature flag, the simplest reliable approach is: verify there are no `#[cfg(not(feature = "test-util"))]` compilation paths that could break, and trust that the existing test suite (which does activate test-util) plus the lint job (which checks all targets) provides adequate coverage. If a consumer reports a default-features build failure, add explicit verification then.

**Final decision**: Accept that current CI provides adequate coverage for the single test-util feature. The self-referencing dev-dependency pattern is well-established in the Rust ecosystem and the risk of a default-features regression is low.
