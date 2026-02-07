# Data Model: Continuous Integration Setup

**Feature**: 007-ci-setup
**Date**: 2026-02-07

## Overview

This feature modifies CI configuration (YAML) rather than application code. There is no traditional data model with entities, fields, or relationships. The "data" is the workflow configuration itself.

## Configuration Entities

### CI Workflow

The single workflow file `.github/workflows/ci.yml` defines all CI behavior.

- **Trigger events**: push to main, pull requests targeting main
- **Environment variables**: CARGO_TERM_COLOR, RUST_BACKTRACE (workflow-level)
- **Jobs**: Independent units of work, each with its own runner and steps

### Jobs (sparse domain)

| Job Name | Purpose | Toolchain | Key Steps |
|----------|---------|-----------|-----------|
| test-sparse | Correctness validation | MSRV (1.87) + stable | `cargo test --all-targets` |
| lint-sparse | Code quality | stable | `cargo fmt --check`, `cargo clippy` |
| doc-sparse | Documentation | stable | `cargo doc --no-deps` |
| bench-sparse | Benchmark compilation | stable | `cargo bench --no-run` |

### Toolchain Matrix

| Dimension | Values | Rationale |
|-----------|--------|-----------|
| Rust version | `[stable, "1.87"]` | MSRV + current stable |
| OS | `[ubuntu-latest]` | Linux-only until platform-specific code exists |

## State Transitions

Not applicable — CI jobs are stateless. Each run starts from a clean checkout with only build caches carried over.

## Validation Rules

- All test jobs must pass before merge (enforced by GitHub branch protection, not CI config itself)
- Clippy warnings are treated as errors (`-D warnings`)
- Rustdoc warnings are treated as errors (`RUSTDOCFLAGS: -D warnings`)
- No retry logic — failures indicate real problems (Constitution I)
