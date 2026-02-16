# Specification Quality Checklist: Triangular Solve & Solver API

**Purpose**: Validate specification completeness and quality before proceeding to planning
**Created**: 2026-02-15
**Feature**: [spec.md](../spec.md)

## Content Quality

- [x] No implementation details (languages, frameworks, APIs)
- [x] Focused on user value and business needs
- [x] Written for non-technical stakeholders
- [x] All mandatory sections completed

## Requirement Completeness

- [x] No [NEEDS CLARIFICATION] markers remain — both resolved via SPRAL reference analysis
- [x] Requirements are testable and unambiguous
- [x] Success criteria are measurable
- [x] Success criteria are technology-agnostic (no implementation details)
- [x] All acceptance scenarios are defined
- [x] Edge cases are identified
- [x] Scope is clearly bounded
- [x] Dependencies and assumptions identified

## Feature Readiness

- [x] All functional requirements have clear acceptance criteria
- [x] User scenarios cover primary flows
- [x] Feature meets measurable outcomes defined in Success Criteria
- [x] No implementation details leak into specification

## Notes

- Both design questions resolved by studying SPRAL's implementation (BSD-3):
  1. **Scaling integration**: Pass scaling factors into factorization, apply at assembly/scatter layer (SPRAL's `add_a_block` pattern). No pre-scaled matrix copy. APTP kernel remains scaling-unaware.
  2. **Rank-deficiency handling**: Follow SPRAL's `action=true` default — zero pivots produce zeroed solution components, detectable via `stats().zero_pivots` and `inertia()`. No error by default.
- The spec intentionally references algorithm references and codebase locations in the "Algorithm References" and "Open Design Questions" sections — these are informational for the planning phase, not implementation details.
- Multi-column RHS and parallelism are explicitly deferred (documented in Assumptions).
