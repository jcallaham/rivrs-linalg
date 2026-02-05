# Feature Specification: Complete Phase 0.1 Reference Library

**Feature Branch**: `001-reference-library`
**Created**: 2026-02-05
**Status**: Draft
**Input**: User description: "Review Phase 0.1 and identify what needs to be done to complete the reference library. Many papers have been compiled in ../references/ssids, and SPRAL and faer source code is also available in ../references/. Let me know if any manual work is necessary to compile other references as well."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Reference All Required Papers (Priority: P1)

A developer starting SSIDS implementation needs to verify that every academic paper listed in Phase 0.1 of the SSIDS plan is present, correctly identified, and organized so they can study the algorithms without hunting for sources.

**Why this priority**: Without the complete paper collection, no informed implementation decisions can be made. This is the foundation for all subsequent work.

**Independent Test**: Can be verified by checking every paper listed in Phase 0.1 of ssids-plan.md against the contents of `references/ssids/` and confirming each has both a readable markdown version and an archival PDF.

**Acceptance Scenarios**:

1. **Given** the Phase 0.1 plan lists specific papers under Core APTP, Multifrontal, Pivoting, and Ordering categories, **When** a developer checks `references/ssids/`, **Then** every listed paper is present with both `.md` and `.pdf` versions.
2. **Given** a paper index file exists, **When** a developer reads it, **Then** each entry maps the filename to its full citation, category, and relevance to the SSIDS project.
3. **Given** the plan lists "RAL Technical Reports on SSIDS version 2.0," **When** a developer checks the references, **Then** either the technical reports are present or the index documents that the peer-reviewed publication (Duff, Hogg, Lopez 2020) supersedes them.

---

### User Story 2 - SPRAL Source Code Review Notes (Priority: P2)

A developer needs annotated notes on SPRAL's SSIDS implementation so they can understand the reference solver's architecture without re-reading thousands of lines of Fortran from scratch.

**Why this priority**: Understanding SPRAL's code structure is essential for designing the Rust port, but the review is interpretive work that builds on the paper collection.

**Independent Test**: Can be verified by confirming that a SPRAL code review document exists covering the key SSIDS source files, their responsibilities, data flow between analysis/factor/solve phases, and key design decisions.

**Acceptance Scenarios**:

1. **Given** SPRAL source is available at `references/spral/src/ssids/`, **When** a developer reads the review document, **Then** they find a module-by-module summary of the SSIDS solver covering at least: core solver entry points, symbolic analysis, numeric factorization, triangular solve, and data structures.
2. **Given** the review document exists, **When** a developer looks for APTP-specific implementation details, **Then** the document identifies which SPRAL files implement APTP logic, how pivoting decisions flow, and how delayed columns are managed.

---

### User Story 3 - faer Integration Points Documentation (Priority: P2)

A developer needs a documented map of faer's sparse infrastructure showing which existing components can be reused directly, which need adaptation, and what must be built from scratch for APTP.

**Why this priority**: Identifying reuse opportunities upfront prevents duplicating existing faer infrastructure and guides architectural decisions.

**Independent Test**: Can be verified by confirming a faer integration notes document exists that maps SSIDS solver needs to specific faer modules with reuse assessments.

**Acceptance Scenarios**:

1. **Given** faer has sparse infrastructure in `references/faer-rs/faer/src/sparse/`, **When** a developer reads the integration notes, **Then** they find a component-by-component mapping: CSC storage, elimination trees, AMD ordering, sparse LDL^T, triangular solve, permutation utilities, and workspace management.
2. **Given** faer's sparse LDL^T is in `cholesky.rs`, **When** the integration notes discuss it, **Then** they clarify whether it handles indefinite systems or only positive definite, and what modifications APTP requires.
3. **Given** the integration notes exist, **When** a developer reviews them, **Then** each faer component has a reuse classification: "direct use," "adapt," or "build new."

---

### User Story 4 - Algorithm Pseudocode Extraction (Priority: P3)

A developer needs the core APTP algorithm pseudocode extracted from the academic papers into a standalone, implementation-ready reference document so they can code from it without flipping between multiple papers.

**Why this priority**: Pseudocode extraction is a synthesis step that makes implementation faster, but developers can also work directly from the papers if needed.

**Independent Test**: Can be verified by confirming an algorithm document exists that presents the APTP factorization, pivoting decision logic, and delayed column handling in clear pseudocode with citations back to the source papers.

**Acceptance Scenarios**:

1. **Given** the APTP algorithm is described across duff2020.md and hogg2016.md, **When** a developer reads the algorithm document, **Then** they find a unified pseudocode for the simplicial APTP factorization including: column processing loop, 1x1 pivot acceptance test, 2x2 pivot fallback (Bunch-Kaufman), and column delay mechanism.
2. **Given** the algorithm document exists, **When** a developer checks the multifrontal and symbolic analysis sections, **Then** they find pseudocode for elimination tree traversal, symbolic factorization (column counts and nonzero structure prediction), and contribution assembly.

---

### Edge Cases

- What happens when a listed paper's markdown conversion is incomplete or garbled? The index should flag quality issues for any papers that need re-conversion.
- What happens when a faer component's API has changed since the reference snapshot? Integration notes should record the faer version reviewed and flag any version-sensitive findings.
- What happens when SPRAL source files reference internal HSL routines? The review document should note these dependencies and identify BSD-licensed alternatives.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Repository MUST contain every academic paper listed in Phase 0.1 of ssids-plan.md in both markdown and PDF format under `references/ssids/`.
- **FR-002**: Repository MUST contain a paper index file (`references/ssids/INDEX.md`) mapping each file to its full bibliographic citation, category (Core APTP, Multifrontal, Pivoting, Ordering), and relevance summary.
- **FR-003**: Index MUST document the status of "RAL Technical Reports on SSIDS version 2.0" — either including them or explaining that the Duff/Hogg/Lopez 2020 journal publication supersedes them.
- **FR-004**: Repository MUST contain a SPRAL code review document (`docs/references/SPRAL-CODE-REVIEW.md`) covering SSIDS module structure, data flow, APTP implementation details, and key design decisions.
- **FR-005**: Repository MUST contain a faer integration notes document (`docs/references/FAER-INTEGRATION-NOTES.md`) mapping each SSIDS solver need to a specific faer component with a reuse classification.
- **FR-006**: Repository MUST contain an algorithm pseudocode document (`docs/references/APTP-ALGORITHM.md`) with implementation-ready pseudocode for APTP factorization, symbolic analysis, and pivoting logic, citing source papers for each section.
- **FR-007**: All documents MUST cite their sources (paper references, specific SPRAL files consulted, faer module paths) to maintain clean room audit trail.

### Key Entities

- **Paper**: An academic publication with citation metadata, categorization, and a relevance summary for the SSIDS project.
- **Reference Component**: A module or function in SPRAL or faer that maps to a capability needed by the SSIDS solver, classified by reuse potential.
- **Algorithm Pseudocode**: A language-agnostic step-by-step procedure extracted from one or more papers, annotated with source citations.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of papers listed in Phase 0.1 of ssids-plan.md are present with both markdown and PDF versions, or the index documents why a specific item is superseded/unavailable.
- **SC-002**: A developer unfamiliar with SPRAL can read the code review document and correctly describe the SSIDS solver's three-phase architecture (analyze, factor, solve) and the role of APTP within 30 minutes.
- **SC-003**: The faer integration notes identify at least 5 directly reusable components with specific module paths and version information.
- **SC-004**: The APTP algorithm pseudocode covers the complete simplicial factorization loop including pivot acceptance, 2x2 fallback, and column delay — sufficient that a developer could implement from pseudocode alone without re-reading the papers.
- **SC-005**: Every document in the reference library includes source citations traceable to specific papers, SPRAL files, or faer modules.

## Assumptions

- The faer reference snapshot in `references/faer-rs/` is recent enough to reflect the current sparse infrastructure. Integration notes will record the version reviewed.
- The SPRAL reference in `references/spral/` contains the complete SSIDS source. If files are missing, the review document will note gaps.
- The Duff/Hogg/Lopez 2020 SIAM publication is assumed to supersede the RAL Technical Reports on SSIDS v2.0, as it is the peer-reviewed version of the same research. This will be noted in the index rather than requiring the user to obtain the technical reports separately.
- Markdown conversions of papers are "good enough" for reference use. If any are garbled or incomplete, the index will flag them for re-conversion but the PDF serves as the authoritative source.

## Current State Assessment

### Papers Already Compiled (14 of ~14 required)

| Plan Category | Paper | File | Status |
|---|---|---|---|
| Core APTP | Hogg, Duff, Lopez (2020) | duff2020.md/.pdf | Present |
| Core APTP | Hogg, Ovtchinnikov, Scott (2016) | hogg2016.md/.pdf | Present |
| Core APTP | RAL Technical Reports on SSIDS v2.0 | — | Superseded by duff2020 |
| Multifrontal | Duff & Reid (1983/1984) | duff1984.md/.pdf | Present |
| Multifrontal | Liu (1992) | liu1992.md/.pdf | Present |
| Multifrontal | Davis survey | davis2016.md/.pdf | Present |
| Pivoting | Schenk & Gartner (2006) | schenk2006.md/.pdf | Present |
| Pivoting | Duff & Pralet (2005) | duff2005.md/.pdf | Present |
| Ordering | George & Liu - nested dissection | george1973.md/.pdf | Present |
| Ordering | George & Liu - minimum degree | george1989.md/.pdf | Present |
| Ordering | Gilbert, Ng, Peyton - symbolic factorization | gilbert1992.md/.pdf, gilbert1994.md/.pdf | Present |
| Ordering | Duff & Koster - matching orderings | duff1999.md/.pdf, duff2001.md/.pdf | Present |
| Supernodal | Ng & Peyton (1993) | ng1993.md/.pdf | Present |
| faer | faer JOSS paper | references/faer-rs/paper.md | Present in faer repo |

### Remaining Deliverables

| Deliverable | Status | Effort |
|---|---|---|
| Paper index (INDEX.md) | Not started | Low — cataloging existing files |
| SPRAL code review document | Not started | Medium — requires reading SSIDS Fortran source |
| faer integration notes | Not started | Medium — requires reviewing faer sparse modules |
| APTP algorithm pseudocode | Not started | Medium — synthesis from duff2020 and hogg2016 |

### Manual Work Required

No additional papers need to be manually obtained. All papers listed in Phase 0.1 are already compiled. The remaining work is synthesis and documentation that can be performed by reviewing the existing references.
