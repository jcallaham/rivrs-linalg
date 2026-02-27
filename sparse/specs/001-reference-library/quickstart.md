# Quickstart: Phase 0.1 Reference Library

## What Was Built

Four documentation deliverables completing the SSIDS Phase 0.1 literature review:

| Document | Location | Purpose |
|----------|----------|---------|
| Paper Index | `references/ssids/INDEX.md` | Find any paper by category, citation, or topic |
| SPRAL Code Review | `sparse/dev/references/SPRAL-CODE-REVIEW.md` | Understand SPRAL SSIDS architecture without reading Fortran |
| faer Integration Notes | `sparse/dev/references/FAER-INTEGRATION-NOTES.md` | Know which faer components to reuse vs build new |
| APTP Algorithm | `sparse/dev/references/APTP-ALGORITHM.md` | Implementation-ready pseudocode with paper citations |

## How to Use

### Starting implementation of a new SSIDS component

1. Read `APTP-ALGORITHM.md` for the algorithm pseudocode
2. Check `FAER-INTEGRATION-NOTES.md` to see if faer already provides infrastructure
3. Consult `SPRAL-CODE-REVIEW.md` to understand how SPRAL implements it
4. Look up the original paper in `INDEX.md` for deeper understanding

### Finding a specific paper

Open `references/ssids/INDEX.md` and search by author, year, or topic category.

### Understanding the SPRAL solver flow

Read `SPRAL-CODE-REVIEW.md` sections 1-3 for the executive summary and data flow diagram.

### Deciding whether to reuse faer code

Check the reuse classification table in `FAER-INTEGRATION-NOTES.md`:
- **Direct use**: Call faer's existing API
- **Adapt**: Use the pattern but modify for APTP
- **Reference only**: Study for understanding, build from scratch

## What's Next

After Phase 0.1, proceed to:
- **Phase 0.2**: Test matrix collection (SuiteSparse, hand-constructed, SPRAL tests)
- **Phase 0.3**: SPRAL golden results generation
- **Phase 0.4**: Final repository setup
