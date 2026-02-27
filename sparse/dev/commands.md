/speckit.specify implement Phase 9.1a in ssids-plan.md.  Treat any code as notional and descriptive of *what* functionality is generally needed, not *how* to implement it. The actual spec should be based on the current codebase status and the high-level progress through the plan - mark anything that's unclear or conflicting for further clarification. Make sure to note algorithm references and the locations of the corresponding Markdown files in /workspace/rivrs-linalg/references/ssids/*.md

Can you review PR #22, looking for code quality, consistency with the plan outlined in docs/ssids-plan.md, Rust best practices, duplication of code or functionality, and code clarity? Also be sure to evaluate test comprehensiveness following the core principle that correctness is the critical goal of the implementation - do we have sufficient test coverage to have high confidence in the current implementation?


/speckit.specify Can you look at items 1-3 in docs/roadmap.md execute the refactor, including related edits to documentation, tests, and examples?  Feel free to suggest changes to the proposed structure as you see fit - it doesn't have to be perfect at this stage but the goal is to reorganize to optimize for future reusability, modularity, and composability as this projects moves towards a full-fledged sparse linear algebra library.

%%

- Review codebase & structure
- Documentation - how do we build and view, and what's the appropriate scope?  All public exports?
