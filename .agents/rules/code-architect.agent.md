# Rule: Architecture, Layering & Modularity

Guidelines for evaluating and structuring Python package dependencies and module organization.

## Scope & Token Efficiency
- **Targeted Scope:** Focus architectural analysis strictly on the specific package, submodule, or design question requested. Do not attempt whole-repository full-text scans.
- **Context Guard:** Never read more than 5 full files into context. Use symbol navigation, AST summaries, or grep to evaluate dependencies efficiently.

## Universal 3-Tier Layering Principles
Organize Python package code into a clean bottom-up dependency hierarchy:
1. **Tier 1: Primitives & Foundations** (system utilities, OS wrappers, base constants, raw data containers). Must NOT import from higher tiers.
2. **Tier 2: Core Domain & Engines** (algorithms, numerical routines, domain data models, transformation engines). Imports Tier 1 only.
3. **Tier 3: Interfaces & Pipelines** (public exports, high-level orchestration, workflow runners, CLI commands). Imports Tier 1 and Tier 2.

## Dependency Rules
- **Strict Directionality:** Higher tiers may import lower tiers; lower tiers must NEVER import higher tiers.
- **Zero Circular Dependencies:** Circular imports are strictly forbidden. Extract shared dependencies into lower-tier utility modules or use deferred typing imports (`if TYPE_CHECKING:`).
- **Loose Coupling:** Minimize inter-module cross-calls. Keep modules focused on a single responsibility.