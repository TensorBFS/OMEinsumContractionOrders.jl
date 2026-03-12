# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository. For full project context, see [CLAUDE.md](./CLAUDE.md).

## Common Commands

```bash
make init          # Initialize project dependencies
make test          # Run full test suite
julia --project -e 'using Pkg; Pkg.test(test_args=["greedy"])'  # Run single test
```

## Architecture Summary

Tensor network contraction order optimizer. Core flow: `EinCode` → `optimize_code()` → `NestedEinsum` (binary contraction tree).

**Key types:** `EinCode`, `NestedEinsum`, `SlicedEinsum` (in `src/Core.jl`). Optimizers (`<: CodeOptimizer`): `GreedyMethod`, `TreeSA`, `Treewidth`, `SABipartite`, `KaHyParBipartite`, `HyperND`. Entry point: `optimize_code()` in `src/interfaces.jl`.

See CLAUDE.md for detailed type hierarchy, data flow, and file-by-file descriptions.

## Coding Conventions

- Julia 1.8+ required
- Tests use `@testset` blocks in `test/runtests.jl` with per-algorithm test files
- Package extensions in `ext/` for optional dependencies (KaHyPar, LuxorGraphPlot)
- Performance-sensitive code — avoid unnecessary allocations (see recent commits for allocation-free patterns in `greedy.jl`)
