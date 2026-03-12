# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OMEinsumContractionOrders.jl is a Julia package for finding optimal contraction orders for tensor networks (Einstein summation expressions). It provides multiple optimization algorithms and is used by OMEinsum.jl, Yao.jl, GenericTensorNetworks.jl, TensorInference.jl, and ITensorNetworks.jl.

## Common Commands

```bash
make init          # Initialize project dependencies
make test          # Run full test suite
make coverage      # Run tests with coverage
make serve         # Live-serve documentation locally
make update        # Update all dependencies
```

Run a single test file:
```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["greedy"])'
```

Run tests interactively:
```bash
julia --project test/runtests.jl
```

## Architecture

### Core Type Hierarchy

**Einsum representations** (`src/Core.jl`):
- `EinCode{LT}` — flat einsum with input index vectors (`ixs`) and output indices (`iy`)
- `NestedEinsum{LT}` — binary tree of contractions (the optimized output)
- `SlicedEinsum{LT,ET}` — wraps a NestedEinsum with sliced indices for memory reduction

**Optimizers** (each `<: CodeOptimizer`):
| Type | File | Algorithm |
|------|------|-----------|
| `GreedyMethod` | `greedy.jl` | Fast greedy pair selection using `IncidenceList` graph |
| `TreeSA` | `treesa.jl` | Simulated annealing over `ExprTree` structures |
| `Treewidth` / `ExactTreewidth` | `treewidth.jl` | Tree decomposition via CliqueTrees.jl |
| `SABipartite` | `sabipartite.jl` | SA-based recursive bipartitioning |
| `KaHyParBipartite` | `kahypar.jl` | KaHyPar graph partitioning (extension) |
| `HyperND` | `hypernd.jl` | Nested dissection using CliqueTrees |

**Slicers** (`<: CodeSlicer`): `TreeSASlicer` in `treesaslicer.jl`

**Simplifiers** (`<: CodeSimplifier`): `MergeVectors`, `MergeGreedy` in `simplify.jl`

**Complexity** (`complexity.jl`): `ContractionComplexity` with `tc` (time), `sc` (space), `rwc` (read-write) fields. Scoring via `ScoreFunction` with configurable weights.

### Data Flow

1. User calls `optimize_code(eincode, size_dict, optimizer)` (`interfaces.jl`)
2. Dispatcher converts `EinCode` → algorithm-specific internal representation
3. Algorithm optimizes contraction order using `size_dict` for complexity estimation
4. Returns `NestedEinsum` (binary contraction tree)
5. Optional: `slice_code` wraps result in `SlicedEinsum`; `optimize_permute` optimizes index ordering

### Key Internal Structures

- `IncidenceList` (`incidencelist.jl`) — sparse hypergraph for greedy algorithm; maps vertices to edges and tracks edge sizes
- `ExprTree` (`treesa.jl`) — mutable binary tree for simulated annealing; supports random initialization and local mutations

### Extensions

- `ext/KaHyParExt.jl` — enables `KaHyParBipartite` optimizer when KaHyPar.jl is loaded
- `ext/LuxorTensorPlot/` — provides `viz_eins()` and `viz_contraction()` when LuxorGraphPlot.jl is loaded

### Serialization

`json.jl` handles round-trip JSON serialization of contraction orders (`NestedEinsum` ↔ JSON).

## Julia Version

Requires Julia 1.8+. Current package version: 1.2.3.
