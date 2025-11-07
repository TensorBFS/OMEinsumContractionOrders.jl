# Benchmark Configuration Details

## Hardware Specifications

- **CPU**: Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz
- **Threading**: Single-threaded execution
- **Test Instance**: sycamore_53_20_0

## Overview

This document provides detailed parameter configurations for benchmarking tensor network contraction order optimizers. The benchmarks compare various optimizer implementations with different parameter settings to evaluate the trade-off between optimization time and solution quality.

## Parameter Configuration

### Julia Optimizers (OMEinsumContractionOrders.jl)

#### Treewidth

| Parameter | Values |
|-----------|--------|
| algorithm | MF, MMD, AMF |

#### KaHyParBipartite

| Parameter | Values |
|-----------|--------|
| sc_target | 25 |
| imbalances | 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 |

#### TreeSA

| Parameter | Values |
|-----------|--------|
| niters | 1, 2, 4, 6, 8, 10, 20, 30, 40, 50 |
| score | TC (tc_weight=1), SC (sc_weight=1) |

#### GreedyMethod

| Parameter | Values |
|-----------|--------|
| α (temperature) | 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 |

#### HyperND

| Parameter | Values |
|-----------|--------|
| variant | base, METISND, KaHyParND |
| imbalances | 100, 110, 120, ..., 800 (METISND/KaHyParND only) |
| score | TC (tc_weight=1), SC (sc_weight=1) |

#### SABipartite

| Parameter | Values |
|-----------|--------|
| imbalances | (parameter scan similar to KaHyParBipartite) |

### Python Optimizers (cotengra)

#### cotengra_greedy

| Parameter | Values |
|-----------|--------|
| max_repeats | 1, 5, 10, 20, 50 |
| minimize | flops, size |

#### cotengra_kahypar

| Parameter | Values |
|-----------|--------|
| imbalance | 0.01, 0.1, 0.3, 0.5, 0.8 |
| minimize | flops, size |

## Benchmark Metrics

For each optimizer configuration, the following metrics are measured:

1. **Contraction Cost**: The computational complexity of the resulting contraction order (typically measured as log2 of FLOPs)
2. **Optimization Time**: Wall-clock time required to find the contraction order
3. **Space Complexity**: Peak memory requirement during contraction
4. **Read-Write Complexity**: Total data movement cost

## Results Summary

The benchmark results are visualized in Figure 1 of the main paper (`paper.md`), showing the Pareto front of multi-objective optimization balancing contraction order quality against optimization runtime.

### Key Findings

- **TreeSA** and **HyperND**: Achieve the lowest contraction costs but require longer optimization times
- **GreedyMethod**: Offers the fastest optimization but with higher contraction costs
- **KaHyParBipartite**: Provides a good balance for large-scale problems
- Performance is highly problem-dependent, with no single optimizer dominating across all metrics

## Reproducing the Benchmarks

The benchmark scripts and data are available in the [OMEinsumContractionOrdersBenchmark](https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark) repository.

To reproduce the benchmarks:

```bash
# Clone the benchmark repository
git clone https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark.git
cd OMEinsumContractionOrdersBenchmark

# Follow the instructions in the repository README
```

## Citation

If you use these benchmarks in your research, please cite:

```bibtex
@article{OMEinsumContractionOrders,
  title={OMEinsumContractionOrders: A Julia package for tensor network contraction order optimization},
  author={Liu, Jin-Guo and Gao, Xuan-Zhao and Samuelson, Richard},
  journal={Journal of Open Source Software},
  year={2025}
}
```

