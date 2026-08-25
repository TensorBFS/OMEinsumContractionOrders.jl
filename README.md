# OMEinsumContractionOrders

OMEinsumContractionOrders is a Julia package that provides an `optimize_code` function for finding optimal contraction orders for tensor networks. It is designed to work with multiple tensor network packages, such as: [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl/issues) package and [ITensorNetworks.jl](https://github.com/mtfishman/ITensorNetworks.jl). A high-performance Rust port with Python bindings is available as [OMECO](https://github.com/GiggleLiu/omeco).

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://tensorbfs.github.io/OMEinsumContractionOrders.jl/dev/)
[![Build Status](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/workflows/CI/badge.svg)](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/actions)
[![codecov](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl/branch/master/graph/badge.svg?token=BwaF0cV6q1)](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.09886/status.svg)](https://doi.org/10.21105/joss.09886)

## Performance

On same-server benchmarks with matched 90 s budgets, contraction orders found by this package beat [cotengra](https://github.com/jcmgray/cotengra)'s hyper-optimizer on every tested instance — e.g. `tc` 39.22 vs 69.04 on the king graph (≈10⁹× fewer flops), 28.03 vs 41.15 on `dbn_13`, and 48.98 vs 54.03 on surface code `d=21` (lower `tc` = fewer flops). The Rust port [OMECO](https://github.com/GiggleLiu/omeco) pushes these records further.

## Installation

To install OMEinsumContractionOrders, please follow these steps:

1. Open Julia's interactive session (known as [REPL](https://docs.julialang.org/en/v1/manual/getting-started/)) by typing `julia` in your terminal.
2. Press the <kbd>]</kbd> key in the REPL to enter the package mode.
3. Type `add OMEinsumContractionOrders` to install the stable release of the package.

For using, please check the [documentation](https://tensorbfs.github.io/OMEinsumContractionOrders.jl/dev/).

## Contributing

For a contributor guide, please refer to [CLAUDE.md](CLAUDE.md) for details.

## References

If you find this package useful in your research, please cite the *relevant* papers in [CITATION.bib](CITATION.bib).