# OMEinsumContractionOrders

OMEinsumContractionOrders is a Julia package that provides an `optimize_code` function for finding optimal contraction orders for tensor networks. It is designed to work with multiple tensor network packages, such as: [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl/issues) package and [ITensorNetworks.jl](https://github.com/mtfishman/ITensorNetworks.jl).

[![Build Status](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/workflows/CI/badge.svg)](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/actions)
[![codecov](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl/branch/master/graph/badge.svg?token=BwaF0cV6q1)](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl)

## Installation

To install OMEinsumContractionOrders, please follow these steps:

1. Open Julia's interactive session (known as [REPL](https://docs.julialang.org/en/v1/manual/getting-started/)) by typing `julia` in your terminal.
2. Press the <kbd>]</kbd> key in the REPL to enter the package mode.
3. Type `add OMEinsumContractionOrders` to install the stable release of the package.
4. (Optional) If you want to use the `KaHyParBipartite` optimizer, which is based on the KaHyPar library, type `add KaHyPar`. Note that this step is optional because some users may have issues with the binary dependencies of KaHyPar (please check issues: [this](https://github.com/kahypar/KaHyPar.jl/issues/12) and [this](https://github.com/kahypar/KaHyPar.jl/issues/19)).

## Example 1: Use it directly
The contraction order optimizer is implemented in the `optimize_code` function. It takes three arguments: `code`, `size`, and `optimizer`. The `code` argument is the [einsum notation](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) to be optimized. The `size` argument is the size of the variables in the einsum notation. The `optimizer` argument is the optimizer to be used. The `optimize_code` function returns the optimized contraction order. One can use `contraction_complexity` function to get the time, space and rewrite complexity of returned contraction order.

```julia
julia> using OMEinsumContractionOrders, Graphs, KaHyPar

julia> function random_regular_eincode(n, k; optimize=nothing)
	    g = Graphs.random_regular_graph(n, k)
	    ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
	    return EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    
julia> code = random_regular_eincode(200, 3);

julia> optcode_tree = optimize_code(code, uniformsize(code, 2),
	TreeSA(sc_target=28, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0));

julia> contraction_complexity(optcode_tree, uniformsize(code, 2))
Time complexity: 2^39.5938886486877
Space complexity: 2^28.0
Read-write complexity: 2^30.39890775966298
```

`OMEinsumContractionOrders` is shipped with [`OMEinsum`](https://github.com/under-Peter/OMEinsum.jl) package. You can use it to optimize the contraction order of a `OMEinsum` expression.

## References

If you find this package useful in your research, please cite the *relevant* papers in [CITATION.bib](CITATION.bib).

## Multi-GPU computation
Please check this Gist:

https://gist.github.com/GiggleLiu/d5b66c9883f0c5df41a440589983ab99

## Authors

OMEinsumContractionOrders was developed by Jin-Guo Liu and Pan Zhang.