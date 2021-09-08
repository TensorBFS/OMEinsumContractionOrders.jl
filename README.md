# OMEinsumContractionOrders
This package provides function `optimize_kahypar` for finding tensor network contraction orders. So far, it is only tested on Ubuntu. Mac OS and Windows users have to wait for the [Binary issue](https://github.com/kahypar/KaHyPar.jl/issues/12) of KaHyPar being solved.

[![Build Status](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/workflows/CI/badge.svg)](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/actions)

## Installation
<p>
OMEinsumContractionOrders is a Julia Language package. To install OMEinsumContractionOrders,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>

For stable release

```julia
pkg> add https://github.com/kahypar/KaHyPar.jl.git#master
pkg> add OMEinsumContractionOrders
```

## Example
Contract a tensor network
```julia
julia> using OMEinsum, OMEinsumContractionOrders, LightGraphs

julia> function random_regular_eincode(n, k; optimize=nothing)
	    g = LightGraphs.random_regular_graph(n, k)
	    ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
	    return EinCode((ixs..., [(i,) for i in     LightGraphs.vertices(g)]...), ())
    end
    
julia> code = random_regular_eincode(200, 3);

julia> optcode = optimize_kahypar(code, uniformsize(code, 2); sc_target=30, max_group_size=50);

julia> optcode_sa = optimize_sa(code, uniformsize(code, 2); sc_target=30, max_group_size=50);

julia> OMEinsum.timespace_complexity(code, uniformsize(code, 2))
(200.0, 0.0)

julia> OMEinsum.timespace_complexity(optcode, uniformsize(code, 2))
(38.0290167456887, 26.0)

julia> OMEinsum.timespace_complexity(optcode_sa, uniformsize(code, 2))
(34.86528023060411, 27.0)
```

## References

If you find this package useful in your research, please cite the following papers

```
@misc{Pan2021,
      title={Simulating the Sycamore quantum supremacy circuits}, 
      author={Feng Pan and Pan Zhang},
      year={2021},
      eprint={2103.03074},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

```
@article{Gray2021,
   title={Hyper-optimized tensor network contraction},
   volume={5},
   ISSN={2521-327X},
   url={http://dx.doi.org/10.22331/q-2021-03-15-410},
   DOI={10.22331/q-2021-03-15-410},
   journal={Quantum},
   publisher={Verein zur Forderung des Open Access Publizierens in den Quantenwissenschaften},
   author={Gray, Johnnie and Kourtis, Stefanos},
   year={2021},
   month={Mar},
   pages={410}
}
```

## Authors
Jin-Guo Liu and Pan Zhang
