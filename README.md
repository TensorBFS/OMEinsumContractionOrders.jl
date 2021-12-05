# OMEinsumContractionOrders
This package provides `optimize_code` function for finding tensor network contraction orders.

[![Build Status](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/workflows/CI/badge.svg)](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/actions)

## Installation
<p>
OMEinsumContractionOrders is a Julia Language package. To install OMEinsumContractionOrders,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>

For stable release

```julia
pkg> add OMEinsumContractionOrders
pkg> add KaHyPar
```
The `KaHyPar` package (used in `KaHyParBipartite` optimizer) is optional because some one may have the binary issue, check [this](https://github.com/kahypar/KaHyPar.jl/issues/12) and [this](https://github.com/kahypar/KaHyPar.jl/issues/19).

## Example
Contract a tensor network
```julia
julia> using OMEinsum, OMEinsumContractionOrders, Graphs, KaHyPar

julia> function random_regular_eincode(n, k; optimize=nothing)
	    g = Graphs.random_regular_graph(n, k)
	    ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
	    return EinCode((ixs..., [(i,) for i in     Graphs.vertices(g)]...), ())
    end
    
julia> code = random_regular_eincode(200, 3);

julia> optcode_tree = optimize_code(code, uniformsize(code, 2), TreeSA(sc_target=28, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0));

julia> optcode_kahypar = optimize_code(code, uniformsize(code, 2), KaHyParBipartite(sc_target=30, max_group_size=50));

julia> optcode_sa = optimize_code(code, uniformsize(code, 2), SABipartite(sc_target=30, max_group_size=50));

julia> OMEinsum.timespace_complexity(code, uniformsize(code, 2))
(200.0, 0.0)

julia> OMEinsum.timespace_complexity(optcode_kahypar, uniformsize(code, 2))
(38.0290167456887, 26.0)

julia> OMEinsum.timespace_complexity(optcode_sa, uniformsize(code, 2))
(34.86528023060411, 27.0)

julia> tc, sc = OMEinsum.timespace_complexity(optcode_tree, uniformsize(code, 2))
(30.541894421918297, 26.0)
```

## References

If you find this package useful in your research, please cite the following papers

To credit the `KaHyParBipartite` and `SABipartite` method,
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

To credit the `KaHyParBipartite` method,
```
@Article{10.21468/SciPostPhys.7.5.060,
	title={{Fast counting with tensor networks}},
	author={Stefanos Kourtis and Claudio Chamon and Eduardo R. Mucciolo and Andrei E. Ruckenstein},
	journal={SciPost Phys.},
	volume={7},
	issue={5},
	pages={60},
	year={2019},
	publisher={SciPost},
	doi={10.21468/SciPostPhys.7.5.060},
	url={https://scipost.org/10.21468/SciPostPhys.7.5.060},
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

To credit the `TreeSA` method,
```
@misc{kalachev2021recursive,
      title={Recursive Multi-Tensor Contraction for XEB Verification of Quantum Circuits}, 
      author={Gleb Kalachev and Pavel Panteleev and Man-Hong Yung},
      year={2021},
      eprint={2108.05665},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

## Authors
Jin-Guo Liu and Pan Zhang
