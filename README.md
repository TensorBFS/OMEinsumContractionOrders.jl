# OMEinsumContractionOrders
This package provides `optimize_code` function for finding tensor network contraction orders.

[![Build Status](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/workflows/CI/badge.svg)](https://github.com/TensorBFS/OMEinsumContractionOrders.jl/actions)
[![codecov](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl/branch/master/graph/badge.svg?token=BwaF0cV6q1)](https://codecov.io/gh/TensorBFS/OMEinsumContractionOrders.jl)

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
Optimize a contraction order
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

julia> optcode_tree_with_slice = optimize_code(code, uniformsize(code, 2),
	TreeSA(sc_target=28, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0, nslices=5));

julia> optcode_kahypar = optimize_code(code, uniformsize(code, 2), 
	KaHyParBipartite(sc_target=30, max_group_size=50));

julia> optcode_sa = optimize_code(code, uniformsize(code, 2),
	SABipartite(sc_target=30, max_group_size=50));

julia> timespace_complexity(code, uniformsize(code, 2))
(200.0, 0.0)

julia> timespace_complexity(optcode_kahypar, uniformsize(code, 2))
(39.00774639886569, 28.0)

julia> timespace_complexity(optcode_sa, uniformsize(code, 2))
(39.09524927371961, 28.0)

julia> timespace_complexity(optcode_tree, uniformsize(code, 2))
(31.13883492534964, 27.0)

julia> timespace_complexity(optcode_tree_with_slice, uniformsize(code, 2))
(31.292828948918775, 22.0)
```

It is already a part of [`OMEinsum`](https://github.com/under-Peter/OMEinsum.jl):

```julia
julia> using OMEinsum

julia> code = ein"ij, jk, kl, il->"
ij, jk, kl, il -> 

julia> optimize_code(code, uniformsize(code, 2), TreeSA())
SlicedEinsum{Char, NestedEinsum{DynamicEinCode{Char}}}(Char[], ki, ki -> 
├─ jk, ij -> ki
│  ├─ jk
│  └─ ij
└─ kl, il -> ki
   ├─ kl
   └─ il
)
```

## Acknowledge OMEinsum/OMEinsumContractionOrders
Since we do not have a paper to cite, we wish anyone who finds this package useful in his research can post a comment under this issue:
https://github.com/TensorBFS/OMEinsumContractionOrders.jl/issues/21

## References

If you find this package useful in your research, please cite the following papers
```
@misc{Liu2022,
  doi = {10.48550/ARXIV.2205.03718},
  url = {https://arxiv.org/abs/2205.03718},
  author = {Liu, Jin-Guo and Gao, Xun and Cain, Madelyn and Lukin, Mikhail D. and Wang, Sheng-Tao},
  keywords = {Statistical Mechanics (cond-mat.stat-mech), FOS: Physical sciences, FOS: Physical sciences},
  title = {Computing solution space properties of combinatorial optimization problems via generic tensor networks},
  publisher = {arXiv},
  year = {2022},
  copyright = {Creative Commons Attribution 4.0 International}
}
```

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

## Multi-GPU computation
Check this Gist:

https://gist.github.com/GiggleLiu/d5b66c9883f0c5df41a440589983ab99

## Authors
Jin-Guo Liu and Pan Zhang
