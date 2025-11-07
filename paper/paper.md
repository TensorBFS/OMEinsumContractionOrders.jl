---
title: 'OMEinsumContractionOrders: A Julia package for tensor network contraction order optimization'
tags:
  - Julia
  - tensor networks
  - contraction order optimization
authors:
  - name: Jin-Guo Liu
    orcid: 0000-0003-1635-2679
    equal-contrib: true
    affiliation: 1
  - name: Xuanzhao Gao
    orcid: 
    equal-contrib: true
    affiliation: 2
  - name: Richard Samuelson
    orcid: 
    equal-contrib: true
    affiliation: 3
affiliations:
 - name: Hong Kong University of Science and Technology (Guangzhou)
   index: 1
 - name: Center of Computational Mathematics, Flatiron Institute
   index: 2
 - name: <Richard fill in>
   index: 3
date: 19 October 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.

# Citations to entries in paper.bib should be in
# [rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
# format.
# 
# If you want to cite a software repository URL (e.g. something on GitHub without a preferred
# citation) then you can do it with the example BibTeX entry below for @fidgit.
# 
# For a quick reference, the following citation commands can be used:
# - `@author:2001`  ->  "Author et al. (2001)"
# - `[@author:2001]` -> "(Author et al., 2001)"
# - `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures can be included like this:
# ![Caption for example figure.\label{fig:example}](figure.png)
# and referenced from text using \autoref{fig:example}.
# 
# Figure sizes can be customized by adding an optional second parameter:
# ![Caption for example figure.](figure.png){ width=20% }
---

# Statement of need

`OMEinsumContractionOrders` (One More Einsum Contraction Orders, or OMECO) is a Julia package [@bezanson2012julia] that implements state-of-the-art algorithms for optimizing tensor network contraction orders. This paper presents its key features, integration with the Julia ecosystem, and performance benchmarks.

![The ecosystem built around `OMEinsumContractionOrders` and its dependencies. OMECO serves as a core component of the tensor network contractor `OMEinsum`, which powers applications including `Yao` (quantum simulation), `TensorQEC` (quantum error correction), `TensorInference` (probabilistic inference), `GenericTensorNetworks` and `TensorBranching` (combinatorial optimization).\label{fig:structure}](figures/structure.pdf){ width=80% }

A _tensor network_ is a mathematical framework that represents multilinear algebra operations as graphical structures, where tensors are nodes and shared indices are edges. This diagrammatic approach transforms complex high-dimensional contractions into visual networks that expose underlying computational structure.

The framework has remarkable universality across diverse domains: _einsum_ notation [@Harris2020] in numerical computing, _factor graphs_ [@Bishop2006] in probabilistic inference, _sum-product networks_ in machine learning, and _junction trees_ [@Villescas2023] in graphical models. Tensor networks have enabled breakthroughs in quantum circuit simulation [@Markov2008], quantum error correction [@Piveteau2024], neural network compression [@Qing2024], strongly correlated quantum materials [@Haegeman2016], and combinatorial optimization problems [@Liu2023]. These applications are reflected in the ecosystem built around OMECO, as illustrated in \autoref{fig:structure}.

The computational cost of tensor network contraction depends critically on the *contraction order*—the sequence in which pairwise tensor multiplications are performed. This order can be represented as a binary tree where leaves correspond to input tensors and internal nodes represent intermediate results. The optimization objective balances multiple complexity measures through the cost function:

$$
\mathcal{L} = w_\text{t} \cdot \text{tc} + w_\text{s} \cdot \max(0, \text{sc} - \text{sc}_{\rm target}) + w_\text{rw} \cdot \text{rwc},
$$
where $w_\text{t}$, $w_\text{s}$, and $w_\text{rw}$ represent weights for time complexity (tc), space complexity (sc), and read-write complexity (rwc), respectively. In practice, memory access costs typically dominate computational costs, motivating $w_\text{rw} > w_\text{t}$. The space complexity penalty activates only when $\text{sc} > \text{sc}_{\rm target}$, allowing unconstrained optimization when memory fits within available device capacity.

Finding the optimal contraction order—even when minimizing only time complexity—is NP-complete [@Markov2008]. This optimization problem has a deep mathematical connection to _tree decomposition_ [@Markov2008] of the tensor network's line graph, where finding the optimal order corresponds to finding a weighted minimal-width tree decomposition. The logarithmic time complexity of the bottleneck contraction step equals the largest bag size in the tree decomposition, while the logarithmic space complexity equals the largest separator size (vertices shared between adjacent bags).

Despite this computational hardness, near-optimal solutions suffice for most practical applications and can be obtained efficiently through heuristic methods. Modern optimization algorithms have achieved remarkable scalability, handling tensor networks with over $10^4$ tensors [@Gray2021; @Roa2024].

OMECO implements several optimization algorithms with complementary performance characteristics:

| Optimizer | Description |
| :----------- | :------------- |
| `GreedyMethod` | Fast greedy heuristic with modest solution quality |
| `TreeSA` | Reliable simulated annealing optimizer [@Kalachev2021] with high-quality solutions |
| `PathSA` | Simulated annealing optimizer for path decomposition |
| `HyperND` | Nested dissection algorithm for hypergraphs, requires `KaHyPar` or `Metis` |
| `KaHyParBipartite` | Graph bipartition method for large tensor networks [@Gray2021], requires `KaHyPar` |
| `SABipartite` | Simulated annealing bipartition method, pure Julia implementation |
| `ExactTreewidth` | Exact algorithm with exponential runtime [@Bouchitte2001], based on `TreeWidthSolver` |
| `Treewidth` | Clique tree elimination methods from `CliqueTrees` package |

The algorithms `HyperND`, `Treewidth`, and `ExactTreewidth` operate on the tensor network's line graph and utilize the `CliqueTrees` and `TreeWidthSolver` packages, as illustrated in \autoref{fig:structure}. Additionally, the `PathSA` optimizer implements path decomposition by constraining contraction orders to path graphs, serving as a variant of `TreeSA`.

These methods balance optimization time against solution quality. \autoref{fig:sycamore} displays benchmark results for the Sycamore quantum supremacy circuit, highlighting the Pareto front where contraction order quality is balanced with optimization runtime.

![Benchmark results for contraction order optimization on the Sycamore quantum circuit tensor network (Intel Xeon Gold 6226R CPU @ 2.90GHz, single-threaded). The $x$-axis shows contraction cost, $y$-axis shows optimization time. Each point represents a different optimizer configuration tested with varying parameters. `TreeSA` and `HyperND` achieve the lowest contraction costs, while `GreedyMethod` offers the fastest optimization time. Detailed parameter configurations are documented in `benchmark-details.md`. \label{fig:sycamore}](figures/sycamore.pdf){ width=80% }


[JG: TODO: Please also cite CliqueTree paper.]


OMECO has been integrated into the `OMEinsum` package and powers several downstream applications: `Yao` [@Luo2020] for quantum circuit simulation, `GenericTensorNetworks` [@Liu2023] and `TensorBranching` (TODO: add citation) for combinatorial optimization, `TensorInference` [@Roa2023] for probabilistic inference, and `TensorQEC` (TODO: add citation) for quantum error correction. This infrastructure is expected to benefit other applications requiring tree or path decomposition, such as polynomial optimization [@Magron2021].

# Usage Example

OMECO provides two main functions: `optimize_code` for finding optimal contraction orders, and `slice_code` for trading time complexity for reduced space complexity through the slicing technique.

To demonstrate basic usage, we generate a random 3-regular graph with 100 vertices using the `Graphs` package, associating each vertex with a binary variable and each edge with a $2 \times 2$ tensor.

```julia
julia> using Graphs: random_regular_graph, edges, vertices

julia> using OMEinsumContractionOrders: EinCode, uniquelabels, contraction_complexity, optimize_code, TreeSA, slice_code, TreeSASlicer, ScoreFunction

julia> function demo_network(n::Int)
           g = random_regular_graph(n, 3)
           code = EinCode([[e.src, e.dst] for e in edges(g)], Int[])
           sizes = Dict(i=>2 for i in uniquelabels(code))
           tensors = [randn([sizes[index] for index in ix]...) for ix in code.ixs]
           return code, tensors, sizes
       end
demo_network (generic function with 1 method)

julia> code, tensors, sizes = demo_network(100);
```
The tensor network topology is represented by an `EinCode` object with two fields: `ixs` (a vector of index vectors for each input tensor) and `iy` (output indices). This structure defines a hypergraph with potentially open edges. Combining this hypergraph with tensor sizes determines the contraction complexity.

```julia
julia> contraction_complexity(code, sizes)
Time complexity: 2^100.0
Space complexity: 2^0.0
Read-write complexity: 2^9.231221180711184
```

The return type contains three fields (`tc`, `sc`, `rwc`) for time, space, and read-write complexity. Without optimization, the time complexity is $2^{100}$, equivalent to brute-force enumeration.

We now use the `TreeSA` optimizer to find an improved contraction order.

```julia
julia> optcode = optimize_code(code, sizes, TreeSA(; score=ScoreFunction(tc_weight=1.0, sc_weight=1.0, rw_weight=10.0)));

julia> cc = contraction_complexity(optcode, sizes)
Time complexity: 2^17.241796993093228
Space complexity: 2^13.0
Read-write complexity: 2^16.360864226366807
```
The `optimize_code` function takes three arguments: the `EinCode` object, tensor size dictionary, and optimizer configuration. It returns a `NestedEinsum` object specifying the contraction tree with three fields: `args` (child nodes), `tensorindex` (input tensor index for leaf nodes), and `eins` (einsum notation for the node). The time complexity $\approx 2^{17.2}$ is dramatically improved from the original $2^{100}$. This result aligns with theory, as the treewidth of a 3-regular graph is approximately upper bounded by $1/6$ of the number of vertices [@Fomin2006]. The `score` keyword argument configures the cost function weights; here we set the read-write weight to 10× the time weight, reflecting the higher cost of memory access.

Space complexity can be further reduced using `slice_code`, which implements the slicing technique to trade time for space.
```julia
julia> sliced_code = slice_code(optcode, sizes, TreeSASlicer(score=ScoreFunction(sc_target=cc.sc-3)));

julia> sliced_code.slicing
3-element Vector{Int64}:
 14
 76
 60

julia> contraction_complexity(sliced_code, sizes)
Time complexity: 2^17.800899899920303
Space complexity: 2^10.0
Read-write complexity: 2^17.199595668955244
```
The `slice_code` function takes the `NestedEinsum` object, tensor sizes, and slicing strategy, returning a `SlicedEinsum` object with two fields: `slicing` (sliced indices) and `eins` (a `NestedEinsum` object). Using `TreeSASlicer`, we reduce space complexity by $2^3$ (from $2^{13}$ to $2^{10}$) with only a modest increase in time complexity. The resulting `SlicedEinsum` object maintains the same interface as `NestedEinsum` for contraction evaluation.

```julia
julia> @assert sliced_code(tensors...) ≈ optcode(tensors...)
```


[JG: TODO: Mention the API to convert between contraction graph and treewidth. (Xuan-Zhao fill in), Remove?]

[JG: TODO: Show a plot about using slicing to reduce the space complexity (based on the above example). (Xuan-Zhao fill in)]

Real-world examples demonstrating applications to quantum circuit simulation, combinatorial optimization, and probabilistic inference are available in the [OMEinsumContractionOrdersBenchmark](https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark) repository. Optimizer performance is highly problem-dependent, with no single algorithm dominating across all metrics and graph topologies.

# Acknowledgments

This work was partially funded by Google Summer of Code 2024 and the Open Source Promotion Plan (OSPP 2023).
We thank Feng Pan for insightful discussions and code contributions on the slicing technique.

# References
