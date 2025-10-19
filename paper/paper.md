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
  - name: Xuan-Zhao Gao
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
 - name: <Xuan-Zhao fill in>
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

A _tensor network_ is a mathematical framework that represents multilinear algebra operations as intuitive graphical structures, where tensors become nodes and shared indices become connecting edges. This diagrammatic approach transforms complex high-dimensional contractions into accessible visual networks that expose underlying computational structure.

The framework exhibits remarkable universality, emerging across diverse domains: _einsum_ notation [@Harris2020] in numerical computing, _factor graphs_ [@Bishop2006] in probabilistic inference, _sum-product networks_ in machine learning, and _junction trees_ [@Villescas2023] in graphical models. Tensor networks have revolutionized quantum circuit simulation [@Markov2008], quantum error correction [@Piveteau2024], neural network compression [@Qing2024], and strongly correlated quantum materials [@Haegeman2016].

The computational cost of tensor network contraction depends critically on the chosen *contraction order*—the sequence in which pairwise tensor multiplications are performed. This order can be represented as a binary tree where leaves correspond to input tensors and internal nodes represent intermediate results.

Consider the contraction `ein"ab,bc,cd->ad"`, which admits multiple valid orderings with dramatically different costs:

The left ordering proves dramatically superior: it achieves $O(n^3)$ time and $O(n^2)$ space complexity by first contracting compatible matrices. The right ordering creates a $O(n^4)$ intermediate tensor through an inefficient Kronecker product, illustrating how ordering choice can determine computational feasibility.

Finding the globally optimal contraction order constitutes an NP-complete optimization problem [@Markov2008]. Fortunately, near-optimal solutions often suffice for practical applications and can be obtained efficiently through sophisticated heuristic methods. Modern optimization algorithms have achieved remarkable scalability, successfully handling tensor networks with over $10^4$ tensors [@Gray2021; @Roa2024].

The optimal contraction order has a deep mathematical connection to the _tree decomposition_ [@Markov2008] of the tensor network's line graph.

**Definition (Tree decomposition and treewidth):** A _tree decomposition_ of a (hyper)graph $G=(V,E)$ is a tree $T=(B,F)$ where each node $B_i \in B$ contains a subset of vertices in $V$ (called a "bag"), satisfying:

1. Every vertex $v \in V$ appears in at least one bag.
2. For each (hyper)edge $e \in E$, there exists a bag containing all vertices in $e$.
3. For each vertex $v \in V$, the bags containing $v$ form a connected subtree of $T$.

The _width_ of a tree decomposition is the size of its largest bag minus one. The _treewidth_ of a graph is the minimum width among all possible tree decompositions.


The line graph of a tensor network is a graph where vertices represent indices and edges represent tensors sharing those indices. The relationship between a tensor network's contraction order and the tree decomposition of its line graph can be understood through several key correspondences:

- Each leg (index) in the tensor network becomes a vertex in the line graph, while each tensor becomes a hyperedge connecting multiple vertices.
- The tree decomposition's first two requirements ensure that all tensors are accounted for in the contraction sequence - each tensor must appear in at least one bag, with each bag representing a contraction step.
- The third requirement of the tree decomposition maps to an important constraint in tensor contraction: an index can only be eliminated after considering all tensors connected to it.
- For tensor networks with varying index dimensions, we can extend this relationship to weighted tree decompositions, where vertex weights correspond to the logarithm of the index dimensions.

The figure below illustrates these concepts with (a) a tensor network containing four tensors $T_1$, $T_2$, $T_3$ and $T_4$ and eight indices labeled $A$ through $H$, (b) its corresponding line graph, and (c) a tree decomposition of that line graph.

![(a) A tensor network. (b) A line graph for the tensor network. Labels are connected if and only if they appear in the same tensor. (c) A tree decomposition (T. D.) of the line graph. \label{fig:mainfig}](figures/mainfig.pdf)

The tree decomposition in \autoref{fig:mainfig}(c) consists of 6 bags, each containing at most 3 indices, indicating that the treewidth of the tensor network is 2. The tensors $T_1$, $T_2$, $T_3$ and $T_4$ are contained in bags $B_1$, $B_5$, $B_6$ and $B_2$ respectively. Following the tree structure, we perform the contraction from the leaves. First, we contract bags $B_1$ and $B_2$ into $B_3$, yielding an intermediate tensor $I_{14} = T_1 * T_4$ (where "$*$" denotes tensor contraction) with indices $B$ and $E$. Next, we contract bags $B_5$ and $B_6$ into $B_4$, producing another intermediate tensor $I_{23} = T_2 * T_3$ also with indices $B$ and $E$. Finally, contracting $B_3$ and $B_4$ yields the desired scalar result.

Finding the optimal contraction order is almost equivalent to finding the minimal-width tree decomposition of the line graph.
The log time complexity for the bottleneck contraction corresponds to the largest bag size in the tree decomposition.
The log space complexity is equivalent to the largest separator (the set of vertices connecting two bags) size in the tree decomposition.


[JG: TODO: We need a benchmark with CoTengra optimizer (Richard fill in)]

# Usage Example
- _Remark_: 1. Basic usage. From contraction pattern representation to optimized contraction order, introduce the data structures and algorithms. Conversion between contraction graph and treewidth.

`OMEinsum` provides multiple heuristic methods for finding the optimal contraction order. They are implemented in the dependency `OMEinsumContractionOrders`. To demonstrate the usage, we first generate a large enough random tensor network with the help of the `Graphs` package.

```julia
julia> using OMEinsum, Graphs

julia> function demo_network(n::Int)
           g = random_regular_graph(n, 3)
           code = EinCode([[e.src, e.dst] for e in edges(g)], Int[])
           sizes = uniformsize(code, 2)
           tensors = [randn([sizes[leg] for leg in ix]...) for ix in getixsv(code)]
           return code, tensors, sizes
       end
demo_network (generic function with 1 method)

julia> code, tensors, sizes = demo_network(100);

julia> contraction_complexity(code, sizes)
Time complexity: 2^100.0
Space complexity: 2^0.0
Read-write complexity: 2^9.231221180711184
```

We first generate a random 3-regular graph with 100 vertices. Then we associate each vertex with a binary variable and each edge with a tensor of size $2 \times 2$. The time complexity without contraction order optimization is $2^{100}$, which is equivalent to brute-force. The order can be optimized with the `optimize_code` function.

```julia
julia> optcode = optimize_code(code, sizes, TreeSA());

julia> cc = contraction_complexity(optcode, sizes)
Time complexity: 2^17.241796993093228
Space complexity: 2^13.0
Read-write complexity: 2^16.360864226366807
```
The `optimize_code` function takes three inputs: the `EinCode` object, the tensor sizes, and the contraction order solver. It returns a `NestedEinsum` object of time complexity $\approx 2^{17.2}$. It is much smaller than the number of vertices. It is a very reasonable number because the treewidth of a 3-regular graph is approximately upper bounded by $1/6$ of the number of vertices [@Fomin2006].

We can use the `slice_code` function to reduce the space complexity.
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
The `slice_code` function takes three inputs: the `NestedEinsum` object, the tensor sizes, and the slicing strategy. Here, we use the `TreeSASlicer` with the `ScoreFunction` to reduce the space complexity by 3. The result type is `SlicedEinsum`, which contains a `slicing` field for storing the slice indices. After slicing, the space complexity is reduced by $3$, while the time complexity is only slightly increased. The usage of `SlicedEinsum` is the same as the `NestedEinsum` object.

```julia
julia> @assert sliced_code(tensors...) ≈ optcode(tensors...)
```

Supported solvers include:

| Optimizer | Description |
| :----------- | :------------- |
| `GreedyMethod` | Fast, but poor contraction order |
| `TreeSA` | Reliable, local search based optimizer [@Kalachev2021], but is a bit slow |
| `HyperND` | Nested dissection algorithm, similar to `KaHyParBipartite`. Requires importing either `KaHyPar` or `Metis`. |
| `KaHyParBipartite` and `SABipartite` | Graph bipartition based, suited for large tensor networks [@Gray2021], requires using `KaHyPar` package. Alternatively, a simulated annealing bipartition method is provided in `SABipartite`. |
| `ExactTreewidth` (alias of `Treewidth{RuleReduction{BT}}`) | Exact, but takes exponential time [@Bouchitte2001], based on package `TreeWidthSolver`. |
| `Treewidth` | Tree width solver based, based on package `CliqueTrees`, performance is elimination algorithm dependent. |

There is a tradeoff between the time and the quality of the contraction order. \autoref{fig:timespace} shows the Pareto front of the multi-objective optimization of the time to optimize the contraction order and the time to contract the tensor network.

![A benchmark result of the time and space trade-off. \label{fig:timespace}](figures/timespace.svg)


[JG: TODO: We should redirectly use the existing materials in the examples folder.]

# Acknowledgments

This work is partially funded by the Google Summer of Code 2024 project.
We extend our gratitude to Feng Pan for his insightful discussions on the slicing technique.

# References
