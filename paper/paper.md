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

# Summary

`OMEinsumContractionOrders` (One More Einsum Contraction Orders, or OMECO) is a Julia package [@bezanson2012julia] that implements state-of-the-art algorithms for optimizing tensor network contraction orders.
OMECO is designed to search for near-optimal contraction orders for exact tensor network contraction, and provides a comprehensive suite of optimization algorithms for tensor network contraction orders, including greedy heuristics, simulated annealing, and tree width solvers.
In this paper, we present the key features of OMECO, its integration with the Julia ecosystem, and performance benchmarks.


# Statement of need

A _tensor network_ is a mathematical framework that represents multilinear algebra operations as graphical structures, where tensors are nodes and shared indices are edges. 
This diagrammatic approach transforms complex high-dimensional contractions into visual networks that expose underlying computational structure.

The framework has remarkable universality across diverse domains: _einsum_ notation [@Harris2020] in numerical computing, _factor graphs_ [@Bishop2006] in probabilistic inference, _sum-product networks_ in machine learning, and _junction trees_ [@Villescas2023] in graphical models. Tensor networks have enabled breakthroughs in quantum circuit simulation [@Markov2008], quantum error correction [@Piveteau2024], neural network compression [@Qing2024], strongly correlated quantum materials [@Haegeman2016], and combinatorial optimization problems [@Liu2023]. 

The computational cost of tensor network contraction depends critically on the *contraction order*—the sequence in which pairwise tensor multiplications are performed. This order can be represented as a binary tree where leaves correspond to input tensors and internal nodes represent intermediate results. The optimization objective balances multiple complexity measures through the cost function:
$$
\mathcal{L} = w_\text{t} \cdot \text{tc} + w_\text{s} \cdot \max(0, \text{sc} - \text{sc}_{\rm target}) + w_\text{rw} \cdot \text{rwc},
$$
where $w_\text{t}$, $w_\text{s}$, and $w_\text{rw}$ represent weights for time complexity (tc), space complexity (sc), and read-write complexity (rwc), respectively. In practice, memory access costs typically dominate computational costs, motivating $w_\text{rw} > w_\text{t}$. The space complexity penalty activates only when $\text{sc} > \text{sc}_{\rm target}$, allowing unconstrained optimization when memory fits within available device capacity.

Finding the optimal contraction order—even when minimizing only time complexity—is NP-complete [@Markov2008]. 
<!-- This optimization problem has a deep mathematical connection to _tree decomposition_ [@Markov2008] of the tensor network's line graph, where finding the optimal order corresponds to finding a weighted minimal-width tree decomposition. The logarithmic time complexity of the bottleneck contraction step equals the largest bag size in the tree decomposition, while the logarithmic space complexity equals the largest separator size (vertices shared between adjacent bags). -->
Algorithms for finding near-optimal contraction orders have been developed and achieve impressive scalability [@Gray2021; @Roa2024], handling tensor networks with over $10^4$ tensors.
However, an efficient and reliable implementation of these methods in Julia is still missing.

OMECO addresses this gap by offering a unified and extensible interface to a comprehensive suite of optimization algorithms for tensor network contraction orders, including greedy heuristics, simulated annealing, and tree-width-based solvers.
OMECO has been integrated into the `OMEinsum` package and powers several downstream applications: `Yao` [@Luo2020] for quantum circuit simulation, `GenericTensorNetworks` [@Liu2023] and `TensorBranching` (TODO: add citation) for combinatorial optimization, `TensorInference` [@Roa2023] for probabilistic inference, and `TensorQEC` (TODO: add citation) for quantum error correction. 
These applications are reflected in the ecosystem built around OMECO, as illustrated in \autoref{fig:structure}.
This infrastructure is expected to benefit other applications requiring tree or path decomposition, such as polynomial optimization [@Magron2021].

![The ecosystem built around `OMEinsumContractionOrders` and its dependencies. OMECO serves as a core component of the tensor network contractor `OMEinsum`, which powers applications including `Yao` (quantum simulation), `TensorQEC` (quantum error correction), `TensorInference` (probabilistic inference), `GenericTensorNetworks` and `TensorBranching` (combinatorial optimization).\label{fig:structure}](figures/structure.pdf){ width=80% }

# Features and benchmarks

The major feature of OMECO is contraction order optimization.
OMECO provides several algorithms with complementary performance characteristics that can be simply called by the `optimize_code` function:

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
Real-world examples demonstrating applications to quantum circuit simulation, combinatorial optimization, and probabilistic inference are available in the [OMEinsumContractionOrdersBenchmark](https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark) repository. Optimizer performance is highly problem-dependent, with no single algorithm dominating across all metrics and graph topologies.

![Time complexity (left) and space complexity (right) benchmark results for contraction order optimization on the Sycamore quantum circuit tensor network (Intel Xeon Gold 6226R CPU @ 2.90GHz, single-threaded). The $x$-axis shows contraction cost, $y$-axis shows optimization time. Each point represents a different optimizer configuration tested with varying parameters. `TreeSA` and `HyperND` achieve the lowest contraction costs, while `GreedyMethod` offers the fastest optimization time.\label{fig:sycamore}](figures/sycamore.pdf){ width=95% }


[JG: TODO: Please also cite CliqueTree paper.]

Optimizers prefixed with `cotengra_` are from the Python package cotengra [@Gray2021]; all others are OMECO implementations. For both optimization objectives (minimizing time and space complexity), OMECO optimizers dominate the Pareto front. Given sufficient optimization time, `TreeSA` consistently achieves the lowest time and space complexity. `GreedyMethod` provides the fastest optimization but yields suboptimal contraction orders, while `HyperND` offers a favorable balance between optimization time and solution quality.


OMECO has been integrated into the `OMEinsum` package and powers several downstream applications: `Yao` [@Luo2020] for quantum circuit simulation, `GenericTensorNetworks` [@Liu2023] and `TensorBranching` (TODO: add citation) for combinatorial optimization, `TensorInference` [@Roa2023] for probabilistic inference, and [`TensorQEC`](https://github.com/TensorBFS/TensorQEC.jl) for quantum error correction. This infrastructure is expected to benefit other applications requiring tree or path decomposition, such as polynomial optimization [@Magron2021].

The second feature of OMECO is index slicing, which is a technique to trade time complexity for reduced space complexity by looping over a subset of indices directly.
In OMECO, the interface to perform index slicing is `slice_code`, and currently we only support one Slicer, `TreeSASlicer`, which is implemented the dynamic slicing algorithm based on `TreeSA`.
An example is shown in \autoref{fig:slicing}, where the network of the Sycamore quantum circuit is sliced to reduce the space complexity from 52 to 31.

![Relation between the time complexity and the target space complexity for the `TreeSASlicer`. The original network is the Sycamore quantum circuit with an original space complexity of 52. \label{fig:slicing}](figures/sycamore_slicing.pdf){ width=40% }

# References

<!-- # Supporting
[JG: we will clean up this part in the future]

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
 -->

