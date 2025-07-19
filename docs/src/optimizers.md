# Choosing Optimizers

Supported solvers include:

| Optimizer | Description |
| :----------- | :------------- |
| [`GreedyMethod`](@ref Sec_GreedyMethod) | Fast, but poor contraction order |
| [`TreeSA`](@ref Sec_TreeSA) | Reliable, local search based optimizer [^Kalachev2021], but is a bit slow |
| [`HyperND`](@ref Sec_HyperND) | Nested dissection algorithm, similar to [`KaHyParBipartite`](@ref). Requires importing either [`KaHyPar`](https://github.com/kahypar/KaHyPar.jl) or [`Metis`](https://github.com/JuliaSparse/Metis.jl). |
| [`KaHyParBipartite` and `SABipartite`](@ref Sec_Bipartite) | Graph bipartition based, suited for large tensor networks [^Gray2021], requires using [`KaHyPar`](https://github.com/kahypar/KaHyPar.jl) package. Alternatively, a simulated annealing bipartition method is provided in [`SABipartite`](@ref). |
| [`ExactTreewidth`](@ref Sec_ExactTreewidth) (alias of `Treewidth{RuleReduction{BT}}`) | Exact, but takes exponential time [^Bouchitté2001], based on package [`TreeWidthSolver`](https://github.com/ArrogantGao/TreeWidthSolver.jl). |
| [`Treewidth`](@ref Sec_Treewidth) | Tree width solver based, based on package [`CliqueTrees`](https://github.com/AlgebraicJulia/CliqueTrees.jl), performance is elimination algorithm dependent. |

There is a tradeoff between the time and the quality of the contraction order. The following figure shows the Pareto front of the multi-objective optimization of the time to optimize the contraction order and the time to contract the tensor network.

![](assets/tradeoff.svg)

Among these methods, the `ExactTreewidth` method produces the lowest treewidth, but it does not scale up to tensor networks with more than 50 tensors. The `TreeSA` is the second best in terms of the treewidth. It works well in most cases, and supports [slicing](@ref "Reduce space complexity by slicing").
The only limitation is that it is a bit slow.
For application sensitive to overhead, the `GreedyMethod` and `Treewidth` method (blue region) are recommended.
The `Treewidth` method is a zoo of methods provided by the package [`CliqueTrees`](https://github.com/AlgebraicJulia/CliqueTrees.jl), which is a collection of methods for finding the approximate tree decomposition of a graph. Most of them have similar performance with the `GreedyMethod`, and most of them are very efficient. The `HyperND` method has a very good overall performance in benchmarks (to be added), and it is much faster than the `TreeSA` method. It relies on the `KaHyPar` package, which is platform picky.

## [`GreedyMethod`](@id Sec_GreedyMethod)
Implemented as [`GreedyMethod`](@ref) in the package.
The Greedy method is one of the simplest and fastest method for optimizing the contraction order. The idea is to greedily select the pair of tensors with the smallest cost to contract at each step.
The cost is defined as:
```math
L = \text{size}(\text{out}) - α \times (\text{size}(\text{in}_1) + \text{size}(\text{in}_2))
```
where $\text{out}$ is the output tensor, and $\text{in}_1$ and $\text{in}_2$ are the input tensors. $α$ is a hyperparameter, which is set to $0.0$ by default, meaning that we greedily select the pair of tensors with the smallest size of the output tensor. For $\alpha = 1$, the size increase in each step is greedily optimized.

The greedy method implemented in this package uses the priority queue to select the pair of tensors with the smallest cost to contract at each step. The time complexity is $O(n^2 \log n)$ for $n$ tensors, since in each of the $n$ steps, we pick the pair with the smallest cost in $O(n \log n)$ time.

## [`TreeSA`](@id Sec_TreeSA)

Implemented as [`TreeSA`](@ref) in the package.
The local search method [^Kalachev2021] is a heuristic method based on the idea of simulated annealing.
The method starts from a random contraction order and then applies the following four possible transforms as shown in the following figure

![](assets/treesa.svg)

They correspond to the different ways to contract three sub-networks:
```math
(A * B) * C = (A * C) * B = (C * B) * A, \\
A * (B * C) = B * (A * C) = C * (B * A),
```
where we slightly abuse the notation "$*$" to denote the tensor contraction, and $A, B, C$ are the sub-networks to be contracted.
Due to the commutative property of the tensor contraction, such transformations do not change the result of the contraction.
Even through these transformations are simple, all possible contraction orders can be reached from any initial contraction order.
The local search method starts from a random contraction tree.
In each step, the above rules are randomly applied to transform the tree and then the cost of the new tree is evaluated, which is defined as
```math
\mathcal{L} = \text{tc} + w_s \times \text{sc} + w_{\text{rw}} \times \text{rwc},
```
where $w_s$ and $w_{\text{rw}}$ are the weights of the space complexity and read-write complexity compared to the time complexity, respectively.
The optimal choice of weights depends on the specific device and tensor contraction algorithm. One can freely tune the weights to achieve a best performance for their specific problem.
Then the transformation is accepted with a probability given by the Metropolis criterion, which is
```math
p_{\text{accept}} = \min(1, e^{-\beta \Delta \mathcal{L}}),
```
where $\beta$ is the inverse temperature, and $\Delta \mathcal{L}$ is the difference of the cost of the new and old contraction trees.
During the process, the temperature is gradually decreased, and the process stop when the temperature is low enough.

This simple algorithm is flexible and works suprisingly well in most cases. Given enough time, it almost always finds the contraction order with the lowest cost.
The only weakness is the runtime. It usually takes minutes to optimize a network with 1k tensors.
Additionally, the `TreeSA` method supports using the [slicing](@ref "Reduce space complexity by slicing") technique to reduce the space complexity.

## [`HyperND`](@id Sec_HyperND)

Implemented as [`HyperND`](@ref) in the package.

## [`KaHyParBipartite` and `SABipartite`](@id Sec_Bipartite)

Implemented as [`KaHyParBipartite`](@ref) and [`SABipartite`](@ref) in the package. These two methods are based on the graph bipartition, 
```math
\min \sum_{e \in \text{cut}}\omega(e)\\
\text{s.t.} \quad c(V_i) \leq (1+\epsilon) \left\lceil \frac{c(V)}{2} \right\rceil.
```
where $\text{cut}$ is the set of hyperedges cut by the partition, $\omega(e)$ is the weight of the hyperedge $e$, $V_i$ is the $i$-th part of the partition, $c(V)$ is the node weight of the partition, and $\epsilon$ is the imbalance parameter.
It tries to minimize the cut crossing two blocks, while the size of each block is balanced.
The algorithm is implemented in the package [KaHyPar.jl](https://github.com/kahypar/KaHyPar.jl)[^Schlag2021], and the implementation in [OMEinsumContractionOrders.jl](https://github.com/TensorBFS/OMEinsumContractionOrders.jl) is mainly based on it.

```@raw html
<img src="https://cloud.githubusercontent.com/assets/484403/25314222/3a3bdbda-2840-11e7-9961-3bbc59b59177.png" alt="alt text" width="50%" height="50%"><img src="https://cloud.githubusercontent.com/assets/484403/25314225/3e061e42-2840-11e7-860c-028a345d1641.png" alt="alt text" width="50%" height="50%">
```

Some platforms may have some issues in installation of `KaHyPar`, please refer to [#12](https://github.com/kahypar/KaHyPar.jl/issues/12) and [#19](https://github.com/kahypar/KaHyPar.jl/issues/19).
The [`SABipartite`](@ref) is a simulated annealing based alternative to [`KaHyParBipartite`](@ref), it can produce similar results while being much more costly.

Note: Benchmarks (to be added) show that the later implementation of [`HyperND`](@ref) method is better and faster. These two methods are no longer the first choice.

## [`ExactTreewidth`](@id Sec_ExactTreewidth)

Implemented as [`ExactTreewidth`](@ref) in the package.
This method is supported by the [Google Summer of Code 2024](https://summerofcode.withgoogle.com) project ["Tensor network contraction order optimization and visualization"](https://summerofcode.withgoogle.com/programs/2024/projects/B8qSy9dO) released by **The Julia Language**.
In this project, we developed [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl), which implements the Bouchitté–Todinca algorithm[^Bouchitté2001]. It later becomes a backend of [OMEinsumContracionOrders.jl](https://github.com/TensorBFS/OMEinsumContractionOrders.jl).

The Bouchitté–Todinca (BT) algorithm [^Bouchitté2001] is a method for calculating the treewidth of a graph exactly. It makes use of the theory of minimal triangulations, characterizing the minimal triangulations of a graph via objects called minimal separators and potential maximal cliques of the graph.
The BT algorithm has a time complexity of $O(|\Pi|nm)$, which are dependent on the graph structure. (TODO: add more details of the algorithm complexity, what is it suited for?).

The blog post [Finding the Optimal Tree Decomposition with Minimal Treewidth - Xuan-Zhao Gao](https://arrogantgao.github.io/blogs/treewidth/) has a more detailed description of this method.

## [`Treewidth`](@id Sec_Treewidth)

Implemented as [`Treewidth`](@ref) in the package.

## Exhaustive Search (planned)

The exhaustive search [^Robert2014] is a method to get the exact optimal contraction complexity.
There are three different ways to implement the exhaustive search:
* **Depth-first constructive approach**: in each step, choose a pair of tensors to contract a new tensor until all tensors are contracted, and then iterate over all possible contraction sequences without duplication. Note the cheapest contraction sequence thus found.
* **Breadth-first constructive approach**: the breadth-first method construct the set of intermediate tensors by contracting $c$ tensors ($c \in [1, n - 1]$, where $n$ is the number of tensors) in each step, and record the optimal cost for constructing each intermediate tensor. Then in the last step, the optimal cost for contracting all $n$ tensors is obtained.
* **Dynamic programming**: in each step, consider all bipartition that split the tensor network into two parts, if the optimal cost for each part is not recorded, further split them until the cost has been already obtained or only one tensor is left. Then combine the two parts and record the optimal cost of contracting the sub-networks. In this end the optimal cost for the whole network is obtained.
In more recent work [^Robert2014], by reordering the search process in favor of cheapest-first and excluding large numbers of outer product contractions which are shown to be unnecessary, the efficiency of the exhaustive search has been greatly improved.
The method has been implemented in [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl).

## Performance Benchmark

The following figure shows the performance of the different optimizers on the Sycamore 53-20-0 benchmark. This network is for computing the expectation value of a quantum circuit. Its optimal space complexity is known to be 52.

![](https://raw.githubusercontent.com/TensorBFS/OMEinsumContractionOrdersBenchmark/refs/heads/main/figures/sycamore_53_20_0.svg)

- The x-axis (contraction cost) here is defined by the space complexity (the size of the largest intermediate tensor).
- The y-axis (time) is the time taken to find the optimal contraction order.

By checking the Pareto front, we can see that the `TreeSA`, `HyperND` and `GreedyMethod` method are among the best. If you want the speed to find the optimal contraction order, the `GreedyMethod` is the best choice. If you want the quality of the contraction order, the `TreeSA` and `HyperND` method are the best choices. `HyperND` is basically a better version of `KaHyParBipartite` method, we are going to deprecate the `KaHyParBipartite` (and also `SABipartite`) method in the future.

More benchmarks, and the platform details can be found in the [benchmark repo](https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark). Or just check this full report: [Benchmark Report](https://raw.githubusercontent.com/TensorBFS/OMEinsumContractionOrdersBenchmark/refs/heads/main/report.pdf).

## References

[^Bouchitté2001]: Bouchitté, V., Todinca, I., 2001. Treewidth and Minimum Fill-in: Grouping the Minimal Separators. SIAM J. Comput. 31, 212–232. https://doi.org/10.1137/S0097539799359683
[^Gray2021]: Gray, Johnnie, and Stefanos Kourtis. "Hyper-optimized tensor network contraction." Quantum 5 (2021): 410.
[^Kalachev2021]: Kalachev, Gleb, Pavel Panteleev, and Man-Hong Yung. "Recursive multi-tensor contraction for XEB verification of quantum circuits." arXiv preprint arXiv:2108.05665 (2021).
[^Robert2014]: Pfeifer, R.N.C., Haegeman, J., Verstraete, F., 2014. Faster identification of optimal contraction sequences for tensor networks. Phys. Rev. E 90, 033315. https://doi.org/10.1103/PhysRevE.90.033315
[^Schlag2021]: Schlag, S., Heuer, T., Gottesbüren, L., Akhremtsev, Y., Schulz, C., Sanders, P., 2021. High-Quality Hypergraph Partitioning. https://doi.org/10.48550/arXiv.2106.08696
