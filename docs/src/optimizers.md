# Choosing Optimizers

Supported solvers include:

| Optimizer | Description |
| :----------- | :------------- |
| [`GreedyMethod`](@ref) | Fast, but poor resulting order |
| [`TreeSA`](@ref) | Reliable, local search based optimizer [^Kalachev2021], but is a bit slow |
| [`KaHyParBipartite`](@ref) and [`SABipartite`](@ref) | Graph bipartition based, suited for large tensor networks [^Gray2021], requires using [`KaHyPar`](https://github.com/kahypar/KaHyPar.jl) package |
| [`Treewidth`](@ref) | Tree width solver based, based on package [`CliqueTrees`](https://github.com/AlgebraicJulia/CliqueTrees.jl), performance is elimination algorithm dependent |
| [`ExactTreewidth`](@ref) (alias of `Treewidth{RuleReduction{BT}}`) | Exact, but takes exponential time [^Bouchitté2001], based on package [`TreeWidthSolver`](https://github.com/ArrogantGao/TreeWidthSolver.jl) |
| [`HyperND`](@ref) | Nested dissection algorithm, similar to [`KaHyParBipartite`](@ref). Requires imporing either [`KaHyPar`](https://github.com/kahypar/KaHyPar.jl) or [`Metis`](https://github.com/JuliaSparse/Metis.jl). |

The `KaHyParBipartite` is implemented as an extension. If you have issues in installing `KaHyPar`, please check these issues: [#12](https://github.com/kahypar/KaHyPar.jl/issues/12) and [#19](https://github.com/kahypar/KaHyPar.jl/issues/19).
Additionally, code simplifiers can be used to preprocess the tensor network to reduce the optimization time:

| Simplifier | Description |
| :----------- | :------------- |
| [`MergeVectors`](@ref) | Merges vector tensors with their neighbors |
| [`MergeGreedy`](@ref) | Greedily merges rank non-increasing tensors |

## References

[^Bouchitté2001]: Bouchitté, V., Todinca, I., 2001. Treewidth and Minimum Fill-in: Grouping the Minimal Separators. SIAM J. Comput. 31, 212–232. https://doi.org/10.1137/S0097539799359683
[^Gray2021]: Gray, Johnnie, and Stefanos Kourtis. "Hyper-optimized tensor network contraction." Quantum 5 (2021): 410.
[^Kalachev2021]: Kalachev, Gleb, Pavel Panteleev, and Man-Hong Yung. "Recursive multi-tensor contraction for XEB verification of quantum circuits." arXiv preprint arXiv:2108.05665 (2021).