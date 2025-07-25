"""
    HyperND(;
        dis = KaHyParND(),
        algs = (MF(), AMF(), MMD()),
        level = 6,
        width = 120,
        imbalances = 130:130,
        score = ScoreFunction(),
    )

Nested-dissection based optimizer. Recursively partitions a tensor network, then calls a
greedy algorithm on the leaves. The optimizer is run a number of times: once for each greedy
algorithm in `algs` and each imbalance value in `imbalances`. The recursion depth is controlled by
the parameters `level` and `width`.

The line graph is partitioned using the algorithm `dis`. OMEinsumContractionOrders currently supports two partitioning
algorithms, both of which require importing an external library.

| type                                                                                             | package                                              |
|:-------------------------------------------------------------------------------------------------|:-----------------------------------------------------|
| [`METISND`](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#CliqueTrees.METISND)     | [Metis.jl](https://github.com/JuliaSparse/Metis.jl)  |
| [`KaHyParND`](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#CliqueTrees.KaHyParND) | [KayHyPar.jl](https://github.com/kahypar/KaHyPar.jl) |

The optimizer is implemented using the tree decomposition library
[CliqueTrees.jl](https://github.com/AlgebraicJulia/CliqueTrees.jl).

# Arguments

- `dis`: [graph partitioning algorithm](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#CliqueTrees.DissectionAlgorithm)
- `algs`: tuple of [elimination algorithms](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#Elimination-Algorithms).
- `level`: maximum level
- `width`: minimum width
- `imbalances`: imbalance parameters 
- `score`: a function to evaluate the quality of the contraction tree. Default is `ScoreFunction()`.
"""
@kwdef struct HyperND{D, A} <: CodeOptimizer
    dis::D = KaHyParND()
    algs::A = (MF(), AMF(), MMD())
    level::Int = 6
    width::Int = 120
    imbalances::StepRange{Int, Int} = 130:1:130
    score::ScoreFunction = ScoreFunction()
end

function optimize_hyper_nd(optimizer::HyperND, code, size_dict)
    dis = optimizer.dis
    algs = optimizer.algs
    level = optimizer.level
    width = optimizer.width
    imbalances = optimizer.imbalances
    score = optimizer.score

    minscore = typemax(Float64)
    local mincode

    for imbalance in imbalances
        curalg = SafeRules(ND(BestWidth(algs), dis; level, width, imbalance))
        curoptimizer = Treewidth(; alg=curalg)
        curcode = _optimize_code(code, size_dict, curoptimizer)
        curtc, cursc, currw = __timespacereadwrite_complexity(curcode, size_dict)

        if score(curtc, cursc, currw) < minscore
            minscore, mincode = score(curtc, cursc, currw), curcode
        end
    end

    return mincode
end

function Base.show(io::IO, ::MIME"text/plain", optimizer::HyperND{D, A}) where {D, A}
    println(io, "HyperND{$D, $A}:")
    show(IOContext(io, :indent => 4), "text/plain", optimizer.dis)

    for alg in optimizer.algs
        show(IOContext(io, :indent => 4), "text/plain", alg)
    end

    println(io, "    level: $(optimizer.level)")
    println(io, "    width: $(optimizer.width)")
    println(io, "    imbalances: $(optimizer.imbalances)")
    println(io, "    target: $(optimizer.target)")
    return
end
