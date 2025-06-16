"""
    HyperND(;
        dis = KaHyParND(),
        algs = (MF(), MMD()),
        level = 6,
        width = 120,
        imbalances = 130:130,
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

"""
@kwdef struct HyperND{D, A} <: CodeOptimizer
    dis::D = KaHyParND()
    algs::A = (MF(), MMD())
    level::Int = 6
    width::Int = 120
    imbalances::StepRange{Int, Int} = 130:1:130
end

function optimize_hyper_nd(optimizer::HyperND, code, size_dict)
    dis = optimizer.dis
    algs = optimizer.algs
    level = optimizer.level
    width = optimizer.width
    imbalances = optimizer.imbalances

    mintc = minsc = minrw = typemax(Float64)
    mincode = nothing

    for alg in algs, imbalance in imbalances
        curalg = SafeRules(ND(alg, dis; imbalance))
        curoptimizer = Treewidth(; alg=curalg)
        curcode = _optimize_code(code, size_dict, curoptimizer)
        curtc, cursc, currw = __timespacereadwrite_complexity(curcode, size_dict)

        if minsc > cursc || (minsc == cursc && mintc > curtc) || (minsc == cursc && mintc == curtc && minrw > currw)
            mintc, minsc, minrw, mincode = curtc, cursc, currw, curcode
        end
    end

    return mincode
end

function Base.show(io::IO, ::MIME"text/plain", optimizer::HyperND{D, A}) where {D, A}
<<<<<<< HEAD
<<<<<<< HEAD
    println(io, "HyperND{$D, $A}:")
=======
    println(io, "HyPar{$D, $A}:")
>>>>>>> e48f22a (Add optimizer `HyperND`.)
=======
    println(io, "HyperND{$D, $A}:")
>>>>>>> 059085c (Fix printing.)
    show(IOContext(io, :indent => 4), "text/plain", optimizer.dis)

    for alg in optimizer.algs
        show(IOContext(io, :indent => 4), "text/plain", alg)
    end

    println(io, "    level: $(optimizer.level)")
    println(io, "    width: $(optimizer.width)")
    println(io, "    imbalances: $(optimizer.imbalances)")
    return
end
