module OMEinsumContractionOrders

using JSON
using SparseArrays
using StatsBase
using Base: RefValue
using Base.Threads
using AbstractTrees
using TreeWidthSolver
using TreeWidthSolver.Graphs

export CodeOptimizer, CodeSimplifier,
    KaHyParBipartite, GreedyMethod, TreeSA, SABipartite, ExactTreewidth,
    MergeGreedy, MergeVectors,
    uniformsize,
    simplify_code, optimize_code, optimize_permute,
    # time space complexity
    peak_memory, flop, contraction_complexity,
    label_elimination_order
    # writejson, readjson are not exported to avoid namespace conflict

# visiualization tools provided by extension `LuxorTensorPlot`
export viz_eins, viz_contraction

include("Core.jl")
include("utils.jl")

# greedy method
include("incidencelist.jl")
include("greedy.jl")

# bipartition based methods
include("sa.jl")
include("kahypar.jl")

# local search method
include("treesa.jl")

# tree width method
include("treewidth.jl")

# simplification passes
include("simplify.jl")

# interfaces
include("complexity.jl")
include("interfaces.jl")

# saveload
include("json.jl")

# extension for visiualization
include("visualization.jl")

@deprecate timespacereadwrite_complexity(code, size_dict::Dict) (contraction_complexity(code, size_dict)...,)
@deprecate timespace_complexity(code, size_dict::Dict) (contraction_complexity(code, size_dict)...,)[1:2]

end
