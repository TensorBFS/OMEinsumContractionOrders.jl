module OMEinsumContractionOrders

using JSON
using SparseArrays
using StatsBase
using Base: RefValue
using Base.Threads
using AbstractTrees
using TreeWidthSolver
using TreeWidthSolver.Graphs
using DataStructures: PriorityQueue, enqueue!, dequeue!, peek, dequeue_pair!
import CliqueTrees
using CliqueTrees: cliquetree, residual, EliminationAlgorithm, MMW, BFS, MCS, LexBFS, RCMMD, RCMGL, MCSM, LexM, AMF, MF, MMD, MF, BT, SafeRules, KaHyParND, METISND, ND

# interfaces
export simplify_code, optimize_code, optimize_permute, label_elimination_order, uniformsize

# optimizers
export CodeOptimizer, KaHyParBipartite, GreedyMethod, TreeSA, SABipartite, Treewidth, ExactTreewidth, HyperND

# preprocessing
export CodeSimplifier, MergeGreedy, MergeVectors

# time space complexity
export peak_memory, flop, contraction_complexity

# Note: writejson, readjson are not exported to avoid namespace conflict

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

# nested dissection method
include("hypernd.jl")

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
