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
using CliqueTrees: cliquetree, cliquetree!, separator, residual, CliqueTree, EliminationAlgorithm, MMW, BFS, MCS, LexBFS, RCMMD, RCMGL, MCSM, LexM, AMF, MF, MMD, MF, BT, SafeRules, KaHyParND, METISND, ND, BestWidth, ConnectedComponents

# interfaces
export simplify_code, optimize_code, slice_code, optimize_permute, label_elimination_order, uniformsize, ScoreFunction

# optimizers
export CodeOptimizer, KaHyParBipartite, GreedyMethod, TreeSA, SABipartite, Treewidth, ExactTreewidth, HyperND, PathSA

# slicers
export CodeSlicer, TreeSASlicer

# decomposition types
export AbstractDecompositionType, TreeDecomp, PathDecomp

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
include("sabipartite.jl")
include("kahypar.jl")

# local search method
include("treesa.jl")
include("treesaslicer.jl")

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
@deprecate optimize_code(code::Union{EinCode, NestedEinsum}, size_dict::Dict, optimizer::CodeOptimizer, simplifier, permute=true) optimize_code(code, size_dict, optimizer; simplifier=simplifier, permute=permute)

end
