module OMEinsumContractionOrders

using Requires
using SparseArrays, Suppressor
using OMEinsum
using OMEinsum: NestedEinsum, getixsv, getiyv, DynamicEinCode, StaticEinCode, isleaf
export MinSpaceDiff, MinSpaceOut
export peak_memory

import .ContractionOrderAlgorithms
using .ContractionOrderAlgorithms: CodeOptimizer, CodeSimplifier,
    KaHyParBipartite, GreedyMethod, TreeSA, SABipartite,
    MinSpaceDiff, MinSpaceOut,
    MergeGreedy, MergeVectors,
    uniformsize,
    simplify_code, optimize_code,
    peak_memory

export CodeOptimizer, CodeSimplifier,
    KaHyParBipartite, GreedyMethod, TreeSA, SABipartite,
    MinSpaceDiff, MinSpaceOut,
    MergeGreedy, MergeVectors,
    uniformsize,
    simplify_code, optimize_code,
    peak_memory, timespace_complexity, timespacereadwrite_complexity

using JSON

include("algorithms/algorithms.jl")
include("omeinsum.jl")
include("interfaces.jl")
include("json.jl")

end
