export simplify_code, optimize_code, GreedyMethod, KaHyParBipartite, SABipartite, TreeSA, MergeGreedy, MergeVectors, uniformsize
using .ContractionOrderAlgorithms: CodeOptimizer, KaHyParBipartite, GreedyMethod, TreeSA, SABipartite, MergeGreedy, MergeVectors, uniformsize, simplify_code, optimize_code
import .ContractionOrderAlgorithms
export peak_memory

abstract type CodeSimplifier end

function ContractionOrderAlgorithms.optimize_code(code::AbstractEinsum, size_dict::Dict, optimizer::CodeOptimizer, simplifier=nothing, permute::Bool=true)
    optimize_code(rawcode(code), size_dict, optimizer, simplifier, permute)
end
ContractionOrderAlgorithms.simplify_code(code::AbstractEinsum, size_dict, simplifier::CodeSimplifier) = simplify_code(rawcode(code), size_dict, simplifier)
ContractionOrderAlgorithms.peak_memory(code::AbstractEinsum, size_dict::Dict) = peak_memory(rawcode(code), size_dict)

# TODO: figure out where is this function
ContractionOrderAlgorithms.uniformsize(code::AbstractEinsum, size) = Dict([l=>size for l in uniquelabels(code)])