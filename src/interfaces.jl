export simplify_code, optimize_code, GreedyMethod, KaHyParBipartite, SABipartite, TreeSA, MergeGreedy, MergeVectors, uniformsize

abstract type CodeSimplifier end
Base.@kwdef struct MergeGreedy <: CodeSimplifier
    threshhold::Float64=-1e-12
end
struct MergeVectors <: CodeSimplifier end

"""
    optimize_code(eincode, size_dict, optimizer = GreedyMethod(), simplifier=nothing)

Optimize the einsum contraction code and reduce the time/space complexity of tensor network contraction.
Returns a `NestedEinsum` instance. Input arguments are

* `eincode` is an einsum contraction code instance, one of `DynamicEinCode`, `StaticEinCode` or `NestedEinsum`.
* `size` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
* `optimizer` is a `CodeOptimizer` instance, should be one of `GreedyMethod`, `KaHyParBipartite`, `SABipartite` or `TreeSA`. Check their docstrings for details.
* `simplifier` is one of `MergeVectors` or `MergeGreedy`.
"""
function optimize_code(code::Union{EinCode, NestedEinsum}, size_dict::Dict, optimizer::CodeOptimizer, simplifier=nothing)
    if simplifier === nothing
        return _optimize_code(code, size_dict, optimizer)
    else
        simpl, code = simplify_code(code, size_dict, simplifier)
        optcode = _optimize_code(code, size_dict, optimizer)
        return embed_simplifier(optcode, simpl)
    end
end

simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, ::MergeVectors) = merge_vectors(code)
simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, method::MergeGreedy) = merge_greedy(code, size_dict; threshhold=method.threshhold)

function _optimize_code(code, size_dict, optimizer::KaHyParBipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::GreedyMethod)
    optimize_greedy(code, size_dict; method=optimizer.method, nrepeat=optimizer.nrepeat)
end
function _optimize_code(code, size_dict, optimizer::SABipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::TreeSA)
    optimize_tree(code, size_dict; sc_target=optimizer.sc_target, βs=optimizer.βs,
        ntrials=optimizer.ntrials, niters=optimizer.niters, nslices=optimizer.nslices,
        sc_weight=optimizer.sc_weight, rw_weight=optimizer.rw_weight, initializer=optimizer.initializer,
        greedy_method=optimizer.greedy_config.method, greedy_nrepeat=optimizer.greedy_config.nrepeat)
end

uniformsize(code::AbstractEinsum, size) = Dict([l=>size for l in uniquelabels(code)])
