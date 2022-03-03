export simplify_code, optimize_code, GreedyMethod, KaHyParBipartite, SABipartite, TreeSA, MergeGreedy, MergeVectors, uniformsize

abstract type CodeSimplifier end

"""
    MergeGreedy <: CodeSimplifier
    MergeGreedy(; threshhold=-1e-12)

Contraction code simplifier (in order to reduce the time of calling optimizers) that
merges tensors greedily if the space complexity of merged tensors is reduced (difference smaller than the `threshhold`).
"""
Base.@kwdef struct MergeGreedy <: CodeSimplifier
    threshhold::Float64=-1e-12
end

"""
    MergeVectors <: CodeSimplifier
    MergeVectors()

Contraction code simplifier (in order to reduce the time of calling optimizers) that merges vectors to closest tensors.
"""
struct MergeVectors <: CodeSimplifier end

"""
    optimize_code(eincode, size_dict, optimizer = GreedyMethod(), simplifier=nothing, permute=true)

Optimize the einsum contraction code and reduce the time/space complexity of tensor network contraction.
Returns a `NestedEinsum` instance. Input arguments are

* `eincode` is an einsum contraction code instance, one of `DynamicEinCode`, `StaticEinCode` or `NestedEinsum`.
* `size` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
* `optimizer` is a `CodeOptimizer` instance, should be one of `GreedyMethod`, `KaHyParBipartite`, `SABipartite` or `TreeSA`. Check their docstrings for details.
* `simplifier` is one of `MergeVectors` or `MergeGreedy`.
* optimize the permutation if `permute` is true.
"""
function optimize_code(code::Union{EinCode, NestedEinsum}, size_dict::Dict, optimizer::CodeOptimizer, simplifier=nothing, permute::Bool=true)
    if simplifier === nothing
        optcode = _optimize_code(code, size_dict, optimizer)
    else
        simpl, code = simplify_code(code, size_dict, simplifier)
        optcode0 = _optimize_code(code, size_dict, optimizer)
        optcode = embed_simplifier(optcode0, simpl)
    end
    if permute
        optimize_permute(optcode, 0)
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

export peak_memory
"""
    peak_memory(code, size_dict::Dict)

Estimate peak memory usage in number of elements.
"""
function peak_memory(code::NestedEinsum, size_dict::Dict)
    ixs = getixsv(code.eins)
    iy = getiyv(code.eins)
    # `largest_size` is the largest size during contraction
    largest_size = 0
    # `tempsize` is the memory to store contraction results from previous branches
    tempsize = 0
    for (i, arg) in enumerate(code.args)
        if isleaf(arg)
            largest_size_i = _mem(ixs[i], size_dict) + tempsize
        else
            largest_size_i = peak_memory(arg, size_dict) + tempsize
        end
        tempsize += _mem(ixs[i], size_dict)
        largest_size = max(largest_size, largest_size_i)
    end
    # compare with currect contraction
    return max(largest_size, tempsize + _mem(iy, size_dict))
end
_mem(iy, size_dict::Dict{LT,VT}) where {LT,VT} = isempty(iy) ? zero(VT) : prod(l->size_dict[l], iy)

function peak_memory(code::EinCode, size_dict::Dict)
    ixs = getixsv(code)
    iy = getiyv(code)
    return sum(ix->_mem(ix, size_dict), ixs) + _mem(iy, size_dict)
end

function peak_memory(code::SlicedEinsum, size_dict::Dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing.legs
        size_dict_sliced[l] = 1
    end
    return peak_memory(code.eins, size_dict_sliced) + _mem(getiyv(code.eins), size_dict)
end