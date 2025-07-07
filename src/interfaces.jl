"""
    optimize_code(eincode, size_dict, optimizer = GreedyMethod(), simplifier=nothing, permute=true) -> optimized_eincode

Optimize the einsum contraction code and reduce the time/space complexity of tensor network contraction.
Returns a `NestedEinsum` instance. Input arguments are

# Arguments
- `eincode` is an einsum contraction code instance, one of `DynamicEinCode`, `StaticEinCode` or `NestedEinsum`.
- `size` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
- `optimizer` is a `CodeOptimizer` instance, should be one of `GreedyMethod`, `Treewidth`, `KaHyParBipartite`, `SABipartite` or `TreeSA`. Check their docstrings for details.
- `simplifier` is one of `MergeVectors` or `MergeGreedy`.
- `permute` is a boolean flag to indicate whether to optimize the permutation of the contraction order.

# Examples

```jldoctest
julia> using OMEinsum

julia> code = ein"ij, jk, kl, il->"
ij, jk, kl, il -> 

julia> optimize_code(code, uniformsize(code, 2), TreeSA());
```
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

function _optimize_code(code, size_dict, optimizer::KaHyParBipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::GreedyMethod)
    optimize_greedy(code, size_dict; α = optimizer.α, temperature = optimizer.temperature, nrepeat=optimizer.nrepeat)
end
function _optimize_code(code, size_dict, optimizer::Treewidth)
    optimize_treewidth(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::SABipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::TreeSA)
    optimize_tree(code, size_dict; sc_target=optimizer.sc_target, βs=optimizer.βs,
        ntrials=optimizer.ntrials, niters=optimizer.niters,
        sc_weight=optimizer.sc_weight, rw_weight=optimizer.rw_weight, initializer=optimizer.initializer,
        greedy_method=optimizer.greedy_config)
end
function _optimize_code(code, size_dict, optimizer::HyperND)
    optimize_hyper_nd(optimizer, code, size_dict)
end

"""
    slice_code(code, size_dict, slicer) -> sliced_code

Slice the einsum contraction code to reduce the space complexity, returns a `SlicedEinsum` instance.

# Arguments
- `code` is a `NestedEinsum` instance.
- `size_dict` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
- `slicer` is a `CodeSlicer` instance, currently only `TreeSASlicer` is supported.
"""
function slice_code(code::NestedEinsum, size_dict::Dict, slicer::CodeSlicer)
    sliced_code = _slice_code(code, size_dict, slicer)
    return sliced_code
end

function _slice_code(code, size_dict, slicer::TreeSASlicer)
    slice_tree(code, size_dict; sc_target=slicer.sc_target, βs=slicer.βs,
        ntrials=slicer.ntrials, niters=slicer.niters,
        sc_weight=slicer.sc_weight, rw_weight=slicer.rw_weight, optimization_ratio=slicer.optimization_ratio)
end