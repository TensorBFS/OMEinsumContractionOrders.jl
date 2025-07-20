"""
    optimize_code(eincode, size_dict, optimizer = GreedyMethod(); slicer=nothing, simplifier=nothing, permute=true) -> optimized_eincode

Optimize the einsum contraction code and reduce the time/space complexity of tensor network contraction.
Returns a `NestedEinsum` instance. Input arguments are

# Arguments
- `eincode` is an einsum contraction code instance, one of `DynamicEinCode`, `StaticEinCode` or `NestedEinsum`.
- `size` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
- `optimizer` is a `CodeOptimizer` instance, should be one of `GreedyMethod`, `Treewidth`, `KaHyParBipartite`, `SABipartite` or `TreeSA`. Check their docstrings for details.

# Keyword Arguments
- `slicer` is for slicing the contraction code to reduce the space complexity, default is nothing. Currently only [`TreeSASlicer`](@ref) is supported.
- `simplifier` is one of `MergeVectors` or `MergeGreedy`. Default is nothing.
- `permute` is a boolean flag to indicate whether to optimize the permutation of the contraction order.

# Examples

```jldoctest
julia> using OMEinsum

julia> code = ein"ij, jk, kl, il->"
ij, jk, kl, il -> 

julia> optimize_code(code, uniformsize(code, 2), TreeSA());
```
"""
function optimize_code(code::Union{EinCode, NestedEinsum}, size_dict::Dict, optimizer::CodeOptimizer; slicer=nothing, simplifier=nothing, permute::Bool=true)
    if simplifier === nothing
        optcode = _optimize_code(code, size_dict, optimizer)
    else
        simpl, code = simplify_code(code, size_dict, simplifier)
        optcode0 = _optimize_code(code, size_dict, optimizer)
        optcode = embed_simplifier(optcode0, simpl)
    end
    if slicer !== nothing
        optcode = slice_code(optcode, size_dict, slicer)
    end
    if permute
        optimize_permute(optcode, 0)
    end
end

function _optimize_code(code, size_dict, optimizer::KaHyParBipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::GreedyMethod)
    optimize_greedy(code, size_dict; α = optimizer.α, temperature = optimizer.temperature)
end
function _optimize_code(code, size_dict, optimizer::Treewidth)
    optimize_treewidth(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::SABipartite)
    recursive_bipartite_optimize(optimizer, code, size_dict)
end
function _optimize_code(code, size_dict, optimizer::TreeSA)
    optimize_tree(code, size_dict; βs=optimizer.βs,
        ntrials=optimizer.ntrials, niters=optimizer.niters,
        initializer=optimizer.initializer, score=optimizer.score)
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
function slice_code(code::NestedEinsum, size_dict, slicer::TreeSASlicer)
    slice_tree(code, size_dict; score=slicer.score, βs=slicer.βs,
        ntrials=slicer.ntrials, niters=slicer.niters,
        optimization_ratio=slicer.optimization_ratio)
end
