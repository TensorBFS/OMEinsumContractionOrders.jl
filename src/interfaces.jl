"""
    optimize_code(eincode, size_dict, optimizer = GreedyMethod(), simplifier=nothing, permute=true) -> optimized_eincode

Optimize the einsum contraction code and reduce the time/space complexity of tensor network contraction.
Returns a `NestedEinsum` instance. Input arguments are

* `eincode` is an einsum contraction code instance, one of `DynamicEinCode`, `StaticEinCode` or `NestedEinsum`.
* `size` is a dictionary of "edge label=>edge size" that contains the size information, one can use `uniformsize(eincode, 2)` to create a uniform size.
* `optimizer` is a `CodeOptimizer` instance, should be one of `GreedyMethod`, `Treewidth`, `KaHyParBipartite`, `SABipartite` or `TreeSA`. Check their docstrings for details.
* `simplifier` is one of `MergeVectors` or `MergeGreedy`.
* optimize the permutation if `permute` is true.

### Examples

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

simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, ::MergeVectors) = merge_vectors(code)
simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, method::MergeGreedy) = merge_greedy(code, size_dict; threshhold=method.threshhold)

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
        ntrials=optimizer.ntrials, niters=optimizer.niters, nslices=optimizer.nslices,
        sc_weight=optimizer.sc_weight, rw_weight=optimizer.rw_weight, initializer=optimizer.initializer,
        greedy_method=optimizer.greedy_config, fixed_slices=optimizer.fixed_slices)
end
