
############### Simplifier and optimizer types #################
"""
    CodeSimplifier

Abstract type for code simplifiers.
"""
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
    NetworkSimplifier{LT}

A network simplifier that contains a list of operations that can be applied to a tensor network to reduce the number of tensors.
It is generated from a proprocessor, such as [`MergeVectors`](@ref) or [`MergeGreedy`](@ref).

# Fields
- `operations`: a list of `NestedEinsum` objects.
"""
struct NetworkSimplifier{LT}
    operations::Vector{NestedEinsum{LT}}
end

"""
    simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, method::CodeSimplifier)

Simplify the contraction code by preprocessing the code with a simplifier.

# Arguments
- `code`: the contraction code to simplify.
- `size_dict`: the size dictionary of the contraction code.
- `method`: the simplifier to use, which can be [`MergeVectors`](@ref) or [`MergeGreedy`](@ref).

# Returns
- A tuple of `(NetworkSimplifier, newcode)`, where `newcode` is a new `EinCode` object.
"""
simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, ::MergeVectors) = merge_vectors(code)
simplify_code(code::Union{EinCode, NestedEinsum}, size_dict, method::MergeGreedy) = merge_greedy(code, size_dict; threshhold=method.threshhold)

function merge_vectors(code::EinCode{LT}) where LT
    ixs = getixsv(code)
    mask = trues(length(ixs))
    ops = [NestedEinsum{LT}(i) for i=1:length(ixs)]
    for i in 1:length(ixs)
        if length(ixs[i]) == 1
            for j in 1:length(ixs)
                if i!=j && mask[j] && ixs[i][1] ∈ ixs[j]  # merge i to j
                    mask[i] = false
                    ops[j] = NestedEinsum([ops[i], ops[j]],
                        EinCode([ixs[i], ixs[j]], ixs[j]))
                    break
                end
            end
        end
    end
    newcode = EinCode(ixs[mask], getiyv(code))
    return NetworkSimplifier(ops[mask]), newcode
end

function merge_greedy(code::EinCode{LT}, size_dict; threshhold=-1e-12) where LT
    ixs, iy = getixsv(code), getiyv(code)
    log2_edge_sizes = Dict{LT,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    n = nv(incidence_list)
    if n == 0
        return nothing
    elseif n == 1
        return collect(vertices(incidence_list))[1]
    end
    tree = Dict{Int,NestedEinsum}([v=>NestedEinsum{LT}(v) for v in vertices(incidence_list)])
    cost_values, cost_graph = evaluate_costs(1.0, incidence_list, log2_edge_sizes)
    while true
        if isempty(cost_values)
            return _buildsimplifier(tree, incidence_list)
        end
        pair, v = peek(cost_values)
        if v <= threshhold
            _, _, c = contract_pair!(incidence_list, pair..., log2_edge_sizes)
            tree[pair[1]] = NestedEinsum([tree[pair[1]], tree[pair[2]]], EinCode([c.first...], c.second))
            if nv(incidence_list) <= 1
                return _buildsimplifier(tree, incidence_list)
            end
            update_costs!(cost_values, cost_graph, pair..., 1.0, incidence_list, log2_edge_sizes)
        else
            return _buildsimplifier(tree, incidence_list)
        end
    end
end
function _buildsimplifier(tree, incidence_list)
    vertices = sort!(collect(keys(incidence_list.v2e)))
    ixs = [incidence_list.v2e[v] for v in vertices]
    iy = incidence_list.openedges
    NetworkSimplifier([tree[v] for v in vertices]), EinCode(ixs, iy)
end

"""
    embed_simplifier(code::NestedEinsum, simplifier::NetworkSimplifier)

Embed the simplifier into the contraction code. A typical workflow is:
(i) generate a simplifier with [`simplify_code`](@ref), (ii) then optimize the simplified code with [`optimize_code`](@ref)
and (iii) post-process the optimized code with [`embed_simplifier`](@ref) to produce correct contraction order for the original code.
This is automatically done in [`optimize_code`](@ref) given the `simplifier` argument is not `nothing`.

# Arguments
- `code`: the contraction code to embed the simplifier into.
- `simplifier`: the simplifier to embed, which is a [`NetworkSimplifier`](@ref) object.

# Returns
- A new `NestedEinsum` object.
"""
function embed_simplifier(code::NestedEinsum, simplifier::NetworkSimplifier)
    if isleaf(code)
        op = simplifier.operations[code.tensorindex]
        return op
    else
        return NestedEinsum(map(code.args) do arg
            embed_simplifier(arg, simplifier)
        end, code.eins)
    end
end
embed_simplifier(code::SlicedEinsum, simplifier::NetworkSimplifier) = SlicedEinsum(code.slicing, embed_simplifier(code.eins, simplifier))

optimize_permute(se::SlicedEinsum, level=0) = SlicedEinsum(se.slicing, se.eins isa EinCode ? se.eins : optimize_permute(se.eins, level))
function optimize_permute(ne::NestedEinsum{LT}, level=0) where LT
    if isleaf(ne)
        return ne
    else
        args = NestedEinsum{LT}[optimize_permute(arg, level+1) for arg in ne.args]
        ixs0 = getixsv(ne.eins)
        ixs = Vector{LT}[isleaf(x) ? ixs0[i] : getiyv(x.eins) for (i, x) in enumerate(args)]
        iy = level == 0 ? getiyv(ne.eins) : optimize_output_permute(ixs, getiyv(ne.eins))
        return NestedEinsum(args, EinCode(ixs, iy))
    end
end

function optimize_output_permute(ixs::AbstractVector{<:AbstractVector{LT}}, iy::AbstractVector{LT}) where LT
    if length(ixs) != 2
        return iy
    else
        iA, iB = ixs
        batchdim = LT[]
        outerA = LT[]
        outerB = LT[]
        bcastdim = LT[]
        for l in iy
            if l ∈ iA
                if l ∈ iB
                    push!(batchdim, l)
                else
                    push!(outerA, l)
                end
            else
                if l ∈ iB
                    push!(outerB, l)
                else
                    push!(bcastdim, l)
                end
            end
        end
        return vcat(
            sort!(outerA, by=l->findfirst(==(l), iA)),
            sort!(outerB, by=l->findfirst(==(l), iB)),
            sort!(batchdim, by=l->findfirst(==(l), iA)),
            bcastdim)
    end
end
