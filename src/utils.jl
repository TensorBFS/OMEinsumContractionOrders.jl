function log2sumexp2(s)
    ms = maximum(s)
    return log2(sum(x->exp2(x - ms), s)) + ms
end

function _log2_size_dict(size_dict::Dict{L, T2}) where {L, T2}
    log2_size_dict = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_size_dict[k] = log2(v)
    end
    return log2_size_dict
end

# Given a vertex-edge incidence matrix `adj` and a subset of vertices `part`,
# group `part` into connected components (vertices sharing an edge are adjacent).
function _connected_components(adj, part::AbstractVector{T}) where T
    A = adj[part,:]
    A = A * A'   # connectivility matrix
    n = length(part)
    visit_mask = zeros(Bool, n)
    groups = Vector{T}[]
    while !all(visit_mask)
        newset = Int[]
        push_connected!(newset, visit_mask, A, findfirst(==(false), visit_mask))
        push!(groups, getindex.(Ref(part), newset))
    end
    return groups
end

function push_connected!(set, visit_mask, adj, i)
    visit_mask[i] = true
    push!(set, i)
    for v = 1:size(adj, 2)
        if !visit_mask[v] && !iszero(adj[i,v])
            push_connected!(set, visit_mask, adj, v)
        end
    end
end

"""
    convert_label(ne::NestedEinsum, labelmap::Dict{T1,T2}) where {T1,T2}

Convert the labels of a `NestedEinsum` object to new labels.
`labelmap` is a dictionary that maps the old labels to the new labels.
"""
function convert_label(ne::NestedEinsum, labelmap::Dict{T1,T2}) where {T1,T2}
    isleaf(ne) && return NestedEinsum{T2}(ne.tensorindex)
    eins = EinCode([getindex.(Ref(labelmap), ix) for ix in ne.eins.ixs], getindex.(Ref(labelmap), ne.eins.iy))
    NestedEinsum([convert_label(arg, labelmap) for arg in ne.args], eins)
end