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