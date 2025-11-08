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