#################### compute peak memory ###########################
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
    for l in code.slicing
        size_dict_sliced[l] = 1
    end
    return peak_memory(code.eins, size_dict_sliced) + _mem(getiyv(code.eins), size_dict)
end

###################### Time space complexity ###################
"""
    timespace_complexity(eincode, size_dict)

Returns the time and space complexity of the einsum contraction.
The time complexity is defined as `log2(number of element multiplication)`.
The space complexity is defined as `log2(size of the maximum intermediate tensor)`.
"""
function timespace_complexity(code, size_dict)
    tc,sc,rw = timespacereadwrite_complexity(code, size_dict)
    return tc, sc
end

"""
    timespacereadwrite_complexity(eincode, size_dict)

Returns the time, space and read-write complexity of the einsum contraction.
The time complexity is defined as `log2(number of element-wise multiplication)`.
The space complexity is defined as `log2(size of the maximum intermediate tensor)`.
The read-write complexity is defined as `log2(the number of read-write operations)`.
"""
function timespacereadwrite_complexity(ei::NestedEinsum, size_dict)
    log2_sizes = Dict([k=>log2(v) for (k,v) in size_dict])
    _timespacereadwrite_complexity(ei, log2_sizes)
end

function timespacereadwrite_complexity(ei::EinCode, size_dict)
    log2_sizes = Dict([k=>log2(v) for (k,v) in size_dict])
    _timespacereadwrite_complexity(getixsv(ei), getiyv(ei), log2_sizes)
end

function _timespacereadwrite_complexity(ei::NestedEinsum, log2_sizes::Dict{L,VT}) where {L,VT}
    isleaf(ei) && return (VT(-Inf), VT(-Inf), VT(-Inf))
    tcs = VT[]
    scs = VT[]
    rws = VT[]
    for arg in ei.args
        tc, sc, rw = _timespacereadwrite_complexity(arg, log2_sizes)
        push!(tcs, tc)
        push!(scs, sc)
        push!(rws, rw)
    end
    tc2, sc2, rw2 = _timespacereadwrite_complexity(getixsv(ei.eins), getiyv(ei.eins), log2_sizes)
    tc = log2sumexp2([tcs..., tc2])
    sc = max(reduce(max, scs), sc2)
    rw = log2sumexp2([rws..., rw2])
    return tc, sc, rw
end

function _timespacereadwrite_complexity(ixs::AbstractVector, iy::AbstractVector{T}, log2_sizes::Dict{L,VT}) where {T, L, VT}
    loop_inds = get_loop_inds(ixs, iy)
    tc = isempty(loop_inds) ? zero(VT) : sum(l->log2_sizes[l], loop_inds)
    sc = isempty(iy) ? zero(VT) : sum(l->log2_sizes[l], iy)
    rw = log2sumexp2([[isempty(ix) ? zero(VT) : sum(l->log2_sizes[l], ix) for ix in ixs]..., sc])
    return tc, sc, rw
end

function get_loop_inds(ixs::AbstractVector, iy::AbstractVector{LT}) where {LT}
    # remove redundant legs
    counts = Dict{LT,Int}()
    for ix in ixs
        for l in ix
            if haskey(counts, l)
                counts[l] += 1
            else
                counts[l] = 1
            end
        end
    end
    for l in iy
        if haskey(counts, l)
            counts[l] += 1
        else
            counts[l] = 1
        end
    end
    loop_inds = LT[]
    for ix in ixs
        for l in ix
            c = count(==(l), ix)
            if counts[l] > c && l ∉ loop_inds
                push!(loop_inds, l)
            end
        end
    end
    return loop_inds
end


"""
    flop(eincode, size_dict)

Returns the number of iterations, which is different with the true floating point operations (FLOP) by a factor of 2.
"""
function flop(ei::EinCode, size_dict::Dict{LT,VT}) where {LT,VT}
    loop_inds = uniquelabels(ei)
    return isempty(loop_inds) ? zero(VT) : prod(l->size_dict[l], loop_inds)
end

function flop(ei::NestedEinsum, size_dict::Dict{L,VT}) where {L,VT}
    isleaf(ei) && return zero(VT)
    return sum(ei.args) do arg
        flop(arg, size_dict)
    end + flop(ei.eins, size_dict)
end


############### Sliced methods   ##################
function timespacereadwrite_complexity(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing
        size_dict_sliced[l] = 1
    end
    tc, sc, rw = timespacereadwrite_complexity(code.eins, size_dict_sliced)
    sliceoverhead = sum(log2.(getindex.(Ref(size_dict), code.slicing)))
    tc + sliceoverhead, sc, rw+sliceoverhead
end
function flop(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing
        size_dict_sliced[l] = 1
    end
    fl = flop(code.eins, size_dict_sliced)
    fl * prod(getindex.(Ref(size_dict), code.slicing))
end

uniformsize(code::AbstractEinsum, size) = Dict([l=>size for l in uniquelabels(code)])

"""
    label_elimination_order(code)

Returns a vector of labels sorted by the order they are eliminated in the contraction tree.
The contraction tree is specified by `code`, which e.g. can be a `NestedEinsum` instance.
"""
label_elimination_order(code::NestedEinsum) = label_elimination_order!(code, labeltype(code)[])
function label_elimination_order!(code, eliminated_vertices)
    isleaf(code) && return eliminated_vertices
    for arg in code.args
        label_elimination_order!(arg, eliminated_vertices)
    end
    append!(eliminated_vertices, setdiff(vcat(getixsv(code.eins)...), getiyv(code.eins)))
    return eliminated_vertices
end
label_elimination_order(code::SlicedEinsum) = label_elimination_order(code.eins)

