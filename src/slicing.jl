export Slicer, SlicedEinsum

struct SlicedEinsum{LT, Ein} <: AbstractEinsum
    slicing::Vector{LT}
    eins::Ein
end
Base.:(==)(se::SlicedEinsum, se2::SlicedEinsum) = se.slicing == se2.slicing && se.eins == se2.eins

# Iterate over tensor network slices, its iterator interface returns `slicemap` as a Dict
# slice and fill tensors with
# * take_slice(x, label_of_x, slicemap)
# * fill_slice!(x, label_of_x, x_slice, slicemap)
struct SliceIterator{LT,IT<:CartesianIndices}
    ixsv::Vector{Vector{LT}}
    iyv::Vector{LT}
    sliced_labels::Vector{LT}
    indices::IT
    size_dict_sliced::Dict{LT,Int}
end

function SliceIterator(se::SlicedEinsum, size_dict::Dict{LT}) where LT
    iyv = OMEinsum.getiyv(se.eins.eins)
    ixsv = OMEinsum.getixsv(se.eins)
    return SliceIterator(ixsv, iyv, se.slicing, size_dict)
end

function SliceIterator(ixsv, iyv, legs, size_dict::Dict{LT}) where LT
    n = length(legs)
    size_dict_sliced = copy(size_dict)
    sliced_sizes = Vector{Int}(undef, n)
    sliced_labels = Vector{LT}(undef, n)
    for i = 1:n
        l = legs[i]
        sliced_sizes[i] = size_dict[l]
        sliced_labels[i] = l
        size_dict_sliced[l] = 1
    end
    indices = CartesianIndices((sliced_sizes...,))
    SliceIterator(ixsv, iyv, sliced_labels, indices, size_dict_sliced)
end
Base.length(si::SliceIterator) = length(si.indices)
Base.eltype(::Type{SliceIterator{LT,IT}}) where {LT,IT} = Dict{LT,Int}

# returns `slicemap` as a Dict
function Base.iterate(si::SliceIterator)
    ci, cistate = iterate(si.indices)
    slicemap = Dict(zip(si.sliced_labels, ones(Int,length(si.sliced_labels))))
    slicemap, (1,(ci,cistate),slicemap)
end
function Base.iterate(si::SliceIterator, state)
    i, (ci,cistate), slicemap = state
    if i >= length(si.indices)
        return nothing  # NOTE: ci is same as cistate
    else
        ci, cistate = iterate(si.indices, cistate)
        for (l, v) in zip(si.sliced_labels, ci.I)
            slicemap[l] = v
        end
        return slicemap, (i+1, (ci,cistate), slicemap)
    end
end
function Base.getindex(si::SliceIterator, indices...)
    ci = si.indices[indices...]
    slicemap = Dict(zip(si.sliced_labels, ci.I))
    return slicemap
end

function take_slice(x, ix, slicemap::Dict)
    slices = map(l->haskey(slicemap, l) ? slicemap[l] : Colon(), ix)
    if all(x->x isa Integer, slices)
        return copy(view(x,slices...))
    else
        return x[slices...]
    end
end
function fill_slice!(x, ix, chunk, slicemap::Dict)
    if ndims(x) == 0
        x .+= chunk  # to avoid CUDA `getindex!`.
    else
        slices = map(l->haskey(slicemap, l) ? slicemap[l] : Colon(), ix)
        view(x, slices...) .+= chunk
    end
    return x
end

function (se::SlicedEinsum{LT,ET})(@nospecialize(xs::AbstractArray...); size_info = nothing, kwargs...) where {LT, ET}
    # get size
    size_dict = size_info===nothing ? Dict{OMEinsum.labeltype(se),Int}() : copy(size_info)
    OMEinsum.get_size_dict!(se, xs, size_dict)
    # compute
    return einsum(se, xs, size_dict; kwargs...)
end

function OMEinsum.einsum(se::SlicedEinsum, @nospecialize(xs::NTuple{N,AbstractArray} where N), size_dict::Dict; kwargs...)
    length(se.slicing) == 0 && return einsum(se.eins, xs, size_dict; kwargs...)
    it = SliceIterator(se, size_dict)
    res = OMEinsum.get_output_array(xs, getindex.(Ref(size_dict), it.iyv))
    eins_sliced = drop_slicedim(se.eins, se.slicing)
    for (k, slicemap) in enumerate(it)
        @debug "computing slice $k/$(length(it))"
        xsi = ntuple(i->take_slice(xs[i], it.ixsv[i], slicemap), length(xs))
        resi = einsum(eins_sliced, xsi, it.size_dict_sliced; kwargs...)
        res = fill_slice!(res, it.iyv, resi, slicemap)
    end
    return res
end

function drop_slicedim(ne::NestedEinsum, slices::Vector)
    isleaf(ne) && return ne
    ixs = map(ix->filter(∉(slices), ix), getixsv(ne.eins))
    iy = filter(∉(slices), getiyv(ne.eins))
    NestedEinsum(map(arg->drop_slicedim(arg, slices), ne.args), similar_eincode(ne.eins, ixs, iy))
end
similar_eincode(::DynamicEinCode, ixs, iy) = DynamicEinCode(ixs, iy)
similar_eincode(::StaticEinCode, ixs, iy) = StaticEinCode{(Tuple.(ixs)...,), (iy...,)}()

# forward some interfaces
function OMEinsum.timespacereadwrite_complexity(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing
        size_dict_sliced[l] = 1
    end
    tc, sc, rw = timespacereadwrite_complexity(code.eins, size_dict_sliced)
    sliceoverhead = sum(log2.(getindex.(Ref(size_dict), code.slicing)))
    tc + sliceoverhead, sc, rw+sliceoverhead
end
function OMEinsum.flop(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing
        size_dict_sliced[l] = 1
    end
    fl = flop(code.eins, size_dict_sliced)
    fl * prod(getindex.(Ref(size_dict), code.slicing))
end

OMEinsum.flatten(se::SlicedEinsum) = OMEinsum.flatten(se.eins)
OMEinsum.labeltype(::SlicedEinsum{LT}) where LT = LT
OMEinsum.get_size_dict!(se::SlicedEinsum, xs, size_info::Dict) = OMEinsum.get_size_dict!(se.eins, xs, size_info)
OMEinsum.getixsv(se::SlicedEinsum) = OMEinsum.getixsv(se.eins)
OMEinsum.getiyv(se::SlicedEinsum) = OMEinsum.getiyv(se.eins)
OMEinsum.label_elimination_order(code::SlicedEinsum) = OMEinsum.label_elimination_order(code.eins)
