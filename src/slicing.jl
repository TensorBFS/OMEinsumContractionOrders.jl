export Slicing, Slicer, SlicedEinsum

struct Slicer
    log2_sizes::Vector{Float64}   # the size dict after slicing
    legs::Dict{Int,Float64}   # sliced leg and its original size
    max_size::Int       # maximum number of sliced legs
end
Slicer(log2_sizes::AbstractVector{Float64}, max_size::Int) = Slicer(collect(log2_sizes), Dict{Int,Float64}(), max_size)
Base.length(s::Slicer) = length(s.legs)
function Base.replace!(slicer::Slicer, pair::Pair)
    worst, best = pair
    @assert haskey(slicer.legs, worst)
    @assert !haskey(slicer.legs, best)
    slicer.log2_sizes[worst] = slicer.legs[worst]       # restore worst size
    slicer.legs[best] = slicer.log2_sizes[best]  # add best to legs
    slicer.log2_sizes[best] = 0.0
    delete!(slicer.legs, worst)                  # remove worst from legs
    return slicer
end

function Base.push!(slicer, best)
    @assert length(slicer) < slicer.max_size
    @assert !haskey(slicer.legs, best)
    slicer.legs[best] = slicer.log2_sizes[best]  # add best to legs
    slicer.log2_sizes[best] = 0.0
    return slicer
end

struct Slicing{LT}
    legs::Vector{LT}   # sliced leg and its original size
end

Slicing(s::Slicer, inverse_map) = Slicing([inverse_map[l] for (l, s) in s.legs])
Base.length(s::Slicing) = length(s.legs)

struct SlicedEinsum{LT, Ein} <: AbstractEinsum
    slicing::Slicing{LT}
    eins::Ein
end

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
    n = length(se.slicing)

    size_dict_sliced = copy(size_dict)
    sliced_sizes = Vector{Int}(undef, n)
    sliced_labels = Vector{LT}(undef, n)
    for i = 1:n
        l = se.slicing.legs[i]
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
    return view(x,slices...)
end
function fill_slice!(x, ix, chunk, slicemap::Dict)
    if ndims(x) == 0
        x .+= chunk  # to avoid CUDA `getindex!`.
    else
        slices = map(l->haskey(slicemap, l) ? slicemap[l] : Colon(), ix)
        view(x, slices...) .+= chunk
    end
end

function (se::SlicedEinsum{LT,ET})(@nospecialize(xs::AbstractArray...); size_info = nothing) where {LT, ET}
    length(se.slicing) == 0 && return se.eins(xs...; size_info=size_info)
    size_dict = size_info===nothing ? Dict{OMEinsum.labeltype(se),Int}() : copy(size_info)
    OMEinsum.get_size_dict!(se, xs, size_dict)

    it = SliceIterator(se, size_dict)
    res = OMEinsum.get_output_array(xs, getindex.(Ref(size_dict), it.iyv))
    eins_sliced = drop_slicedim(se.eins, se.slicing)
    for slicemap in it
        xsi = ntuple(i->take_slice(xs[i], it.ixsv[i], slicemap), length(xs))
        resi = einsum(eins_sliced, xsi, it.size_dict_sliced)
        fill_slice!(res, it.iyv, resi, slicemap)
    end
    return res
end

function drop_slicedim(ne::NestedEinsum, slicing::Slicing)
    isleaf(ne) && return ne
    ixs = map(ix->filter(∉(slicing.legs), ix), getixsv(ne.eins))
    iy = filter(∉(slicing.legs), getiyv(ne.eins))
    NestedEinsum(map(arg->drop_slicedim(arg, slicing), ne.args), similar_eincode(ne.eins, ixs, iy))
end
similar_eincode(::DynamicEinCode, ixs, iy) = DynamicEinCode(ixs, iy)
similar_eincode(::StaticEinCode, ixs, iy) = StaticEinCode{(Tuple.(ixs)...,), (iy...,)}()

# forward some interfaces
function OMEinsum.timespace_complexity(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing.legs
        size_dict_sliced[l] = 1
    end
    tc, sc = timespace_complexity(code.eins, size_dict_sliced)
    tc + sum(log2.(getindex.(Ref(size_dict), code.slicing.legs))), sc
end
function OMEinsum.flop(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing.legs
        size_dict_sliced[l] = 1
    end
    fl = flop(code.eins, size_dict_sliced)
    fl * prod(getindex.(Ref(size_dict), code.slicing.legs))
end

OMEinsum.flatten(se::SlicedEinsum) = OMEinsum.flatten(se.eins)
OMEinsum.labeltype(::SlicedEinsum{LT}) where LT = LT
OMEinsum.get_size_dict!(se::SlicedEinsum, xs, size_info::Dict) = OMEinsum.get_size_dict!(se.eins, xs, size_info)
OMEinsum.getixsv(se::SlicedEinsum) = OMEinsum.getixsv(se.eins)
OMEinsum.getiyv(se::SlicedEinsum) = OMEinsum.getiyv(se.eins)
OMEinsum.label_elimination_order(code::SlicedEinsum) = OMEinsum.label_elimination_order(code.eins)
