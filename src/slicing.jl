export Slicing, Slicer, SlicedEinsum

struct Slicer
    log2_sizes::Vector{Float64}   # the size dict after slicing
    legs::Dict{Int,Float64}   # sliced leg and its original size
    max_size::Int       # maximum number of sliced legs
    fixed_slices::Vector{Int}     # number of fixed legs
end
function Slicer(log2_sizes::AbstractVector{Float64}, max_size::Int, fixed_slices::AbstractVector)
    slicer = Slicer(collect(log2_sizes), Dict{Int,Float64}(), max_size, collect(Int,fixed_slices))
    for l in fixed_slices
        push!(slicer, l)
    end
    return slicer
end
Base.length(s::Slicer) = length(s.legs)
function Base.replace!(slicer::Slicer, pair::Pair)
    worst, best = pair
    @assert worst ∉ slicer.fixed_slices
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
Base.:(==)(se::Slicing, se2::Slicing) = se.legs == se2.legs

function Slicing(s::Slicer, inverse_map::Dict{Int,LT}) where LT
    # we want to keep the order of input fixed slices!
    Slicing(LT[[inverse_map[l] for l in s.fixed_slices]..., [inverse_map[l] for (l, sz) in s.legs if l ∉ s.fixed_slices]...])
end
Base.length(s::Slicing) = length(s.legs)

struct SlicedEinsum{LT, Ein} <: AbstractEinsum
    slicing::Slicing{LT}
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
    return SliceIterator(ixsv, iyv, se.slicing.legs, size_dict)
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
end

function (se::SlicedEinsum{LT,ET})(@nospecialize(xs::AbstractArray...); size_info = nothing, kwargs...) where {LT, ET}
    length(se.slicing) == 0 && return se.eins(xs...; size_info=size_info, kwargs...)
    size_dict = size_info===nothing ? Dict{OMEinsum.labeltype(se),Int}() : copy(size_info)
    OMEinsum.get_size_dict!(se, xs, size_dict)

    it = SliceIterator(se, size_dict)
    res = OMEinsum.get_output_array(xs, getindex.(Ref(size_dict), it.iyv))
    eins_sliced = drop_slicedim(se.eins, se.slicing)
    for (k, slicemap) in enumerate(it)
        @debug "computing slice $k/$(length(it))"
        xsi = ntuple(i->take_slice(xs[i], it.ixsv[i], slicemap), length(xs))
        resi = einsum(eins_sliced, xsi, it.size_dict_sliced; kwargs...)
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
function OMEinsum.timespacereadwrite_complexity(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing.legs
        size_dict_sliced[l] = 1
    end
    tc, sc, rw = timespacereadwrite_complexity(code.eins, size_dict_sliced)
    sliceoverhead = sum(log2.(getindex.(Ref(size_dict), code.slicing.legs)))
    tc + sliceoverhead, sc, rw+sliceoverhead
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
