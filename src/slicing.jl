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

struct SlicedEinsum{LT, Ein}
    slicing::Slicing{LT}
    eins::Ein
end

function (se::SlicedEinsum{LT,ET})(@nospecialize(xs::AbstractArray...); size_info = nothing) where {LT, ET}
    length(se.slicing) == 0 && return se.eins(xs...; size_info=size_info)
    size_dict = size_info===nothing ? Dict{OMEinsum.labeltype(se),Int}() : copy(size_info)
    OMEinsum.get_size_dict!(se, xs, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in se.slicing.legs
        size_dict_sliced[l] = 1
    end

    # TODO: fix this interface, do not assume `eins` being NestedEinsum
    iy = OMEinsum.getiy(se.eins.eins)
    ixs = collect_ixs(se.eins)
    res = OMEinsum.get_output_array(xs, getindex.(Ref(size_dict), iy))

    sliced_sizes = Int[]
    sliced_labels = LT[]
    for l in se.slicing.legs
        push!(sliced_sizes, size_dict[l])
        push!(sliced_labels, l)
    end
    for ci in CartesianIndices((sliced_sizes...,))
        slicemap = Dict(zip(sliced_labels, ci.I))
        xsi = ntuple(i->take_slice(xs[i], ixs[i], slicemap), length(xs))
        resi = einsum(se.eins, xsi, size_dict_sliced)
        fill_slice!(res, iy, resi, slicemap)
    end
    return res
end

function take_slice(x, ix, slicemap::Dict)
    slices = map(l->haskey(slicemap, l) ? (slicemap[l]:slicemap[l]) : Colon(), ix)
    return x[slices...]
end
function fill_slice!(x, ix, chunk, slicemap::Dict)
    if ndims(x) == 0
        x[] += chunk[]
    else
        slices = map(l->haskey(slicemap, l) ? (slicemap[l]:slicemap[l]) : Colon(), ix)
        x[slices...] .+= chunk
    end
end

# forward some interfaces
function OMEinsum.timespace_complexity(code::SlicedEinsum, size_dict)
    size_dict_sliced = copy(size_dict)
    for l in code.slicing.legs
        size_dict_sliced[l] = 1
    end
    tc, sc = timespace_complexity(code.eins, size_dict_sliced)
    tc + sum(log2.(getindex.(Ref(size_dict), code.slicing.legs))), sc
end

OMEinsum.flatten(se::SlicedEinsum) = OMEinsum.flatten(se.eins)
OMEinsum.labeltype(::SlicedEinsum{LT}) where LT = LT
OMEinsum.get_size_dict!(se::SlicedEinsum, xs, size_info::Dict) = OMEinsum.get_size_dict!(se.eins, xs, size_info)
OMEinsum.getixsv(se::SlicedEinsum) = OMEinsum.getixsv(se.eins)
OMEinsum.getiyv(se::SlicedEinsum) = OMEinsum.getiyv(se.eins)
