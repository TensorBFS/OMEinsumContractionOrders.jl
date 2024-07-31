################### The data types in OMEinsum ###################
abstract type AbstractEinsum end

struct EinCode{LT} <: AbstractEinsum
    ixs::Vector{Vector{LT}}
    iy::Vector{LT}
end
getixsv(rc::EinCode) = rc.ixs
getiyv(rc::EinCode) = rc.iy
Base.:(==)(a::EinCode, b::EinCode) = a.ixs == b.ixs && a.iy == b.iy

struct NestedEinsum{LT} <: AbstractEinsum
    args::Vector{NestedEinsum}
    tensorindex::Int  # -1 if not leaf
    eins::EinCode{LT}
    NestedEinsum(args::Vector{NestedEinsum{LT}}, eins::EinCode) where LT = new{LT}(args, -1, eins)
    NestedEinsum{LT}(arg::Int) where LT = new{LT}(NestedEinsum{LT}[], arg)
end
function Base.:(==)(a::NestedEinsum, b::NestedEinsum)
    return a.args == b.args && a.tensorindex == b.tensorindex && if isdefined(a, :eins)
        isdefined(b, :eins) && a.eins == b.eins
    else
        !isdefined(b, :eins)
    end
end
isleaf(ne::NestedEinsum) = ne.tensorindex != -1
function getixsv(ne::NestedEinsum{LT}) where LT
    d = collect_ixs!(ne, Dict{Int,Vector{LT}}())
    ks = sort!(collect(keys(d)))
    return @inbounds [d[i] for i in ks]
end
function collect_ixs!(ne::NestedEinsum, d::Dict{Int,Vector{LT}}) where LT
    @inbounds for i=1:length(ne.args)
        arg = ne.args[i]
        if isleaf(arg)
            d[arg.tensorindex] = getixsv(ne.eins)[i]
        else
            collect_ixs!(arg, d)
        end
    end
    return d
end
getiyv(ne::NestedEinsum) = getiyv(ne.eins)

struct SlicedEinsum{LT,ET<:Union{EinCode{LT},NestedEinsum{LT}}} <: AbstractEinsum
    slicing::Vector{LT}
    eins::ET
end
Base.:(==)(a::SlicedEinsum, b::SlicedEinsum) = a.slicing == b.slicing && a.eins == b.eins
getixsv(ne::SlicedEinsum) = getixsv(ne.eins)
getiyv(ne::SlicedEinsum) = getiyv(ne.eins)
uniquelabels(code::AbstractEinsum) = unique!(vcat(getixsv(code)..., getiyv(code)))
labeltype(code::AbstractEinsum) = eltype(getiyv(code))

# Better printing
struct LeafString
    str::String
end
function AbstractTrees.children(ne::NestedEinsum)
    [isleaf(item) ? LeafString(_join(getixsv(ne.eins)[k])) : item for (k,item) in enumerate(ne.args)]
end
function AbstractTrees.printnode(io::IO, x::NestedEinsum)
    isleaf(x) ? print(io, x.tensorindex) : print(io, x.eins)
end
AbstractTrees.printnode(io::IO, e::LeafString) = print(io, e.str)
function Base.show(io::IO, e::EinCode)
    s = join([_join(ix) for ix in getixsv(e)], ", ") * " -> " * _join(getiyv(e))
    print(io, s)
end
function Base.show(io::IO, e::NestedEinsum)
    print_tree(io, e)
end
Base.show(io::IO, ::MIME"text/plain", e::NestedEinsum) = show(io, e)
Base.show(io::IO, ::MIME"text/plain", e::EinCode) = show(io, e)
_join(ix) = isempty(ix) ? "" : join(ix, connector(eltype(ix)))
connector(::Type{Char}) = ""
connector(::Type{Int}) = "∘"
connector(::Type) = "-"

function is_binary_tree(code::NestedEinsum)
    if isleaf(code) return true end
    if length(code.args) > 2 return false end
    return all(is_binary_tree, code.args)
end


# reformulate the nested einsum, removing a given tensor without change the space complexity
# consider only binary contraction tree with no openedges
function tree_reformulate(code::NestedEinsum, removed_tensor_id::Int)

    try @assert is_binary_tree(code) catch 
        error("The contraction tree is not binary") 
    end

    try @assert isempty(getiyv(code)) catch
        error("The contraction tree has open edges")
    end

    path = path_to_tensor(code, removed_tensor_id)

    right = popfirst!(path)
    left = right == 1 ? 2 : 1

    if isleaf(code.args[right])
        return NestedEinsum([code.args[left].args...], EinCode(getixsv(code.args[left].eins), getixsv(code.eins)[right]))
    else
        # update the ein code to make sure the root of the left part and the right part are the same
        left_code = code.args[left]
        right_code = NestedEinsum([code.args[right].args...], EinCode(getixsv(code.args[right].eins), getixsv(code.eins)[left]))
    end
    tree = _tree_reformulate!(left_code, right_code, path)

    return tree
end


function _tree_reformulate!(left_code::NestedEinsum{LT}, right_code::NestedEinsum{LT}, path::Vector{Int}) where{LT}
    if !isleaf(right_code)
        right = popfirst!(path)
        left = right == 1 ? 2 : 1
        if length(right_code.args) == 1
            # orign: left: a, right: b -> a
            # reformulated: left: a -> b, right: b
            new_eins = EinCode([getiyv(right_code.eins)], getixsv(right_code.eins)[1])
            left_code = NestedEinsum([left_code], new_eins)
            left_code = _tree_reformulate!(left_code, right_code.args[1], path)
        elseif length(right_code.args) == 2
            # origin: left: a, right: b, c -> a
            # reformulated: left: a, b -> c, right: c
            new_eins = EinCode([getiyv(right_code.eins), getixsv(right_code.eins)[left]], getixsv(right_code.eins)[right])
            left_code = NestedEinsum([left_code, right_code.args[left]], new_eins)
            left_code = _tree_reformulate!(left_code, right_code.args[right], path)
        else
            error("The contraction tree is not binary")
        end
    end
    return left_code
end

# find the path to a given tensor in a nested einsum
function path_to_tensor(code::NestedEinsum, index::Int)
    path = Vector{Int}()
    _find_root!(code, index, path)
    return path
end

function _find_root!(code::NestedEinsum, index::Int, path::Vector{Int})
    if isleaf(code) return code.tensorindex == index end

    for (i, arg) in enumerate(code.args)
        if _find_root!(arg, index, path)
            pushfirst!(path, i)
            return true
        end
    end
    return false
end

############### Simplifier and optimizer types #################
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

# code optimizer
abstract type CodeOptimizer end