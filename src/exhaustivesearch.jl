# Exact ("netcon") contraction-order optimizer.
#
# This is a port of the `optimaltree` routine from TensorOperations.jl
# (https://github.com/Jutho/TensorOperations.jl, MIT licensed), which implements
# the cost-capped breadth-first dynamic program of
#
#   Pfeifer, R.N.C., Haegeman, J., Verstraete, F., 2014.
#   Faster identification of optimal contraction sequences for tensor networks.
#   Phys. Rev. E 90, 033315. https://doi.org/10.1103/PhysRevE.90.033315
#
# The algorithm finds the binary contraction tree that minimizes the *total*
# contraction cost (the sum of FLOPs over all pairwise contractions), which is
# exactly the time complexity `tc`. Costs are accumulated in `Float64` to avoid
# integer overflow on large networks.
#
# Compared with the original, this version is adapted to OMEinsum's
# representation: it tracks each label's occurrence set explicitly (so it
# supports hyperedges and batch/diagonal output indices), seeds the cost cap
# per connected component from a greedy contraction, and builds the output
# `NestedEinsum` directly from the search state.

"""
    ExhaustiveSearch <: CodeOptimizer
    ExhaustiveSearch(; verbose = false)

Exact contraction-order optimizer based on the cost-capped breadth-first dynamic
program of Pfeifer, Haegeman and Verstraete (2014) (the "netcon" algorithm),
ported from [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl).

It returns the contraction tree that **globally minimizes the total contraction
cost** (the sum of FLOPs over all pairwise contractions, i.e. the time complexity
`tc`). This is complementary to [`ExactTreewidth`](@ref), which instead minimizes
the treewidth (the size of the largest intermediate tensor, i.e. the space
complexity `sc`).

The search is exponential in the number of tensors; it is intended for small to
moderate networks (a few tens of tensors). The cost cap of each connected
component is seeded with a greedy contraction, so the search typically finds the
optimum in a single cost-capped pass.

# Scope
Hyperedges (a label shared by more than two tensors) and batch/diagonal output
indices are supported. Two cases are rejected with an `ArgumentError`: a label
that appears more than once within a single tensor (a partial trace), and a
label that appears once but is not an output (a dangling summed index). Simplify
those away first (e.g. with a simplifier) or use another optimizer.

# Fields
- `verbose::Bool`: print progress of the search. Default `false`.
"""
Base.@kwdef struct ExhaustiveSearch <: CodeOptimizer
    verbose::Bool = false
end

############################ set representations ############################
# Subsets of tensors (DP keys), subsets of labels (open-index sets) and the
# per-label occurrence sets are all stored as bit sets, using the narrowest
# unsigned integer that fits, falling back to `BitVector`.
function _storeset(::Type{BitVector}, ints, maxint)
    set = falses(maxint)
    for i in ints
        set[i] = true
    end
    return set
end
function _storeset(::Type{T}, ints, maxint) where {T <: Unsigned}
    set = zero(T)
    u = one(T)
    for i in ints
        set |= (u << (i - 1))
    end
    return set
end

_intersect(s1::T, s2::T) where {T <: Unsigned} = s1 & s2
_intersect(s1::BitVector, s2::BitVector) = s1 .& s2
_union(s1::T, s2::T) where {T <: Unsigned} = s1 | s2
_union(s1::BitVector, s2::BitVector) = s1 .| s2
_setdiff(s1::T, s2::T) where {T <: Unsigned} = s1 & (~s2)
_setdiff(s1::BitVector, s2::BitVector) = s1 .& (.~s2)
_isemptyset(s::Unsigned) = iszero(s)
_isemptyset(s::BitVector) = !any(s)
_issubset(s1, s2) = _isemptyset(_setdiff(s1, s2))  # s1 ⊆ s2
_hasbit(s::Unsigned, i) = !iszero(s & (one(s) << (i - 1)))
_hasbit(s::BitVector, i) = s[i]

############################ cost model (Float64) ############################
# Product of the dimensions of every label in a set.
function _setcost(allcosts::Vector{Float64}, set::T) where {T <: Unsigned}
    cost = 1.0
    ind = set
    n = 1
    @inbounds while !iszero(ind)
        if isodd(ind)
            cost *= allcosts[n]
        end
        ind >>= 1
        n += 1
    end
    return cost
end
function _setcost(allcosts::Vector{Float64}, set::BitVector)
    cost = 1.0
    @inbounds for n in findall(set)
        cost *= allcosts[n]
    end
    return cost
end

# A loose upper bound on the cost of contracting a component, used only as a
# safe ceiling for the cost-cap loop (the real cap is seeded from greedy). The
# cost of contracting two tensors is the product of the dimensions of every
# label on either tensor, i.e. the FLOP count of that pairwise contraction.
function _maxcost(allcosts::Vector{Float64}, indexsets)
    length(indexsets) <= 1 && return 1.0
    maxcost = 0.0
    s1 = indexsets[1]
    for n in 2:length(indexsets)
        s2 = indexsets[n]
        maxcost += _setcost(allcosts, _union(s1, s2))
        s1 = _setdiff(_union(s1, s2), _intersect(s1, s2))
    end
    return maxcost
end

# Labels of `allopen` that become fully contracted when the tensor-subset `s` is
# formed: those not in the output whose every occurrence lies inside `s`.
function _closedmask(allopen::T, s::T, tensormask::Vector{T}, iymask::T) where {T <: Unsigned}
    closed = zero(T)
    ind = allopen
    n = 1
    u = one(T)
    @inbounds while !iszero(ind)
        if isodd(ind) && iszero(iymask & (u << (n - 1))) && _issubset(tensormask[n], s)
            closed |= (u << (n - 1))
        end
        ind >>= 1
        n += 1
    end
    return closed
end
function _closedmask(allopen::BitVector, s::BitVector, tensormask::Vector{BitVector}, iymask::BitVector)
    closed = falses(length(allopen))
    @inbounds for n in findall(allopen)
        if !iymask[n] && _issubset(tensormask[n], s)
            closed[n] = true
        end
    end
    return closed
end

############################ core dynamic program ############################
# Find the optimal contraction tree of one connected component, returning the
# tree as a nested `Any[left, right]` of integer tensor positions, the open
# (uncontracted) label set of the whole component, and the optimal total cost.
function _solve_component(
        ::Type{T}, component, singleopen::Vector{T}, tensormask::Vector{T}, iymask::T,
        numtensors::Int, allcosts::Vector{Float64}, initialcost::Float64; verbose::Bool = false
    ) where {T}
    componentsize = length(component)
    if componentsize == 1
        i = component[1]
        return (i, singleopen[i], 0.0)
    end

    costdict = [Dict{T, Float64}() for _ in 1:componentsize]
    # treedict stores only the bipartition (the left-child subset key); the tree
    # is reconstructed once at the end, avoiding a nested `Any[]` alloc per merge.
    treedict = [Dict{T, T}() for _ in 1:componentsize]
    indexdict = [Dict{T, T}() for _ in 1:componentsize]
    for i in component
        s = _storeset(T, [i], numtensors)
        costdict[1][s] = 0.0
        indexdict[1][s] = singleopen[i]
    end

    costfac = maximum(allcosts)
    maxcost = max(_maxcost(allcosts, @view(singleopen[component])), initialcost)
    currentcost = min(initialcost, maxcost)
    previouscost = 0.0
    while currentcost <= maxcost
        nextcost = maxcost
        for n in 2:componentsize
            verbose && println("Component $(component): subsets of size $n with cost cap $currentcost")
            for k in 1:div(n - 1, 2)
                for s1 in keys(costdict[k]), s2 in keys(costdict[n - k])
                    nextcost = _trymerge!(costdict, treedict, indexdict, n, k, s1, s2,
                        tensormask, iymask, allcosts, currentcost, previouscost, nextcost)
                end
            end
            if iseven(n) # k = n/2: pair distinct subsets of the same dict once
                k = div(n, 2)
                it = keys(costdict[k])
                st1 = iterate(it)
                while st1 !== nothing
                    s1, ns1 = st1
                    st2 = iterate(it, ns1)
                    while st2 !== nothing
                        s2, ns2 = st2
                        nextcost = _trymerge!(costdict, treedict, indexdict, n, k, s1, s2,
                            tensormask, iymask, allcosts, currentcost, previouscost, nextcost)
                        st2 = iterate(it, ns2)
                    end
                    st1 = iterate(it, ns1)
                end
            end
        end
        !isempty(costdict[componentsize]) && break
        previouscost = currentcost
        # grow the cap; guarantee progress even when costfac == 1 (all dims 1)
        currentcost = min(maxcost, max(nextcost * costfac, nextcost))
        currentcost <= previouscost && break
    end
    if isempty(costdict[componentsize])
        error("ExhaustiveSearch: cost cap $maxcost reached without a solution") # should be impossible
    end
    s = _storeset(T, component, numtensors)
    return (_reconstruct_tree(s, treedict), indexdict[componentsize][s], costdict[componentsize][s])
end

# Rebuild the nested `Any[left, right]`/Int tree from the stored bipartitions.
function _reconstruct_tree(s::T, treedict::Vector{Dict{T, T}}) where {T}
    sz = _setsize(s)
    sz == 1 && return _singleton_index(s)
    s1 = treedict[sz][s]
    s2 = _setdiff(s, s1)
    return Any[_reconstruct_tree(s1, treedict), _reconstruct_tree(s2, treedict)]
end

# Try to record the contraction of subsets `s1` (size k) and `s2` (size n-k).
# Returns the updated `nextcost` (the smallest over-cap cost seen this pass).
@inline function _trymerge!(costdict, treedict, indexdict, n, k, s1::T, s2::T,
        tensormask::Vector{T}, iymask::T, allcosts, currentcost, previouscost, nextcost) where {T}
    _isemptyset(_intersect(s1, s2)) || return nextcost  # subsets must be disjoint
    s = _union(s1, s2)
    existing = get(costdict[n], s, currentcost)
    existing > previouscost || return nextcost
    ind1 = indexdict[k][s1]
    ind2 = indexdict[n - k][s2]
    _isemptyset(_intersect(ind1, ind2)) && return nextcost  # skip outer products within a component
    allopen = _union(ind1, ind2)
    cost = costdict[k][s1] + costdict[n - k][s2] + _setcost(allcosts, allopen)
    if cost <= existing
        costdict[n][s] = cost
        indexdict[n][s] = _setdiff(allopen, _closedmask(allopen, s, tensormask, iymask))
        treedict[n][s] = s1  # store the bipartition; partner is s ∖ s1
    elseif currentcost < cost < nextcost
        nextcost = cost
    end
    return nextcost
end

############################ output construction ############################
# The labels left uncontracted after forming tensor-subset `S`: those touching
# `S` that are outputs or still appear on a tensor outside `S`.
function _open_labels(S::T, tensormask::Vector{T}, iymask::T, idlabels::Vector{L}) where {T, L}
    res = L[]
    @inbounds for l in 1:length(tensormask)
        if !_isemptyset(_intersect(tensormask[l], S)) && (_hasbit(iymask, l) || !_issubset(tensormask[l], S))
            push!(res, idlabels[l])
        end
    end
    return res
end

############################ scope validation ############################
function _check_exhaustive_scope(ixs::AbstractVector{<:AbstractVector}, iy)
    iyset = Set(iy)
    counts = Dict{eltype(iy), Int}()
    for ix in ixs
        seen = Set{eltype(iy)}()
        for l in ix
            if l in seen
                throw(ArgumentError("ExhaustiveSearch does not support partial traces: index $l appears more than once within a single tensor. Simplify it away or use another optimizer."))
            end
            push!(seen, l)
            counts[l] = get(counts, l, 0) + 1
        end
    end
    for (l, c) in counts
        if c == 1 && !(l in iyset)
            throw(ArgumentError("ExhaustiveSearch does not support a dangling index that is summed over: index $l appears once but is not an output. Simplify the network first or use another optimizer."))
        end
    end
    return nothing
end

############################ driver ############################
"""
    optimize_exhaustive(code::EinCode, size_dict; verbose=false) -> NestedEinsum

Optimize the contraction order of `code` exactly, minimizing the total
contraction cost (time complexity). See [`ExhaustiveSearch`](@ref).
"""
function optimize_exhaustive(code::EinCode{L}, size_dict::Dict{L}; verbose::Bool = false) where {L}
    ixs, iy = getixsv(code), getiyv(code)
    numtensors = length(ixs)
    numtensors <= 2 && return NestedEinsum(NestedEinsum{L}.(1:numtensors), EinCode(ixs, iy))
    _check_exhaustive_scope(ixs, iy)

    # intern labels to 1..numlabels and pick a bit type covering both the tensor
    # and label universes, then dispatch into a type-stable barrier on that type.
    idlabels = unique!(reduce(vcat, ixs))
    T = _choose_settype(max(numtensors, length(idlabels)))
    return _optimize_exhaustive(T, ixs, iy, idlabels, size_dict, verbose)
end

function _optimize_exhaustive(::Type{T}, ixs::Vector{Vector{L}}, iy::Vector{L}, idlabels::Vector{L}, size_dict, verbose::Bool) where {T, L}
    numtensors = length(ixs)
    numlabels = length(idlabels)
    lab2id = Dict(idlabels[i] => i for i in 1:numlabels)
    allcosts = Float64[Float64(get(size_dict, l, 1)) for l in idlabels]

    # per-tensor open-label set (= its full leg set; scope guarantees no dangling
    # or repeated legs), per-label occurrence set over tensors, and the
    # tensor-label incidence matrix used to find connected components
    singleopen = Vector{T}(undef, numtensors)
    tensormask = T[_zero_set(T, numtensors) for _ in 1:numlabels]
    incidence = falses(numtensors, numlabels)
    for i in 1:numtensors
        ids = [lab2id[l] for l in ixs[i]]
        singleopen[i] = _storeset(T, ids, numlabels)
        for id in ids
            tensormask[id] = _addbit!(tensormask[id], i)
            incidence[i, id] = true
        end
    end
    iymask = _storeset(T, [lab2id[l] for l in iy if haskey(lab2id, l)], numlabels)

    # connected components (tensors adjacent iff they share a label)
    componentlist = _connected_components(incidence, collect(1:numtensors))

    # solve each component, seeding its cost cap with a per-component greedy run
    ntrees = Vector{Any}(undef, length(componentlist))
    opens = Vector{T}(undef, length(componentlist))
    costs = Vector{Float64}(undef, length(componentlist))
    for (c, component) in enumerate(componentlist)
        initialcost = if length(component) == 1
            0.0
        else
            sub_ixs = [ixs[i] for i in component]
            complabels = Set(l for ix in sub_ixs for l in ix)
            sub_iy = L[l for l in iy if l in complabels]
            greedy = optimize_greedy(EinCode(sub_ixs, sub_iy), size_dict; α = 0.0, temperature = 0.0)
            Float64(flop(greedy, size_dict))
        end
        ntrees[c], opens[c], costs[c] = _solve_component(
            T, component, singleopen, tensormask, iymask, numtensors, allcosts, initialcost; verbose
        )
    end

    # combine disconnected components cheapest-first via outer products
    p = sortperm(costs)
    tree = ntrees[p[1]]
    open = opens[p[1]]
    for c in 2:length(componentlist)
        tree = Any[tree, ntrees[p[c]]]
        open = _union(open, opens[p[c]])
    end

    # build the NestedEinsum directly from the search state
    fullset = _storeset(T, 1:numtensors, numtensors)
    _, ne, _ = _build_nested(tree, ixs, iy, fullset, tensormask, iymask, idlabels, numtensors, T)
    return ne
end

function optimize_exhaustive(::NestedEinsum, size_dict; verbose::Bool = false)
    throw(ArgumentError("ExhaustiveSearch currently supports flat `EinCode` input only. Flatten the `NestedEinsum` (e.g. with `flatten`) before optimizing, or use another optimizer."))
end

# recursively turn the nested `Any[left,right]`/Int tree into a `NestedEinsum`
function _build_nested(node, ixs, iy, fullset::T, tensormask, iymask, idlabels::Vector{L}, numtensors, ::Type{T}) where {T, L}
    if node isa Integer
        i = Int(node)
        return (_storeset(T, [i], numtensors), NestedEinsum{L}(i), ixs[i])
    end
    s1, ne1, legs1 = _build_nested(node[1], ixs, iy, fullset, tensormask, iymask, idlabels, numtensors, T)
    s2, ne2, legs2 = _build_nested(node[2], ixs, iy, fullset, tensormask, iymask, idlabels, numtensors, T)
    S = _union(s1, s2)
    out = (S == fullset) ? iy : _open_labels(S, tensormask, iymask, idlabels)
    ne = NestedEinsum([ne1, ne2], EinCode([legs1, legs2], out))
    return (S, ne, out)
end

############################ small set utilities ############################
function _choose_settype(nbits::Int)
    nbits <= 32 && return UInt32
    nbits <= 64 && return UInt64
    (nbits <= 128 && !(Int === Int32 && Sys.iswindows())) && return UInt128
    return BitVector
end
_zero_set(::Type{T}, _maxint) where {T <: Unsigned} = zero(T)
_zero_set(::Type{BitVector}, maxint) = falses(maxint)
_addbit!(s::T, i) where {T <: Unsigned} = s | (one(T) << (i - 1))
function _addbit!(s::BitVector, i)
    s[i] = true
    return s
end
_setsize(s::Unsigned) = count_ones(s)
_setsize(s::BitVector) = count(s)
_singleton_index(s::Unsigned) = trailing_zeros(s) + 1
_singleton_index(s::BitVector) = findfirst(s)::Int
