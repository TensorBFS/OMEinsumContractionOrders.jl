################### expression tree ###################
# `ExprInfo` stores the node information.
# * `out_dims` is the output dimensions of this tree/subtree.
# * `tensorid` specifies the tensor index for leaf nodes. It is `-1` is for non-leaf node.
struct ExprInfo
    out_dims::Vector{Int}
    tensorid::Int
end
ExprInfo(out_dims::Vector{Int}) = ExprInfo(out_dims, -1)

# `ExprTree` is the expression tree for tensor contraction (or contraction tree), it is a binary tree (including leaf nodes without siblings).
# `left` and `right` are left and right branches, they are either both specified (non-leaf) or both unspecified (leaf), see [`isleaf`](@ref) function.
# `ExprTree()` for constructing a leaf node,
# `ExprTree(left, right, info)` for constructing a non-leaf node.
mutable struct ExprTree
    left::ExprTree
    right::ExprTree
    info::ExprInfo
    ExprTree(info) = (res = new(); res.info=info; res)
    ExprTree(left, right, info) = new(left, right, info)
end
function print_expr(io::IO, expr::ExprTree, level=0)
    isleaf(expr) && return print(io, " "^(2*level), labels(expr), " ($(expr.info.tensorid))")
    print(io, " "^(2*level), "(\n")
    print_expr(io, expr.left, level+1)
    print("\n")
    print_expr(io, expr.right, level+1)
    print("\n")
    print(io, " "^(2*level), ") := ", labels(expr))
end
# if `expr` is a leaf, it should have `left` and `right` fields both unspecified.
isleaf(expr::ExprTree) = !isdefined(expr, :left)
Base.show(io::IO, expr::ExprTree) = print_expr(io, expr, 0)
Base.show(io::IO, ::MIME"text/plain", expr::ExprTree) = show(io, expr)
siblings(t::ExprTree) = isleaf(t) ? ExprTree[] : ExprTree[t.left, t.right]
Base.copy(t::ExprTree) = isleaf(t) ? ExprTree(t.info) : ExprTree(copy(t.left), copy(t.right), copy(t.info))
Base.copy(info::ExprInfo) = ExprInfo(copy(info.out_dims), info.tensorid)
# output tensor labels
labels(t::ExprTree) = t.info.out_dims
# find the maximum label recursively, this is a helper function for converting an expression tree back to einsum.
maxlabel(t::ExprTree) = isleaf(t) ? maximum(isempty(labels(t)) ? 0 : labels(t)) : max(isempty(labels(t)) ? 0 : maximum(labels(t)), maxlabel(t.left), maxlabel(t.right))
# comparison between `ExprTree`s, mainly for testing
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims) && t1.tensorid == t2.tensorid
function _equal(t1::ExprTree, t2::ExprTree)
    isleaf(t1) != isleaf(t2) && return false
    isleaf(t1) ? t1.info == t2.info : _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
end
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)

############# Slicer ######################
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

function Base.push!(slicer::Slicer, best)
    @assert length(slicer) < slicer.max_size
    @assert !haskey(slicer.legs, best)
    slicer.legs[best] = slicer.log2_sizes[best]  # add best to legs
    slicer.log2_sizes[best] = 0.0
    return slicer
end

# convert the slicer to a vector of sliced labels
function get_slices(s::Slicer, inverse_map::Dict{Int,LT}) where LT
    # we want to keep the order of input fixed slices!
    LT[[inverse_map[l] for l in s.fixed_slices]..., [inverse_map[l] for (l, sz) in s.legs if l ∉ s.fixed_slices]...]
end

############# random expression tree ###############
function random_exprtree(code::EinCode)
    ixs, iy = getixsv(code), getiyv(code)
    labels = _label_dict(ixs, iy)
    return random_exprtree([Int[labels[l] for l in ix] for ix in ixs], Int[labels[l] for l in iy], length(labels))
end

function random_exprtree(ixs::Vector{Vector{Int}}, iy::Vector{Int}, nedge::Int)
    outercount = zeros(Int, nedge)
    allcount = zeros(Int, nedge)
    for l in iy
        outercount[l] += 1
        allcount[l] += 1
    end
    for ix in ixs
        for l in ix
            allcount[l] += 1
        end
    end
    _random_exprtree(ixs, collect(1:length(ixs)), outercount, allcount)
end
function _random_exprtree(ixs::Vector{Vector{Int}}, xindices, outercount::Vector{Int}, allcount::Vector{Int})
    n = length(ixs)
    if n == 1
        return ExprTree(ExprInfo(ixs[1], xindices[1]))
    end
    mask = rand(Bool, n)
    if all(mask) || !any(mask)  # prevent invalid partition
        i = rand(1:n)
        mask[i] = ~(mask[i])
    end
    info = ExprInfo(Int[i for i=1:length(outercount) if outercount[i]!=allcount[i] && outercount[i]!=0])
    outercount1, outercount2 = copy(outercount), copy(outercount)
    for i=1:n
        counter = mask[i] ? outercount2 : outercount1
        for l in ixs[i]
            counter[l] += 1
        end
    end
    return ExprTree(_random_exprtree(ixs[mask], xindices[mask], outercount1, allcount), _random_exprtree(ixs[(!).(mask)], xindices[(!).(mask)], outercount2, allcount), info)
end


##################### convert a contraction tree back to a nested einsum ####################
NestedEinsum(expr::ExprTree) = _nestedeinsum(expr, 1:maxlabel(expr))
NestedEinsum(expr::ExprTree, labelmap) = _nestedeinsum(expr, labelmap)
function _nestedeinsum(tree::ExprTree, lbs::Union{AbstractVector{LT}, Dict{Int,LT}}) where LT
    isleaf(tree) && return NestedEinsum{LT}(tree.info.tensorid)
    eins = EinCode([getindex.(Ref(lbs), labels(tree.left)), getindex.(Ref(lbs), labels(tree.right))], getindex.(Ref(lbs), labels(tree)))
    NestedEinsum([_nestedeinsum(tree.left, lbs), _nestedeinsum(tree.right, lbs)], eins)
end


##################### The main program ##############################
"""
    TreeSA{RT,IT,GM}
    TreeSA(; sc_target=20, βs=collect(0.01:0.05:15), ntrials=10, niters=50,
        sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_config=GreedyMethod(; nrepeat=1))

Optimize the einsum contraction pattern using the simulated annealing on tensor expression tree.

* `sc_target` is the target space complexity,
* `ntrials`, `βs` and `niters` are annealing parameters, doing `ntrials` indepedent annealings, each has inverse tempteratures specified by `βs`, in each temperature, do `niters` updates of the tree.
* `sc_weight` is the relative importance factor of space complexity in the loss compared with the time complexity.
* `rw_weight` is the relative importance factor of memory read and write in the loss compared with the time complexity.
* `initializer` specifies how to determine the initial configuration, it can be `:greedy` or `:random`. If it is using `:greedy` method to generate the initial configuration, it also uses two extra arguments `greedy_method` and `greedy_nrepeat`.
* `nslices` is the number of sliced legs, default is 0.
* `fixed_slices` is a vector of sliced legs, default is `[]`.

### References
* [Recursive Multi-Tensor Contraction for XEB Verification of Quantum Circuits](https://arxiv.org/abs/2108.05665)
"""
Base.@kwdef struct TreeSA{RT,IT,GM,LT} <: CodeOptimizer
    sc_target::RT = 20
    βs::IT = 0.01:0.05:15
    ntrials::Int = 10
    niters::Int = 50
    sc_weight::Float64 = 1.0
    rw_weight::Float64 = 0.2
    initializer::Symbol = :greedy
    nslices::Int = 0
    fixed_slices::Vector{LT} = []
    # configure greedy method
    greedy_config::GM = GreedyMethod(nrepeat=1)
end

# this is the main function
"""
    optimize_tree(code, size_dict; sc_target=20, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_method=MinSpaceOut(), fixed_slices=[])

Optimize the einsum contraction pattern specified by `code`, and edge sizes specified by `size_dict`.
Check the docstring of [`TreeSA`](@ref) for detailed explaination of other input arguments.
"""
function optimize_tree(code::AbstractEinsum, size_dict; nslices::Int=0, sc_target=20, βs=0.1:0.1:10, ntrials=20, niters=100, sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_method=GreedyMethod(nrepeat = 1), fixed_slices=[])
    LT = labeltype(code)
    if nslices < length(fixed_slices)
        @warn("Number of slices: $(nslices) is smaller than the number of fixed slices, setting it to: $(length(fixed_slices)).")
        nslices = length(fixed_slices)
    end
    # get input labels (`getixsv`) and output labels (`getiyv`) in the einsum code.
    ixs, iy = getixsv(code), getiyv(code)
    ninputs = length(ixs)  # number of input tensors
    if ninputs <= 2  # number of input tensors ≤ 2, can not be optimized
        return SlicedEinsum(LT[], NestedEinsum(NestedEinsum{LT}.(1:ninputs), EinCode(ixs, iy)))
    end
    ###### Stage 1: preprocessing ######
    labels = _label_dict(ixs, iy)  # map labels to integers
    inverse_map = Dict([v=>k for (k,v) in labels])  # the inverse transformation, map integers to labels
    log2_sizes = [log2.(size_dict[inverse_map[i]]) for i=1:length(labels)]   # use `log2` sizes in computing time
    if ntrials <= 0  # no optimization at all, then 1). initialize an expression tree and 2). convert back to nested einsum.
        best_tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method)
        return SlicedEinsum(LT[], NestedEinsum(best_tree, inverse_map))
    end
    ###### Stage 2: computing ######
    # create vectors to store optimized 1). expression tree, 2). time complexities, 3). space complexities, 4). read-write complexities and 5). slicing information.
    trees, tcs, scs, rws, slicers = Vector{ExprTree}(undef, ntrials), zeros(ntrials), zeros(ntrials), zeros(ntrials), Vector{Slicer}(undef, ntrials)
    @threads for t = 1:ntrials  # multi-threading on different trials, use `JULIA_NUM_THREADS=5 julia xxx.jl` for setting number of threads.
        # 1). random/greedy initialize a contraction tree.
        tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method)
        # 2). optimize the `tree` and `slicer` in a inplace manner.
        slicer = Slicer(log2_sizes, nslices, Int[labels[l] for l in fixed_slices])
        optimize_tree_sa!(tree, log2_sizes, slicer; sc_target=sc_target, βs=βs, niters=niters, sc_weight=sc_weight, rw_weight=rw_weight)
        # 3). evaluate time-space-readwrite complexities.
        tc, sc, rw = tree_timespace_complexity(tree, slicer.log2_sizes)
        @debug "trial $t, time complexity = $tc, space complexity = $sc, read-write complexity = $rw."
        trees[t], tcs[t], scs[t], rws[t], slicers[t] = tree, tc, sc, rw, slicer
    end
    ###### Stage 3: postprocessing ######
    # compare and choose the best solution
    best_tree, best_tc, best_sc, best_rw, best_slicer = first(trees), first(tcs), first(scs), first(rws), first(slicers)
    for i=2:ntrials
        if scs[i] < best_sc || (scs[i] == best_sc && exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(best_tc) + rw_weight * exp2(rws[i]))
            best_tree, best_tc, best_sc, best_rw, best_slicer = trees[i], tcs[i], scs[i], rws[i], slicers[i]
        end
    end
    @debug "best space complexities = $best_tc, time complexity = $best_sc, read-write complexity $best_rw."
    if best_sc > sc_target
        @warn "target space complexity not found, got: $best_sc, with time complexity $best_tc, read-write complexity $best_rw."
    end
    # returns a sliced einsum we need to map the sliced dimensions back from integers to labels.
    return SlicedEinsum(get_slices(best_slicer, inverse_map), NestedEinsum(best_tree, inverse_map))
end

# initialize a contraction tree
function _initializetree(code, size_dict, method; greedy_method)
    if method == :greedy
        labels = _label_dict(code)  # label to int
        return _exprtree(optimize_greedy(code, size_dict; α = greedy_method.α, temperature = greedy_method.temperature, nrepeat=greedy_method.nrepeat), labels)
    elseif method == :random
        return random_exprtree(code)
    elseif method == :specified
        labels = _label_dict(code)  # label to int
        return _exprtree(code, labels)
    else
        throw(ArgumentError("intializier `$method` is not defined!"))
    end
end

# use simulated annealing to optimize a contraction tree
function optimize_tree_sa!(tree::ExprTree, log2_sizes, slicer::Slicer; βs, niters, sc_target, sc_weight, rw_weight)
    @assert rw_weight >= 0
    @assert sc_weight >= 0
    log2rw_weight = log2(rw_weight)
    for β in βs
        @debug begin
            tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
            "β = $β, tc = $tc, sc = $sc, rw = $rw"
        end
        ###### Stage 1: add one slice at each temperature  ######
        if slicer.max_size > length(slicer.fixed_slices)  # `max_size` specifies the maximum number of sliced dimensions.
            # 1). find legs that reduce the dimension the most
            scs, lbs = Float64[], Vector{Int}[]
            # space complexities and labels of all intermediate tensors
            tensor_sizes!(tree, slicer.log2_sizes, scs, lbs)
            # the set of (intermediate) tensor labels that producing maximum space complexity
            best_labels = _best_labels(scs, lbs)

            # 2). slice the best not sliced label (it must appear in largest tensors)
            best_not_sliced_labels = filter(x->!haskey(slicer.legs, x), best_labels)
            if !isempty(best_not_sliced_labels)
                #best_not_sliced_label = rand(best_not_sliced_labels)  # random or best
                best_not_sliced_label = best_not_sliced_labels[findmax(l->count(==(l), best_not_sliced_labels), best_not_sliced_labels)[2]]
                if length(slicer) < slicer.max_size  # if has not reached maximum number of slices, add one slice
                    push!(slicer, best_not_sliced_label)
                else                                 # otherwise replace one slice
                    legs = [l for l in keys(slicer.legs) if l ∉ slicer.fixed_slices]  # only slice over not fixed legs
                    score = [count(==(l), best_labels) for l in legs]
                    replace!(slicer, legs[argmin(score)]=>best_not_sliced_label)
                end
            end
            @debug begin
                tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
                "after slicing: β = $β, tc = $tc, sc = $sc, rw = $rw"
            end
        end
        ###### Stage 2: sweep and optimize the contraction tree for `niters` times  ######
        for _ = 1:niters
            optimize_subtree!(tree, β, slicer.log2_sizes, sc_target, sc_weight, log2rw_weight)  # single sweep
        end
    end
    return tree, slicer
end

# here "best" means giving maximum space complexity
function _best_labels(scs, lbs)
    max_sc = maximum(scs)
    return vcat(lbs[scs .> max_sc-0.99]...)
end

# find tensor sizes and their corresponding labels of all intermediate tensors
function tensor_sizes!(tree::ExprTree, log2_sizes, scs, lbs)
    sc = isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree))
    push!(scs, sc)
    push!(lbs, labels(tree))
    isleaf(tree) && return
    tensor_sizes!(tree.left, log2_sizes, scs, lbs)
    tensor_sizes!(tree.right, log2_sizes, scs, lbs)
end

# the time-space-readwrite complexity of a contraction tree
function tree_timespace_complexity(tree::ExprTree, log2_sizes)
    isleaf(tree) && return (-Inf, isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree)), -Inf)
    tcl, scl, rwl = tree_timespace_complexity(tree.left, log2_sizes)
    tcr, scr, rwr = tree_timespace_complexity(tree.right, log2_sizes)
    tc, sc, rw = tcscrw(labels(tree.left), labels(tree.right), labels(tree), log2_sizes, true)
    return (fast_log2sumexp2(tc, tcl, tcr), max(sc, scl, scr), fast_log2sumexp2(rw, rwl, rwr))
end

# returns time complexity, space complexity and read-write complexity (0 if `compute_rw` is false)
# `ix1` and `ix2` are vectors of labels for the first and second input tensors.
# `iy` is a vector of labels for the output tensors.
# `log2_sizes` is the log2 size of labels (note labels are integers, we do not need dict to index label sizes).\
@inline function tcscrw(ix1, ix2, iy, log2_sizes::Vector{T}, compute_rw) where T
    l1, l2, l3 = ix1, ix2, iy
    sc1 = (!compute_rw || isempty(l1)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l1)
    sc2 = (!compute_rw || isempty(l2)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l2)
    sc = isempty(l3) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l3)
    tc = sc
    # Note: assuming labels in `l1` being unique
    @inbounds for l in l1
        if l ∈ l2 && l ∉ l3
            tc += log2_sizes[l]
        end
    end
    rw = compute_rw ? fast_log2sumexp2(sc, sc1, sc2) : 0.0
    return tc, sc, rw
end

# optimize a contraction tree recursively
function optimize_subtree!(tree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
    # find appliable local rules, at most 4 rules can be applied.
    # Sometimes, not all rules are applicable because either left or right sibling do not have siblings.
    rst = ruleset(tree)
    if !isempty(rst)
        # propose a random update rule, TODO: can we have a better selector?
        rule = rand(rst)
        optimize_rw = log2rw_weight != -Inf
        # difference in time, space and read-write complexity if the selected rule is applied
        tc0, tc1, dsc, rw0, rw1, subout = tcsc_diff(tree, rule, log2_sizes, optimize_rw)
        dtc = optimize_rw ? fast_log2sumexp2(tc1, log2rw_weight + rw1) - fast_log2sumexp2(tc0, log2rw_weight + rw0) : tc1 - tc0
        sc = _sc(tree, rule, log2_sizes)  # current space complexity

        # update the loss function
        dE = (max(sc, sc+dsc) > sc_target ? sc_weight : 0) * dsc + dtc
        if rand() < exp(-β*dE)  # ACCEPT
            update_tree!(tree, rule, subout)
        end
        for subtree in siblings(tree)  # RECURSE
            optimize_subtree!(subtree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
        end
    end
end
# if rule ∈ [1, 2], left sibling will be updated, otherwise, right sibling will be updated.
# we need to compute space complexity for current node and the updated sibling, and return the larger one.
_sc(tree, rule, log2_sizes) = max(__sc(tree, log2_sizes), __sc((rule == 1 || rule == 2) ? tree.left : tree.right, log2_sizes))
__sc(tree, log2_sizes) = length(labels(tree))==0 ? 0.0 : sum(l->log2_sizes[l], labels(tree)) # space complexity of current node

@inline function ruleset(tree::ExprTree)
    if isleaf(tree) || (isleaf(tree.left) && isleaf(tree.right))
        return 1:0
    elseif isleaf(tree.right)
        return 1:2
    elseif isleaf(tree.left)
        return 3:4
    else
        return 1:4
    end
end

function tcsc_diff(tree::ExprTree, rule, log2_sizes, optimize_rw)
    if rule == 1 # (a,b), c -> (a,c),b
        return abcacb(labels(tree.left.left), labels(tree.left.right), labels(tree.left), labels(tree.right), labels(tree), log2_sizes, optimize_rw)
    elseif rule == 2 # (a,b), c -> (c,b),a
        return abcacb(labels(tree.left.right), labels(tree.left.left), labels(tree.left), labels(tree.right), labels(tree), log2_sizes, optimize_rw)
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        return abcacb(labels(tree.right.right), labels(tree.right.left), labels(tree.right), labels(tree.left), labels(tree), log2_sizes, optimize_rw)
    else  # a,(b,c) -> c,(b,a)
        return abcacb(labels(tree.right.left), labels(tree.right.right), labels(tree.right), labels(tree.left), labels(tree), log2_sizes, optimize_rw)
    end
end

# compute the time complexity, space complexity and read-write complexity information for the contraction update rule "((a,b),c) -> ((a,c),b)"
function abcacb(a, b, ab, c, d, log2_sizes, optimize_rw)
    tc0, sc0, rw0 = _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
    ac = Int[]  # labels for contraction result of (a, c)
    for l in a
        if l ∈ b || l ∈ d  # suppose no repeated indices
            push!(ac, l)
        end
    end
    for l in c
        if l ∉ a && (l ∈ b || l ∈ d)  # suppose no repeated indices
            push!(ac, l)
        end
    end
    tc1, sc1, rw1 = _tcsc_merge(a, c, ac, b, d, log2_sizes, optimize_rw)
    return tc0, tc1, sc1-sc0, rw0, rw1, ac
end

# compute complexity for a two-step contraction: (a, b) -> ab, (ab, c) -> d
function _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
    tcl, scl, rwl = tcscrw(a, b, ab, log2_sizes, optimize_rw)  # this is correct
    tcr, scr, rwr = tcscrw(ab, c, d, log2_sizes, optimize_rw)
    fast_log2sumexp2(tcl, tcr), max(scl, scr), (optimize_rw ? fast_log2sumexp2(rwl, rwr) : 0.0)
end

# apply the update rule
function update_tree!(tree::ExprTree, rule::Int, subout)
    if rule == 1 # (a,b), c -> (a,c),b
        b, c = tree.left.right, tree.right
        tree.left.right = c
        tree.right = b
        tree.left.info = ExprInfo(subout)
    elseif rule == 2 # (a,b), c -> (c,b),a
        a, c = tree.left.left, tree.right
        tree.left.left = c
        tree.right = a
        tree.left.info = ExprInfo(subout)
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        a, b = tree.left, tree.right.left
        tree.left = b
        tree.right.left = a
        tree.right.info = ExprInfo(subout)
    else  # a,(b,c) -> c,(b,a)
        a, c = tree.left, tree.right.right
        tree.left = c
        tree.right.right = a
        tree.right.info = ExprInfo(subout)
    end
    return tree
end

# map labels to integers.
_label_dict(code) = _label_dict(getixsv(code), getiyv(code))
function _label_dict(ixsv::AbstractVector{<:AbstractVector{LT}}, iyv::AbstractVector{LT}) where LT
    v = unique(vcat(ixsv..., iyv))
    labels = Dict{LT,Int}(zip(v, 1:length(v)))
    return labels
end

# construct the contraction tree recursively from a nested einsum.
function ExprTree(code::NestedEinsum)
    _exprtree(code, _label_dict(code))
end
function _exprtree(code::NestedEinsum, labels)
    @assert length(code.args) == 2 "einsum contraction not in the binary form, got number of arguments: $(length(code.args))"
    ExprTree(map(enumerate(code.args)) do (i,arg)
        if isleaf(arg)  # leaf nodes
            ExprTree(ExprInfo(getindex.(Ref(labels), getixsv(code.eins)[i]), arg.tensorindex))
        else
            res = _exprtree(arg, labels)
        end
    end..., ExprInfo(Int[labels[i] for i=getiyv(code.eins)]))
end

@inline function fast_log2sumexp2(a, b)
    mm, ms = minmax(a, b)
    return log2(exp2(mm - ms) + 1) + ms
end

@inline function fast_log2sumexp2(a, b, c)
    if a > b
        if a > c
            m1, m2, ms = b, c, a
        else
            m1, m2, ms = a, b, c
        end
    else
        if b > c
            m1, m2, ms = c, a, b
        else
            m1, m2, ms = b, a, c
        end
    end
    return Base.FastMath.log2(Base.FastMath.exp2(m1 - ms) + Base.FastMath.exp2(m2 - ms) + 1) + ms
end

